#include "chrono_planet/lod/TileMeshBuilder.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>

#include "chrono_planet/core/BenchProfiler.h"
#include "chrono_planet/core/FilterChain.h"
#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/core/Planet.h"
#include "chrono_planet/core/SphereMath.h"
#include "chrono_planet/lod/SphericalCoordinates.h"

using namespace qtplanet;

int MeshTopology::divisions(int level) { return std::min(1 << (level + 1), kMaxDivisions); }

double MeshTopology::vertexSpacingDeg(int level, double rootDeg) {
    return std::ldexp(rootDeg, -std::max(0, level)) / divisions(std::max(0, level));
}

size_t MeshTopology::vertexCount(int level) {
    const size_t n = gridSide(level);
    return n * n + 4 * n;
}

namespace {

// Bake sampling stays this far from the poles, where the lon/lat lattice degenerates.
constexpr double kBakeLatLimitDeg = 89.99;

// Grid triangles followed by one skirt strip per border, for a grid of n vertices per side.
std::vector<unsigned int> buildIndices(int n) {
    const int divisions = n - 1;
    std::vector<unsigned int> out;
    out.reserve(static_cast<size_t>(divisions) * divisions * 6 + 4 * (n - 1) * 6);
    for (int j = 0; j < divisions; ++j) {
        for (int i = 0; i < divisions; ++i) {
            const unsigned int tl = j * n + i, tr = tl + 1, bl = (j + 1) * n + i, br = bl + 1;
            out.insert(out.end(), {tl, bl, tr, tr, bl, br});
        }
    }
    // Skirts, one strip per border, hanging off the skirt rows after the grid.
    unsigned int base = static_cast<unsigned int>(n) * n;
    auto strip = [&](auto&& topIndex) {
        for (int k = 0; k + 1 < n; ++k) {
            const unsigned int t0 = topIndex(k), t1 = topIndex(k + 1);
            const unsigned int b0 = base + k, b1 = base + k + 1;
            out.insert(out.end(), {t0, b0, t1, t1, b0, b1});
        }
        base += n;
    };
    strip([n](int k) { return static_cast<unsigned int>(k); });                 // south row
    strip([n](int k) { return static_cast<unsigned int>((n - 1) * n + k); });   // north row
    strip([n](int k) { return static_cast<unsigned int>(k * n); });             // west column
    strip([n](int k) { return static_cast<unsigned int>(k * n + (n - 1)); });   // east column
    return out;
}

}   // namespace

const std::vector<unsigned int>& MeshTopology::indices(int level) {
    // One table per distinct topology, all built on first use so any thread may read them afterwards.
    static const auto tables = [] {
        std::array<std::vector<unsigned int>, kMaxDistinctLevel + 1> t;
        for (int l = 0; l <= kMaxDistinctLevel; ++l) {
            t[l] = buildIndices(gridSide(l));
        }
        return t;
    }();
    return tables[std::clamp(level, 0, kMaxDistinctLevel)];
}

template <typename C>
Mesh TileMeshBuilder<C>::build(const Boundary& bounds, int level) const {
    BENCH_SCOPE("generate");
    return filterChain(BuildRequest{bounds, level, loader_.get()})
        .then(sampleTerrain)
        .then(measureGeometry)
        .then(packMesh)
        .then(bakeSurface)
        .take();
}

template <typename C>
typename TileMeshBuilder<C>::SampledTile TileMeshBuilder<C>::sampleTerrain(BuildRequest request) {
    const int divisions = MeshTopology::divisions(request.level);
    auto fine = sampleFineGrid(request.bounds, request.level, divisions, request.loader);
    auto parent = sampleParentLattice(request.bounds, request.level, divisions, request.loader);
    return {request, std::move(fine), std::move(parent)};
}

template <typename C>
typename TileMeshBuilder<C>::PreparedTile TileMeshBuilder<C>::measureGeometry(SampledTile sampled) {
    const GeometrySummary geometry = summarizeGeometry(sampled.fine.positionsM);
    const double skirtDepthM = skirtDepthForLevel(sampled.request.level);
    return {std::move(sampled), geometry, skirtDepthM};
}

template <typename C>
typename TileMeshBuilder<C>::PackedTile TileMeshBuilder<C>::packMesh(PreparedTile prepared) {
    const auto& request = prepared.sampled.request;
    const auto& bounds = request.bounds;
    const auto& fine = prepared.sampled.fine;
    const auto& geometry = prepared.geometry;
    Mesh mesh;
    mesh.level = request.level;
    mesh.centerX = geometry.centerM.x;
    mesh.centerY = geometry.centerM.y;
    mesh.centerZ = geometry.centerM.z;
    mesh.minElevation = geometry.minElevationM;
    mesh.maxElevation = geometry.maxElevationM;
    mesh.radius = boundingRadiusM(fine.positionsM, geometry.centerM, prepared.skirtDepthM);
    mesh.uvMin[0] = (bounds.centerLonDeg - bounds.halfWidthDeg + 180.0) / 360.0;
    mesh.uvMin[1] = (bounds.centerLatDeg - bounds.halfHeightDeg + 90.0) / 180.0;
    mesh.uvSize[0] = 2.0 * bounds.halfWidthDeg / 360.0;
    mesh.uvSize[1] = 2.0 * bounds.halfHeightDeg / 180.0;
    mesh.vertexData = packVertices(fine, prepared.sampled.parent, geometry.centerM, prepared.skirtDepthM);
    return {request, std::move(mesh)};
}

template <typename C>
Mesh TileMeshBuilder<C>::bakeSurface(PackedTile packed) {
    const auto& request = packed.request;
    packed.mesh.bake = bakeShading(request.bounds, request.level, request.loader);
    return std::move(packed.mesh);
}

template <typename C>
double TileMeshBuilder<C>::skirtDepthForLevel(int level) {
    // Walls stay deeper than the sag a coarser neighbor's T-junctions open.
    static constexpr double depthByLevel[] = {2000.0, 1000.0, 500.0, 250.0, 150.0, 80.0, 70.0, 50.0, 30.0, 20.0};
    constexpr int count = sizeof(depthByLevel) / sizeof(depthByLevel[0]);
    static const bool noSkirt = std::getenv("QT_NOSKIRT") != nullptr;
    if (noSkirt) {
        return 0.0;
    }
    return depthByLevel[std::clamp(level, 0, count - 1)];
}

template <typename C>
typename TileMeshBuilder<C>::FineGrid TileMeshBuilder<C>::sampleFineGrid(
    const Boundary& bounds, int level, int divisions, const MultiGeoTIFFManager* loader) {
    FineGrid grid;
    grid.side = divisions + 1;
    {
        BENCH_SCOPE("grid");
        grid.positionsM = CoordinateTraits<C>::cartesianGrid(bounds, divisions, loader, level);
        grid.uv.resize(static_cast<size_t>(grid.side) * grid.side * 2);
        for (int row = 0; row < grid.side; ++row) {
            for (int column = 0; column < grid.side; ++column) {
                const size_t vertex = static_cast<size_t>(row) * grid.side + column;
                grid.uv[vertex * 2] = static_cast<float>(static_cast<double>(column) / divisions);
                grid.uv[vertex * 2 + 1] = static_cast<float>(static_cast<double>(row) / divisions);
            }
        }
    }
    {
        BENCH_SCOPE("slopes");
        grid.slopes = slopesFor(grid.positionsM, grid.side);
    }
    return grid;
}

// The parent's lattice, one ring of cells wider than this tile, so every fine vertex has a parent cell.
template <typename C>
std::optional<typename TileMeshBuilder<C>::ParentLattice> TileMeshBuilder<C>::sampleParentLattice(
    const Boundary& bounds, int level, int divisions, const MultiGeoTIFFManager* loader) {
    BENCH_SCOPE("lattice");
    const int parentDivisions = level > 0 ? MeshTopology::divisions(level - 1) / 2 : 0;   // parent cells across this tile
    if (parentDivisions < 1) {
        return std::nullopt;
    }
    ParentLattice parent;
    parent.divisions = parentDivisions;
    parent.fineCellsPerCell = divisions / parentDivisions;   // fine cells per parent cell
    const double parentStepLonDeg = 2.0 * bounds.halfWidthDeg / parent.divisions;
    Boundary ring = bounds;
    ring.halfWidthDeg += parentStepLonDeg;
    ring.halfHeightDeg += 2.0 * bounds.halfHeightDeg / parent.divisions;
    parent.side = parent.divisions + 3;   // coarse grid side incl. ring
    parent.positionsM = CoordinateTraits<C>::cartesianGrid(ring, parent.divisions + 2, loader, level - 1);
    parent.slopes = slopesFor(parent.positionsM, parent.side);
    return parent;
}

template <typename C>
typename TileMeshBuilder<C>::GeometrySummary TileMeshBuilder<C>::summarizeGeometry(const std::vector<double>& positionsM) {
    const size_t gridCount = positionsM.size() / 3;
    double cx = 0, cy = 0, cz = 0, maxElev = -1e30, minElev = 1e30;
    for (size_t v = 0; v < gridCount; ++v) {
        cx += positionsM[v * 3];
        cy += positionsM[v * 3 + 1];
        cz += positionsM[v * 3 + 2];
        const double e = CoordinateTraits<C>::elevationOf(positionsM[v * 3], positionsM[v * 3 + 1], positionsM[v * 3 + 2]);
        maxElev = std::max(maxElev, e);
        minElev = std::min(minElev, e);
    }
    return {{cx / double(gridCount), cy / double(gridCount), cz / double(gridCount)}, minElev, maxElev};
}

template <typename C>
std::vector<float> TileMeshBuilder<C>::packVertices(const FineGrid& fine, const std::optional<ParentLattice>& parent,
                                                    const Vec3& centerM, double skirtDepthM) {
    const size_t gridCount = static_cast<size_t>(fine.side) * fine.side;
    std::vector<float> vertices((gridCount + 4 * fine.side) * Mesh::kFloatsPerVertex);
    {
        BENCH_SCOPE("vertices");
        for (int j = 0; j < fine.side; ++j) {
            for (int i = 0; i < fine.side; ++i) {
                const size_t v = static_cast<size_t>(j) * fine.side + i;
                const Vec3 positionM = Vec3::from(&fine.positionsM[v * 3]);
                writeVertex(vertices, v, {positionM - centerM, {fine.slopes[v * 2], fine.slopes[v * 2 + 1]}, {fine.uv[v * 2], fine.uv[v * 2 + 1]}, parentPositionAt(fine, parent, j, i) - centerM, parentSlopeAt(fine, parent, j, i)});
            }
        }
    }
    {
        BENCH_SCOPE("skirts");
        size_t dst = gridCount;
        const int n = fine.side;
        auto skirt = [&](auto&& gridIndexAt) {
            for (int k = 0; k < n; ++k) {
                const size_t src = gridIndexAt(k);
                const Vec3 skirtPositionM = (CoordinateTraits<C>::skirtVertex(Vec3::from(&fine.positionsM[src * 3]), skirtDepthM) - centerM);
                const Slope2 slope{fine.slopes[src * 2], fine.slopes[src * 2 + 1]};
                writeVertex(vertices, dst++, {skirtPositionM, slope, {fine.uv[src * 2], fine.uv[src * 2 + 1]}, skirtPositionM, slope});
            }
        };
        skirt([n](int k) { return static_cast<size_t>(k); });                 // south row
        skirt([n](int k) { return static_cast<size_t>((n - 1) * n + k); });   // north row
        skirt([n](int k) { return static_cast<size_t>(k * n); });             // west column
        skirt([n](int k) { return static_cast<size_t>(k * n + (n - 1)); });   // east column
    }
    return vertices;
}

template <typename C>
double TileMeshBuilder<C>::boundingRadiusM(const std::vector<double>& positionsM, const Vec3& centerM, double skirtDepthM) {
    double radiusSquaredM2 = 0.0;
    for (size_t vertex = 0; vertex < positionsM.size() / 3; ++vertex) {
        const double dx = positionsM[vertex * 3] - centerM.x;
        const double dy = positionsM[vertex * 3 + 1] - centerM.y;
        const double dz = positionsM[vertex * 3 + 2] - centerM.z;
        radiusSquaredM2 = std::max(radiusSquaredM2, dx * dx + dy * dy + dz * dz);
    }
    return std::sqrt(radiusSquaredM2) + skirtDepthM;
}

// Height above the sphere of the packed cartesian points, as floats. Rows go to the team on the bake-sized grids.
template <typename C>
std::vector<float> TileMeshBuilder<C>::elevationsM(const std::vector<double>& positionsM) {
    const size_t count = positionsM.size() / 3;
    std::vector<float> heightsM(count);
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(static) if (qtplanet::parallelGrid(count))
    for (size_t v = 0; v < count; ++v) {
        heightsM[v] = static_cast<float>(CoordinateTraits<C>::elevationOf(positionsM[v * 3], positionsM[v * 3 + 1], positionsM[v * 3 + 2]));
    }
    return heightsM;
}

// Unnormalized normal map over tile+margin at kBakeZoom and height map over tile+kHeightMargin.
template <typename C>
std::shared_ptr<const ShadingBake> TileMeshBuilder<C>::bakeShading(
    const Boundary& bounds, int level, const MultiGeoTIFFManager* loader) {
    BENCH_SCOPE("bake");
    if (level < Mesh::kBakeMinLevel) {
        return nullptr;
    }
    Boundary normalBakeBounds = bounds;
    normalBakeBounds.halfWidthDeg *= (1.0 + 2.0 * Mesh::kBakeMargin);
    normalBakeBounds.halfHeightDeg *= (1.0 + 2.0 * Mesh::kBakeMargin);
    const double bMinLat = normalBakeBounds.centerLatDeg - normalBakeBounds.halfHeightDeg, bMaxLat = normalBakeBounds.centerLatDeg + normalBakeBounds.halfHeightDeg;
    if (bMinLat <= -90.0 || bMaxLat >= 90.0) {
        return nullptr;
    }

    const int res = Mesh::bakeRes(level), samplesPerTexelAxis = Mesh::bakeSupersample(level);
    const int normalSampleSide = res * samplesPerTexelAxis;   // supersampled samples per side
    auto bake = std::make_shared<ShadingBake>();
    bake->res = res;
    bake->minLon = normalBakeBounds.centerLonDeg - normalBakeBounds.halfWidthDeg;
    bake->maxLon = normalBakeBounds.centerLonDeg + normalBakeBounds.halfWidthDeg;
    bake->minLat = bMinLat;
    bake->maxLat = bMaxLat;
    // Sample at (sub)texel centers. Texel i = mean of its ss x ss samples.
    const double normalStepLonDeg = 2.0 * normalBakeBounds.halfWidthDeg / normalSampleSide, normalStepLatDeg = 2.0 * normalBakeBounds.halfHeightDeg / normalSampleSide;
    Boundary normalSampleBounds = normalBakeBounds;
    normalSampleBounds.halfWidthDeg -= 0.5 * normalStepLonDeg;
    normalSampleBounds.halfHeightDeg -= 0.5 * normalStepLatDeg;
    std::vector<float> normalSampleHeightsM;
    {
        BENCH_SCOPE("bake_dem");
        const std::vector<double> normalSamplePositionsM = CoordinateTraits<C>::cartesianGrid(normalSampleBounds, normalSampleSide - 1, loader, Mesh::kBakeZoom);
        BENCH_SCOPE("dem_elev");
        normalSampleHeightsM = elevationsM(normalSamplePositionsM);
    }
    BENCH_SCOPE("bake_normals");
    const double dyM = metresPerDegLat(kRadiusM) * normalStepLatDeg;
    const size_t count = static_cast<size_t>(res) * res;
    bake->normals.resize(count * 3);
    // A sample row of unit normals at a time, then the ss x ss mean per texel.
#pragma omp parallel num_threads(qtplanet::parallelThreads()) if (qtplanet::parallelGrid(count))
    {
        std::vector<double> nx(normalSampleSide), ny(normalSampleSide), nz(normalSampleSide);
        std::vector<Vec3> acc(res);
#pragma omp for schedule(static)
        for (int ty = 0; ty < res; ++ty) {
            std::fill(acc.begin(), acc.end(), Vec3{});
            for (int sy = 0; sy < samplesPerTexelAxis; ++sy) {
                const int y = ty * samplesPerTexelAxis + sy;
                const double lat = bMinLat + (y + 0.5) * normalStepLatDeg;
                const double dxM = std::max(1.0, metresPerDegLon(kRadiusM, lat) * normalStepLonDeg);
                const float* rowM = &normalSampleHeightsM[static_cast<size_t>(std::max(y - 1, 0)) * normalSampleSide];
                const float* rowP = &normalSampleHeightsM[static_cast<size_t>(std::min(y + 1, normalSampleSide - 1)) * normalSampleSide];
                const float* row = &normalSampleHeightsM[static_cast<size_t>(y) * normalSampleSide];
                auto at = [&](int x, int xm, int xp) {
                    const float dhx = row[xp] - row[xm], dhy = rowP[x] - rowM[x];
                    const double gx = dhx / (2.0 * dxM), gy = dhy / (2.0 * dyM);
                    const double len = std::sqrt(gx * gx + gy * gy + 1.0);
                    nx[x] = -gx / len;
                    ny[x] = -gy / len;
                    nz[x] = 1.0 / len;
                };
                at(0, 0, 1);
                for (int x = 1; x < normalSampleSide - 1; ++x) {
                    at(x, x - 1, x + 1);
                }
                at(normalSampleSide - 1, normalSampleSide - 2, normalSampleSide - 1);
                for (int tx = 0; tx < res; ++tx) {
                    for (int sx = 0; sx < samplesPerTexelAxis; ++sx) {
                        const int x = tx * samplesPerTexelAxis + sx;
                        acc[tx] += Vec3{nx[x], ny[x], nz[x]};
                    }
                }
            }
            for (int tx = 0; tx < res; ++tx) {
                Vec3 nsum = acc[tx];
                nsum /= double(samplesPerTexelAxis * samplesPerTexelAxis);
                const std::array<unsigned short, 3> packed = packSnorm16(nsum);
                std::copy(packed.begin(), packed.end(), &bake->normals[(static_cast<size_t>(ty) * res + tx) * 3]);
            }
        }
    }

    {
        BENCH_SCOPE("bake_height");
        Boundary heightBakeBounds = bounds;
        heightBakeBounds.halfWidthDeg *= (1.0 + 2.0 * Mesh::kHeightMargin);
        const double hMinLat = std::max(bounds.centerLatDeg - bounds.halfHeightDeg * (1.0 + 2.0 * Mesh::kHeightMargin), -kBakeLatLimitDeg);
        const double hMaxLat = std::min(bounds.centerLatDeg + bounds.halfHeightDeg * (1.0 + 2.0 * Mesh::kHeightMargin), kBakeLatLimitDeg);
        heightBakeBounds.centerLatDeg = 0.5 * (hMinLat + hMaxLat);
        heightBakeBounds.halfHeightDeg = 0.5 * (hMaxLat - hMinLat);
        const int heightSide = Mesh::heightRes(level);
        const double heightStepLonDeg = 2.0 * heightBakeBounds.halfWidthDeg / heightSide, heightStepLatDeg = 2.0 * heightBakeBounds.halfHeightDeg / heightSide;
        Boundary heightSampleBounds = heightBakeBounds;   // texel centers
        heightSampleBounds.halfWidthDeg -= 0.5 * heightStepLonDeg;
        heightSampleBounds.halfHeightDeg -= 0.5 * heightStepLatDeg;
        const std::vector<double> heightSamplePositionsM = CoordinateTraits<C>::cartesianGrid(heightSampleBounds, heightSide - 1, loader, Mesh::kBakeZoom, normalStepLonDeg);
        bake->hRes = heightSide;
        bake->hMinLon = heightBakeBounds.centerLonDeg - heightBakeBounds.halfWidthDeg;
        bake->hMaxLon = heightBakeBounds.centerLonDeg + heightBakeBounds.halfWidthDeg;
        bake->hMinLat = hMinLat;
        bake->hMaxLat = hMaxLat;
        BENCH_SCOPE("height_elev");
        bake->heights = elevationsM(heightSamplePositionsM);
    }
    return bake;
}

// Area-weighted face normals accumulated on a (side x side) grid.
template <typename C>
std::vector<double> TileMeshBuilder<C>::smoothNormals(const std::vector<double>& positionsM, int side) {
    std::vector<double> out;
    out.assign(static_cast<size_t>(side) * side * 3, 0.0);
    for (int j = 0; j + 1 < side; ++j) {
        for (int i = 0; i + 1 < side; ++i) {
            const size_t tl = static_cast<size_t>(j) * side + i, tr = tl + 1, bl = tl + side, br = bl + 1;
            const size_t tris[2][3] = {{tl, bl, tr}, {tr, bl, br}};
            for (const auto& t : tris) {
                const Vec3 a = Vec3::from(&positionsM[t[0] * 3]);
                // j grows north, i grows east.
                const Vec3 faceNormal = cross(Vec3::from(&positionsM[t[2] * 3]) - a, Vec3::from(&positionsM[t[1] * 3]) - a);
                for (size_t k : t) {
                    (Vec3::from(&out[k * 3]) + faceNormal).store(&out[k * 3]);
                }
            }
        }
    }
    for (size_t v = 0; v < out.size() / 3; ++v) {
        normalize(Vec3::from(&out[v * 3])).store(&out[v * 3]);
    }
    return out;
}

// A normal expressed as east/north rise over its own up component.
template <typename C>
typename TileMeshBuilder<C>::Slope2 TileMeshBuilder<C>::slopeAt(const Vec3& positionM, const Vec3& normal) {
    const EnuFrame frame = enuAlong(positionM);
    const double upComponent = std::max(dot(normal, frame.up), 0.05);
    return {dot(normal, frame.east) / upComponent, dot(normal, frame.north) / upComponent};
}

// Per-vertex slope of the smoothed grid normal.
template <typename C>
std::vector<double> TileMeshBuilder<C>::slopesFor(const std::vector<double>& positionsM, int side) {
    std::vector<double> outSlope2;
    const auto normals = smoothNormals(positionsM, side);
    const size_t count = static_cast<size_t>(side) * side;
    outSlope2.resize(count * 2);
    for (size_t v = 0; v < count; ++v) {
        const Slope2 slope = slopeAt(Vec3::from(&positionsM[v * 3]), Vec3::from(&normals[v * 3]));
        outSlope2[v * 2] = slope.east;
        outSlope2[v * 2 + 1] = slope.north;
    }
    return outSlope2;
}

template <typename C>
Vec3 TileMeshBuilder<C>::parentPositionAt(const FineGrid& fine, const std::optional<ParentLattice>& parent, int row, int column) {
    if (!parent) {
        return Vec3::from(&fine.positionsM[(static_cast<size_t>(row) * fine.side + column) * 3]);
    }
    return {parentComponentAt(*parent, parent->positionsM, 3, row, column, 0),
            parentComponentAt(*parent, parent->positionsM, 3, row, column, 1),
            parentComponentAt(*parent, parent->positionsM, 3, row, column, 2)};
}

template <typename C>
typename TileMeshBuilder<C>::Slope2 TileMeshBuilder<C>::parentSlopeAt(const FineGrid& fine, const std::optional<ParentLattice>& parent, int row, int column) {
    if (!parent) {
        const size_t offset = (static_cast<size_t>(row) * fine.side + column) * 2;
        return {fine.slopes[offset], fine.slopes[offset + 1]};
    }
    return {parentComponentAt(*parent, parent->slopes, 2, row, column, 0),
            parentComponentAt(*parent, parent->slopes, 2, row, column, 1)};
}

// Packed parent samples include one extra cell around the tile.
template <typename C>
double TileMeshBuilder<C>::parentComponentAt(const ParentLattice& parent, const std::vector<double>& samples, int componentsPerSample,
                                             int row, int column, int component) {
    const int southRow = row / parent.fineCellsPerCell,
              westColumn = column / parent.fineCellsPerCell;   // parent cell
    const int northRow = std::min(southRow + 1, parent.divisions),
              eastColumn = std::min(westColumn + 1, parent.divisions);
    const double rowFraction = double(row - southRow * parent.fineCellsPerCell) / parent.fineCellsPerCell,
                 columnFraction = double(column - westColumn * parent.fineCellsPerCell) / parent.fineCellsPerCell;
    auto at = [&](int r, int c) {
        return samples[(static_cast<size_t>(r + 1) * parent.side + (c + 1)) * componentsPerSample + component];
    };
    const double southValue = at(southRow, westColumn) * (1.0 - columnFraction) + at(southRow, eastColumn) * columnFraction;
    const double northValue = at(northRow, westColumn) * (1.0 - columnFraction) + at(northRow, eastColumn) * columnFraction;
    return southValue * (1.0 - rowFraction) + northValue * rowFraction;
}

template <typename C>
void TileMeshBuilder<C>::writeVertex(std::vector<float>& vertices, size_t vertexIndex, const VertexSample& sample) {
    float* output = &vertices[vertexIndex * Mesh::kFloatsPerVertex];
    auto writePosition = [&output](const Vec3& positionM) {
        *output++ = static_cast<float>(positionM.x);
        *output++ = static_cast<float>(positionM.y);
        *output++ = static_cast<float>(positionM.z);
    };
    auto writeSlope = [&output](Slope2 slope) {
        *output++ = static_cast<float>(slope.east);
        *output++ = static_cast<float>(slope.north);
    };
    writePosition(sample.positionRelativeM);
    writeSlope(sample.slope);
    *output++ = sample.tileUV[0];
    *output++ = sample.tileUV[1];
    writePosition(sample.parentPositionRelativeM);
    writeSlope(sample.parentSlope);
}

template class TileMeshBuilder<Spherical>;
