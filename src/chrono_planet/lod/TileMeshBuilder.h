#ifndef QTPLANET_TILEMESHBUILDER_H
#define QTPLANET_TILEMESHBUILDER_H

#include "chrono_planet/ChApiPlanet.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/lod/CoordinateSystems.h"
#include "chrono_planet/lod/TileMetadata.h"

class MultiGeoTIFFManager;

// One tile's renderable geometry, plus its shading bake.
struct CH_PLANET_API Mesh {
    // Vertex layout is pos(3) slope(2) uv(2) parentPos(3) parentSlope(2).
    static constexpr int kFloatsPerVertex = 12;
    static constexpr int kBakeMinLevel = 2;                            // no shading bake below this level
    static constexpr double kBakeMargin = 1.0 / 16.0;                  // normal map filtering room past the edge
    static constexpr double kHeightMargin = 1.0;                       // height map reach past every edge, in tile widths
    static constexpr int kBakeZoom = 4;                                // DEM zoom the bakes sample at
    static int bakeRes(int level) { return level <= 2 ? 512 : 256; }   // normal map texels per side
    // Supersamples per texel per axis for the normal map.
    static int bakeSupersample(int level) { return level <= 2 ? 1 : 2; }
    // Height map texels per side. Must be 3 * 2^k so neighboring maps share one lattice and mip grouping.
    static int heightRes(int level) { return level <= 2 ? 768 : 384; }

    std::shared_ptr<const ShadingBake> bake;
    std::uint64_t id = ShadingBake::nextId();   // process-unique, in build order
    std::vector<float> vertexData;              // vertexCount(level) * kFloatsPerVertex, center-relative
    int level = 0;
    double uvMin[2] = {0, 0}, uvSize[2] = {1, 1};               // tile extent in global UV
    double centerX = 0, centerY = 0, centerZ = 0, radius = 0;   // bounding sphere, skirts included
    // Elevation range of the grid, meters above the sphere.
    double maxElevation = 0;
    double minElevation = 0;
};

// Grid dimensions and the shared index buffer per level.
struct CH_PLANET_API MeshTopology {
    static constexpr int kMaxDivisions = 64;      // cells per side cap
    static constexpr int kMaxDistinctLevel = 5;   // divisions() saturates here, deeper levels share its topology

    static int divisions(int level);   // cells per side, capped at kMaxDivisions
    static int gridSide(int level) { return divisions(level) + 1; }
    static double vertexSpacingDeg(int level, double rootDeg);    // degrees between grid vertices
    static size_t vertexCount(int level);                         // grid plus the four skirt rows
    static const std::vector<unsigned int>& indices(int level);   // grid + skirt triangles, built once, thread-safe
};

// Builds one tile's Mesh from the DEM stack. Instantiated for Spherical only.
template <typename CoordSystem>
class TileMeshBuilder {
public:
    using Boundary = typename CoordSystem::Boundary;

    explicit TileMeshBuilder(const std::shared_ptr<const MultiGeoTIFFManager>& loader) : loader_(loader) {}

    // Builds an independent mesh. No per-build state is retained; repeated calls may
    // use different bounds and levels. Allocates mesh/bake IDs and records profiling.
    [[nodiscard]] Mesh build(const Boundary& bounds, int level) const;

    // Skirt depth in meters at level. QT_NOSKIRT=1 disables skirts.
    static double skirtDepthForLevel(int level);

private:
    struct FineGrid {
        int side = 0;
        std::vector<double> positionsM, slopes;   // packed xyz positions and east/north slopes
        std::vector<float> uv;                    // packed tile-local UV pairs
    };
    struct ParentLattice {
        int side = 0, divisions = 0, fineCellsPerCell = 1;
        std::vector<double> positionsM, slopes;
    };
    struct GeometrySummary {
        qtplanet::Vec3 centerM;
        double minElevationM, maxElevationM;
    };

    // Stage types encode which data is available. All buffers move forward;
    // the DEM pointer is borrowed only for the synchronous build call.
    struct BuildRequest {
        Boundary bounds;
        int level;
        const MultiGeoTIFFManager* loader;
    };
    struct SampledTile {
        BuildRequest request;
        FineGrid fine;
        std::optional<ParentLattice> parent;
    };
    struct PreparedTile {
        SampledTile sampled;
        GeometrySummary geometry;
        double skirtDepthM;
    };
    struct PackedTile {
        BuildRequest request;
        Mesh mesh;
    };

    [[nodiscard]] static SampledTile sampleTerrain(BuildRequest request);
    [[nodiscard]] static PreparedTile measureGeometry(SampledTile sampled);
    [[nodiscard]] static PackedTile packMesh(PreparedTile prepared);
    [[nodiscard]] static Mesh bakeSurface(PackedTile packed);

    [[nodiscard]] static FineGrid sampleFineGrid(const Boundary& bounds, int level, int divisions, const MultiGeoTIFFManager* loader);
    [[nodiscard]] static std::optional<ParentLattice> sampleParentLattice(const Boundary& bounds, int level, int divisions, const MultiGeoTIFFManager* loader);
    [[nodiscard]] static GeometrySummary summarizeGeometry(const std::vector<double>& positionsM);
    [[nodiscard]] static double boundingRadiusM(const std::vector<double>& positionsM, const qtplanet::Vec3& centerM, double skirtDepthM);
    [[nodiscard]] static std::vector<float> packVertices(const FineGrid& fine, const std::optional<ParentLattice>& parent,
                                                         const qtplanet::Vec3& centerM, double skirtDepthM);
    // Returns no bake below the bake level or when the normal-map bounds reach a pole.
    [[nodiscard]] static std::shared_ptr<const ShadingBake> bakeShading(const Boundary& bounds, int level, const MultiGeoTIFFManager* loader);

    // The normal's east/north components divided by its up component.
    struct Slope2 {
        double east, north;
    };
    struct VertexSample {
        qtplanet::Vec3 positionRelativeM;
        Slope2 slope;
        std::array<float, 2> tileUV;
        qtplanet::Vec3 parentPositionRelativeM;
        Slope2 parentSlope;
    };

    // Height above the sphere of every packed xyz position, in meters.
    [[nodiscard]] static std::vector<float> elevationsM(const std::vector<double>& positionsM);
    [[nodiscard]] static std::vector<double> smoothNormals(const std::vector<double>& positionsM, int side);
    [[nodiscard]] static Slope2 slopeAt(const qtplanet::Vec3& positionM, const qtplanet::Vec3& normal);
    [[nodiscard]] static std::vector<double> slopesFor(const std::vector<double>& positionsM, int side);
    // Bilinear parent samples at a fine-grid vertex; absent parent uses fine samples.
    [[nodiscard]] static qtplanet::Vec3 parentPositionAt(const FineGrid& fine, const std::optional<ParentLattice>& parent, int row, int column);
    [[nodiscard]] static Slope2 parentSlopeAt(const FineGrid& fine, const std::optional<ParentLattice>& parent, int row, int column);
    [[nodiscard]] static double parentComponentAt(const ParentLattice& parent, const std::vector<double>& samples, int componentsPerSample,
                                                  int row, int column, int component);
    // Overwrites one vertex in the explicitly supplied packed output buffer.
    static void writeVertex(std::vector<float>& vertices, size_t vertexIndex, const VertexSample& sample);

    // Shared read-only DEM source; every build owns its intermediate data.
    std::shared_ptr<const MultiGeoTIFFManager> loader_;
};

#endif   // QTPLANET_TILEMESHBUILDER_H
