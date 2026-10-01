// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Florian Reinle
// =============================================================================
// Vulkan RT scene staging layer.
// =============================================================================

#ifndef CH_VULKAN_RT_SCENE_H
#define CH_VULKAN_RT_SCENE_H

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "chrono/assets/ChColor.h"
#include "chrono/assets/ChVisualMaterial.h"

#include "chrono/assets/ChVisualModel.h"
#include "chrono/assets/ChVisualShapes.h"

#include "chrono/core/ChFrame.h"
#include "chrono/physics/ChSystem.h"
#include "chrono_sensor/ChApiSensor.h"
#include "chrono_sensor/ChSensorRenderTypes.h"

namespace chrono {
namespace sensor {

/// @addtogroup sensor_vulkan
/// @{

struct CH_SENSOR_API ChVulkanRTSceneStats {
    uint32_t bodies = 0;
    uint32_t other_items = 0;
    uint32_t visible_shapes = 0;
    uint32_t boxes = 0;
    uint32_t spheres = 0;
    uint32_t cylinders = 0;
    uint32_t triangle_meshes = 0;
    uint32_t unsupported_shapes = 0;
};

enum class ChVulkanRTPrimitiveType { BOX, SPHERE, CYLINDER, TRIANGLE_MESH, MESH_PROXY };

struct CH_SENSOR_API ChVulkanRTTexCoord {
    float u = 0.f;
    float v = 0.f;
};

struct CH_SENSOR_API ChVulkanRTMaterial {
    /// Construct the implicit Vulkan material from Chrono's canonical visual default.
    /// Keeping this conversion centralized prevents the backend from drifting when
    /// ChVisualMaterial::Default() changes.
    ChVulkanRTMaterial();

    ChVector3f diffuse;
    ChVector3f ambient;
    ChVector3f specular;
    ChVector3f emissive;

    // Mirrors the fields consumed by the OptiX camera shaders. OptiX names this
    // value "transparency", but it is used as opacity/surface weight: 1 is
    // opaque, 0 is fully transparent and traces through the surface.
    float opacity;
    float roughness;
    float metallic;
    float emissive_power;
    float shininess;
    bool use_specular_workflow;

    // Non-camera sensor response parameters. OptiX keeps these in its material
    // record; mirror them here so LiDAR/Radar parity does not depend on OptiX.
    float lidar_intensity = 1.f;
    float radar_backscatter = 1.f;

    // Camera BSDF. Vulkan RT shades HAPKE with the Hapke model and every other type
    // with the OptiX LEGACY model, as the OptiX camera shader does.
    BSDFType bsdf_type;
    float hapke_w = 0.f;        // single scattering albedo
    float hapke_b = 0.f;        // phase function shape
    float hapke_c = 0.f;        // backward/forward scattering weight
    float hapke_B_s0 = 0.f;     // shadow hiding opposition amplitude
    float hapke_h_s = 0.f;      // shadow hiding opposition width
    float hapke_phi = 0.f;      // filling factor
    float hapke_theta_p = 0.f;  // macroscopic roughness (rad)

    float tex_scale_u;
    float tex_scale_v;
    std::string diffuse_texture;
    std::string specular_texture;
    std::string emissive_texture;
    std::string normal_texture;
    std::string roughness_texture;
    std::string metallic_texture;
    std::string opacity_texture;
    std::string weight_texture;

    unsigned short int class_id;
    unsigned short int instance_id;
};

struct CH_SENSOR_API ChVulkanRTTriangle {
    ChVector3d v0 = ChVector3d(0.0, 0.0, 0.0);
    ChVector3d v1 = ChVector3d(0.0, 0.0, 0.0);
    ChVector3d v2 = ChVector3d(0.0, 0.0, 0.0);
    ChVector3d normal = ChVector3d(0.0, 0.0, 1.0);
    ChVector3d n0 = ChVector3d(0.0, 0.0, 1.0);
    ChVector3d n1 = ChVector3d(0.0, 0.0, 1.0);
    ChVector3d n2 = ChVector3d(0.0, 0.0, 1.0);
    ChVulkanRTTexCoord uv0;
    ChVulkanRTTexCoord uv1;
    ChVulkanRTTexCoord uv2;
    ChVector3d tangent = ChVector3d(1.0, 0.0, 0.0);
    bool has_vertex_normals = false;
    bool has_uvs = false;

    /// Index of the triangle's material in ChVulkanRTPrimitive::materials.
    uint32_t material_index = 0;

    // Source-mesh indices of each corner's position, normal and UV, or -1 where the corner's value does not
    // come from the mesh (a face normal, or no UVs). Corners with equal indices have equal positions,
    // normals and UVs, which lets the GPU path share their vertex.
    int32_t position_index[3] = {-1, -1, -1};
    int32_t normal_index[3] = {-1, -1, -1};
    int32_t uv_index[3] = {-1, -1, -1};
};

struct CH_SENSOR_API ChVulkanRTPrimitive {
    ChVulkanRTPrimitiveType type = ChVulkanRTPrimitiveType::BOX;
    ChFrame<double> frame;
    ChVector3d scale = ChVector3d(1.0, 1.0, 1.0);
    ChVulkanRTMaterial material;

    // Triangle-mesh data in primitive-local coordinates. This mirrors the OptiX
    // triangle GAS input: vertices are kept untransformed except for the visual
    // shape scale; the body/asset transform remains in frame. The triangles are
    // immutable once staged and shared between syncs while the mesh is unchanged.
    std::shared_ptr<const std::vector<ChVulkanRTTriangle>> triangles;

    /// The staged triangles (empty for a primitive that is not a triangle mesh).
    const std::vector<ChVulkanRTTriangle>& Triangles() const {
        static const std::vector<ChVulkanRTTriangle> none;
        return triangles ? *triangles : none;
    }

    /// Materials of a triangle mesh's triangles, indexed by ChVulkanRTTriangle::material_index. Entry 0 is
    /// the primitive's own material, for triangles without a valid material index of their own.
    std::vector<ChVulkanRTMaterial> materials;

    /// The material a triangle of this primitive is shaded with.
    const ChVulkanRTMaterial& TriangleMaterial(const ChVulkanRTTriangle& tri) const {
        return tri.material_index < materials.size() ? materials[tri.material_index] : material;
    }

    ChVector3d aabb_min = ChVector3d(0.0, 0.0, 0.0);
    ChVector3d aabb_max = ChVector3d(0.0, 0.0, 0.0);
    bool has_aabb = false;
    bool backface_cull = false;

    // Body motion staged for radar Doppler/velocity returns.  Static/non-body
    // items stay at zero, matching the OptiX convention for stationary hits.
    ChVector3d translational_velocity = ChVector3d(0.0, 0.0, 0.0);
    ChVector3d angular_velocity = ChVector3d(0.0, 0.0, 0.0);
    float object_id = 0.f;
};

struct CH_SENSOR_API ChVulkanRTLight {
    LightType type = LightType::POINT_LIGHT;
    ChVector3f pos = ChVector3f(0.f, 0.f, 0.f);
    ChVector3f dir = ChVector3f(0.f, 0.f, -1.f);
    ChVector3f color = ChVector3f(1.f, 1.f, 1.f);
    float range = 100.f;
    float angle = 0.f;
    bool const_color = true;
    float atten_scale = 1.f;
    float angle_falloff_start = 0.f;
    float angle_atten_rate = -1.f;

    // Area-light parameters. Vulkan RT currently evaluates area lights from
    // their center point, but keeping the full public data here preserves the
    // OptiX ChOptixScene API and allows a renderer upgrade without API changes.
    ChVector3f length_vec = ChVector3f(0.f, 0.f, 0.f);
    ChVector3f width_vec = ChVector3f(0.f, 0.f, 0.f);
    float radius = 0.f;
    float area = 0.f;
    /// A directional light's disk, as the cameras see it where they look toward the light past everything else: its
    /// angular radius (rad; 0, the default, for none). The disk is as bright as the light's color over its solid angle
    /// makes it, darkened toward its limb as the Sun's is, and is seen only along a camera's own rays (not in mirrors,
    /// whose highlights the light's shading gives already). Its light on the scene is unchanged: a direction still
    /// casts sharp shadows.
    float disk_radius = 0.f;

    std::string texture;
};

/// A participating medium on a voxel grid, such as lofted dust, drawn by camera sensors with single scattering of
/// the scene's lights. Voxel (i, j, k) is centered at origin + voxel * (i + 1/2, j + 1/2, k + 1/2) and stored at
/// index i + nx * (j + ny * k); values are trilinear between voxel centers.
///
/// Scattering follows the particle part of the Hapke model, so a medium made of the same grains as a Hapke surface
/// looks like it: single-scattering albedo w and the double Henyey-Greenstein phase function with shape b and
/// back/forward weight c, as ChVisualMaterial::SetHapkeParameters takes them, tinted by color.
///
/// Light reaching the medium is checked against the scene's geometry with a shadow ray, and the light from one
/// directional light, the one along sun_dir, can also be given per voxel: its transmittance through the medium (the
/// medium's own shadow, which the shadow rays do not see) and its visibility past other geometry, such as terrain
/// outside the scene, which also spares the shadow rays where it is zero. Surfaces are shadowed by the medium from
/// every light.
struct ChVulkanRTVolume {
    ChVector3f origin = ChVector3f(0.f, 0.f, 0.f);  ///< corner of the grid (m)
    float voxel = 0.1f;                             ///< voxel edge (m)
    unsigned int nx = 0, ny = 0, nz = 0;            ///< voxels along x, y and z
    std::vector<float> extinction;                  ///< extinction coefficient (1/m)
    std::vector<float> sun_transmittance;           ///< toward sun_dir, through the medium; empty for none
    std::vector<float> sun_visibility;              ///< toward sun_dir, past other geometry; empty for full
    ChVector3f sun_dir = ChVector3f(0.f, 0.f, 1.f);  ///< unit direction toward the light the two fields are for
    float albedo = 0.3f;                            ///< Hapke single-scattering albedo w
    float phase_b = 0.25f;                          ///< Hapke phase function shape b
    float phase_c = 0.3f;                           ///< Hapke phase function back/forward weight c
    ChVector3f color = ChVector3f(1.f, 1.f, 1.f);   ///< tint of the scattered light
    /// The grid's axes in the scene, so a thin layer such as dust over sloping ground needs no grid as tall as it is
    /// wide. `origin` is its corner in the scene; voxel (i, j, k) is at origin + rotation (voxel (i + 1/2, ...)).
    ChQuaternionf rotation = ChQuaternionf(1.f, 0.f, 0.f, 0.f);
};

/// A planet's atmosphere, drawn by the cameras as light scattered along their rays (Rayleigh, Mie and ozone
/// absorption, single scattering, lit by the scene's directional lights and shadowed by the planet), and the frame a
/// PLANET material is shaded in: its center and axes, so that its clouds' shadows and its night side are found where
/// they are. Coefficients are per meter at the ground, falling off exponentially with height (ozone: a tent about
/// its center). The defaults are the Earth's (Bruneton & Neyret 2008; Hillaire 2020).
struct ChVulkanRTAtmosphere {
    ChVector3d center = ChVector3d(0, 0, 0);          ///< the planet's center, in the scene (m)
    ChQuaterniond rotation = QUNIT;                   ///< the planet's axes in the scene's
    double planet_radius = 6371e3;                    ///< the ground (m)
    double top_radius = 6471e3;                       ///< the top of the air (m)
    ChVector3f rayleigh = ChVector3f(5.802e-6f, 13.558e-6f, 33.1e-6f);  ///< Rayleigh scattering (1/m), red green blue
    float rayleigh_height = 8000.f;                   ///< its scale height (m)
    float mie_scattering = 3.996e-6f;                 ///< Mie (aerosol) scattering (1/m)
    float mie_extinction = 4.44e-6f;                  ///< Mie extinction (1/m)
    float mie_height = 1200.f;                        ///< its scale height (m)
    float mie_g = 0.8f;                               ///< its Henyey-Greenstein asymmetry
    ChVector3f ozone = ChVector3f(0.650e-6f, 1.881e-6f, 0.085e-6f);  ///< ozone absorption at its peak (1/m)
    float ozone_center = 25000.f;                     ///< its peak's height (m)
    float ozone_width = 15000.f;                      ///< its half width (m)
    float cloud_height = 6000.f;                      ///< the height a PLANET material's clouds cast their shadows from (m)
};

/// The stars, as the cameras see them past everything else: each a point of light from a direction in the scene's
/// axes, baked into a map of the sky (equirectangular, `width` × width / 2 texels, each star's light spread over the
/// four texels about it as its radiance there), which the cameras' rays that meet nothing take their light from. A
/// star's irradiance is in the units the lights' colors are in, times `scale`: with a directional light for the Sun
/// of color 1, a star of visual magnitude m gives 10^(-0.4 (m + 26.74)), and `scale` is then the Sun's color, so the
/// stars are as bright beside the sunlit ground as they are. A camera exposed for sunlit ground sees the brightest few
/// at most, as a real one does; one exposed for the stars (a star tracker's) sees them all.
struct ChVulkanRTStarField {
    std::vector<ChVector3f> directions;  ///< unit, in the scene's axes
    std::vector<ChVector3f> irradiance;  ///< rgb, before `scale`
    float scale = 1.f;
    unsigned int width = 8192;           ///< the map's texels about the sky: 8192, 0.044° each
};

/// Staging scene for the Vulkan backend.
///
/// This object mirrors the public ChOptixScene methods used by existing Sensor demos
/// while collecting a compact renderable representation. The representation is
/// intentionally independent from OptiX/CUDA and can feed either the Vulkan RT
/// BLAS/TLAS builder or the host fallback used for bring-up/testing.
class CH_SENSOR_API ChVulkanRTScene {
  public:
    ChVulkanRTScene() = default;

    void SyncFromSystem(ChSystem* system);

    const ChVulkanRTSceneStats& GetStats() const { return m_stats; }
    uint64_t GetRevision() const { return m_revision; }
    const std::vector<ChVulkanRTPrimitive>& GetPrimitives() const { return m_primitives; }

    void SetAmbientLight(const ChVector3f& color);
    const ChVector3f& GetAmbientLight() const { return m_ambient_light; }

    void SetBackground(const Background& background);
    const Background& GetBackground() const { return m_background; }

    unsigned int AddPointLight(ChVector3f pos, ChColor color, float max_range, bool const_color = true);

    unsigned int AddDirectionalLight(const ChVector3f& dir, const ChVector3f& color);
    /// A directional light from spherical angles; `disk_radius` the angular radius of its disk as the cameras see it (rad; see
    /// ChVulkanRTLight::disk_radius)
    unsigned int AddDirectionalLight(ChColor color, float elevation, float azimuth, float disk_radius = 0.f);

    unsigned int AddSpotLight(const ChVector3f& pos, const ChVector3f& dir, const ChVector3f& color, float range, float angle);
    unsigned int AddSpotLight(ChVector3f pos,
                              ChColor color,
                              float max_range,
                              ChVector3f light_dir,
                              float angle_falloff_start,
                              float angle_range,
                              bool const_color = true);

    unsigned int AddRectangleLight(ChVector3f pos,
                                   ChColor color,
                                   float max_range,
                                   ChVector3f length_vec,
                                   ChVector3f width_vec,
                                   bool const_color = true);

    unsigned int AddDiskLight(ChVector3f pos,
                              ChColor color,
                              float max_range,
                              ChVector3f light_dir,
                              float radius,
                              bool const_color = true);

    // Redefine an existing light in place. Each takes the ID returned by the matching
    // Add*Light plus that function's own parameters, mirroring ChOptixScene's backend-neutral
    // Modify*Light overloads. Out-of-range IDs are ignored.
    void ModifyPointLight(unsigned int light_ID, ChVector3f pos, ChColor color, float max_range, bool const_color = true);

    void ModifyDirectionalLight(unsigned int light_ID, ChColor color, float elevation, float azimuth, float disk_radius = 0.f);

    void ModifySpotLight(unsigned int light_ID,
                         ChVector3f pos,
                         ChColor color,
                         float max_range,
                         ChVector3f light_dir,
                         float angle_falloff_start,
                         float angle_range,
                         bool const_color = true);

    void ModifyRectangleLight(unsigned int light_ID, ChVector3f pos, ChColor color, float max_range, ChVector3f length_vec, ChVector3f width_vec, bool const_color = true);

    void ModifyDiskLight(unsigned int light_ID, ChVector3f pos, ChColor color, float max_range, ChVector3f light_dir, float radius, bool const_color = true);

    unsigned int AddEnvironmentLight(const std::string& env_tex, const ChVector3f& color = ChVector3f(1.f, 1.f, 1.f));
    unsigned int AddEnvironmentLight(std::string env_tex_path, float intensity_scale);

    void SetLights(const std::vector<ChVulkanRTLight>& lights);

    /// Set the participating medium the cameras draw (null for none). The renderers upload it again whenever it
    /// is set, without rebuilding the scene's geometry, so set a new one each time it changes.
    void SetVolume(std::shared_ptr<const ChVulkanRTVolume> volume) {
        m_volumes.clear();
        if (volume)
            m_volumes.push_back(std::move(volume));
        ++m_volume_revision;
    }
    /// Set several media, each on its own grid with its own scattering, such as a jet's gas and the dust it raises:
    /// up to four, the rest ignored. Each is marched on its own along a camera ray and they are laid over each other,
    /// near enough where thin media overlap. Surfaces are shadowed by all of them.
    void SetVolumes(std::vector<std::shared_ptr<const ChVulkanRTVolume>> volumes) {
        m_volumes = std::move(volumes);
        ++m_volume_revision;
    }
    const std::vector<std::shared_ptr<const ChVulkanRTVolume>>& GetVolumes() const { return m_volumes; }
    uint64_t GetVolumeRevision() const { return m_volume_revision; }
    /// Set the planet's atmosphere the cameras draw (null for none)
    void SetAtmosphere(std::shared_ptr<const ChVulkanRTAtmosphere> atmosphere) {
        m_atmosphere = std::move(atmosphere);
        Touch();
    }
    const std::shared_ptr<const ChVulkanRTAtmosphere>& GetAtmosphere() const { return m_atmosphere; }

    /// Set the stars the cameras see (null for none). The map is baked once for each star field set: to change only
    /// their brightness, set the same stars with another scale through SetStarScale
    void SetStars(std::shared_ptr<const ChVulkanRTStarField> stars) {
        m_stars = std::move(stars);
        m_star_scale = m_stars ? m_stars->scale : 1.f;
        Touch();
    }
    void SetStarScale(float scale) {
        m_star_scale = scale;
        Touch();
    }
    const std::shared_ptr<const ChVulkanRTStarField>& GetStars() const { return m_stars; }
    float GetStarScale() const { return m_star_scale; }
    void ClearLights() { if (!m_lights.empty()) { m_lights.clear(); Touch(); } }
    const std::vector<ChVulkanRTLight>& GetLights() const { return m_lights; }

  private:
    void Touch() { ++m_revision; }

    /// Append a light and return its ID (its index).
    unsigned int Append(const ChVulkanRTLight& light);

    /// Overwrite the light at `id`; out-of-range IDs are ignored.
    void Replace(unsigned int id, const ChVulkanRTLight& light);

    ChVulkanRTMaterial ExtractMaterial(const std::shared_ptr<ChVisualShape>& shape) const;

    /// Stage a triangle mesh into primitive, reusing the triangles staged for it in an earlier sync when
    /// the mesh is not mutable and its scale, face culling and materials are unchanged. Mutable meshes
    /// are edited in place, so they are staged again every time. Returns false for an empty mesh.
    bool StageMesh(ChVulkanRTPrimitive& primitive,
                   const std::shared_ptr<ChTriangleMeshConnected>& mesh,
                   const ChVector3d& scale,
                   bool backface_cull,
                   const std::vector<ChVulkanRTMaterial>& materials,
                   bool is_mutable);

    /// Triangles staged for a mesh, kept while the mesh stays in the scene.
    struct MeshCacheEntry {
        std::shared_ptr<ChTriangleMeshConnected> mesh;  ///< held, so the mesh's address is not reused while cached
        ChVector3d scale;
        bool backface_cull = false;
        std::vector<ChVulkanRTMaterial> materials;  ///< the staged primitive's material table
        std::shared_ptr<const std::vector<ChVulkanRTTriangle>> triangles;
        ChVector3d aabb_min;
        ChVector3d aabb_max;
        bool has_aabb = false;
        uint64_t last_sync = 0;  ///< last sync that used the entry
    };
    std::unordered_map<const ChTriangleMeshConnected*, std::vector<MeshCacheEntry>> m_mesh_cache;
    uint64_t m_sync_count = 0;

    ChVulkanRTSceneStats m_stats;
    std::vector<ChVulkanRTPrimitive> m_primitives;
    std::vector<ChVulkanRTLight> m_lights;
    ChVector3f m_ambient_light = ChVector3f(0.1f, 0.1f, 0.1f);
    Background m_background;
    uint64_t m_revision = 1;
    uint64_t m_system_signature = 0;
    std::vector<std::shared_ptr<const ChVulkanRTVolume>> m_volumes;
    uint64_t m_volume_revision = 1;
    std::shared_ptr<const ChVulkanRTAtmosphere> m_atmosphere;
    std::shared_ptr<const ChVulkanRTStarField> m_stars;
    float m_star_scale = 1.f;
};

/// @} sensor_vulkan

}  // namespace sensor
}  // namespace chrono

#endif
