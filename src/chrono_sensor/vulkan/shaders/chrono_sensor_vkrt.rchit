#version 460
#extension GL_EXT_ray_tracing : require
#extension GL_EXT_buffer_reference : require
#extension GL_EXT_buffer_reference_uvec2 : require

struct GpuMaterial {
    vec4 diffuse;
    vec4 specular;
    vec4 emissive;
    vec4 params;
    uvec4 ids;
    vec4 sensor;
    uvec4 texture0;
    uvec4 texture1;
    vec4 tex_scale;
    vec4 hapke0;     // Hapke w, b, c, B_s0
    vec4 hapke1;     // Hapke h_s, phi, theta_p, reserved
};

struct GpuVertex {
    vec4 pos;
    vec4 normal;
    vec4 uv;
    vec4 tangent;
};

struct GpuTriangle {
    uvec4 index_material;
};

struct HitPayload {
    vec4 hit_tuv;        // hit distance, interpolated uv, reserved
    vec4 normal_object;  // normal xyz, object id
    vec4 tangent_flags;  // tangent xyz, has_uv
    uvec4 ids_hit;       // class id, instance id, hit flag, material index
};

layout(std430, set = 0, binding = 1) readonly buffer Materials {
    GpuMaterial materials[];
};
// A mesh instance: device addresses of its vertices and triangles, which are in the mesh's own coordinates, and
// the first of its materials (ChVulkanRTGpuInstance on the host)
struct GpuInstance {
    uvec2 vertices;
    uvec2 triangles;
    uvec4 material_base;
};
layout(buffer_reference, std430, buffer_reference_align = 16) readonly buffer VertexRef {
    GpuVertex v[];
};
layout(buffer_reference, std430, buffer_reference_align = 16) readonly buffer TriangleRef {
    GpuTriangle t[];
};
layout(std430, set = 0, binding = 2) readonly buffer Instances {
    GpuInstance instances[];
};

hitAttributeEXT vec2 attribs;
layout(location = 0) rayPayloadInEXT HitPayload payload;

void main() {
    GpuInstance inst = instances[gl_InstanceCustomIndexEXT];
    VertexRef vertices = VertexRef(inst.vertices);
    GpuTriangle tri = TriangleRef(inst.triangles).t[gl_PrimitiveID];
    uint i0 = tri.index_material.x;
    uint i1 = tri.index_material.y;
    uint i2 = tri.index_material.z;
    uint mat_id = inst.material_base.x + tri.index_material.w;

    float b1 = attribs.x;
    float b2 = attribs.y;
    float b0 = 1.0 - b1 - b2;

    GpuVertex v0 = vertices.v[i0];
    GpuVertex v1 = vertices.v[i1];
    GpuVertex v2 = vertices.v[i2];
    // Instances place meshes by rotation and translation only (a shape's scale is in its vertices), so the
    // rotation turns normals and tangents to world coordinates.
    mat3 to_world = mat3(gl_ObjectToWorldEXT);
    vec3 n = normalize(to_world * (v0.normal.xyz * b0 + v1.normal.xyz * b1 + v2.normal.xyz * b2));
    vec3 t = normalize(to_world * (v0.tangent.xyz * b0 + v1.tangent.xyz * b1 + v2.tangent.xyz * b2));
    vec2 uv = v0.uv.xy * b0 + v1.uv.xy * b1 + v2.uv.xy * b2;
    float has_uv = max(v0.uv.z, max(v1.uv.z, v2.uv.z));

    // Do not face-forward the staged shading normal based on triangle winding.
    // Chrono/OptiX keeps object/mesh normals in material space and uses them
    // directly for light sampling. The generated Vulkan primitive triangles
    // intentionally stage explicit normals, and several generated faces have a
    // winding opposite to that staged normal. Flipping here makes lit surfaces
    // point away from Point/Spot/Directional lights, producing visible material
    // colors from ambient/background only but no direct illumination or shadows.

    GpuMaterial mat = materials[mat_id];
    payload.hit_tuv = vec4(gl_HitTEXT, uv, 0.0);
    payload.normal_object = vec4(n, mat.sensor.z);
    payload.tangent_flags = vec4(t, has_uv);
    payload.ids_hit = uvec4(mat.ids.x, mat.ids.y, 1u, mat_id);
}
