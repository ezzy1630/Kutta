#version 440
// Velocity vectors: an arrow mesh along +x, instanced per sample point.
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 base;     // per instance
layout(location = 3) in vec3 vector;   // per instance
layout(location = 4) in float value;   // per instance (colour scalar)

layout(location = 0) out vec3 vNormal;
layout(location = 1) out float vValue;

layout(std140, binding = 0) uniform Draw {
    mat4 mvp;
    mat4 view;
    vec4 range;
    vec4 solid;
    vec4 params;
    vec4 viewport;   // .w = arrow length in world units for the reference magnitude
};

void main()
{
    float mag = length(vector);
    vec3 x = mag > 1e-20 ? vector / mag : vec3(1.0, 0.0, 0.0);
    vec3 helper = abs(x.z) < 0.9 ? vec3(0.0, 0.0, 1.0) : vec3(0.0, 1.0, 0.0);
    vec3 y = normalize(cross(helper, x));
    vec3 z = cross(x, y);
    // Length scales with magnitude relative to the range maximum; width is fixed.
    float len = viewport.w * clamp(mag / max(range.y, 1e-30), 0.15, 1.5);
    float wid = viewport.w * 0.16;
    vec3 p = base + x * (position.x * len) + y * (position.y * wid) + z * (position.z * wid);
    vNormal = mat3(view) * (x * normal.x + y * normal.y + z * normal.z);
    vValue = value;
    gl_Position = mvp * vec4(p, 1.0);
}
