#version 440
// Triangles coloured by a per-vertex scalar. Shared by zones, slices and iso-surfaces.
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in float value;

layout(location = 0) out vec3 vNormal;
layout(location = 1) out float vValue;

layout(std140, binding = 0) uniform Draw {
    mat4 mvp;
    mat4 view;
    vec4 range;    // min, max, bands (0 = smooth), mode (0 solid colour, 1 field)
    vec4 solid;    // rgb, opacity
    vec4 params;   // lighting strength, band-line strength, missing grey, object id
    vec4 viewport; // width px, height px, line width px, glyph length
};

void main()
{
    vNormal = mat3(view) * normal;
    vValue = value;
    gl_Position = mvp * vec4(position, 1.0);
}
