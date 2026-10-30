#version 440
// Pathlines as screen-space ribbons of constant pixel width, plus arrowheads.
// Each vertex sits on a line point; `other` is the neighbouring point that
// gives the direction, and `offset` moves it across (x) and along (y) that
// direction in half-widths. Ribbon quads use offset (±1, 0); arrowheads use
// a tip ahead of the point and two corners behind it.
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 other;
layout(location = 2) in vec2 offset;
layout(location = 3) in float value;

layout(location = 0) out vec3 vNormal;
layout(location = 1) out float vValue;

layout(std140, binding = 0) uniform Draw {
    mat4 mvp;
    mat4 view;
    vec4 range;
    vec4 solid;
    vec4 params;
    vec4 viewport;
};

void main()
{
    vec4 a = mvp * vec4(position, 1.0);
    vec4 b = mvp * vec4(other, 1.0);
    vec2 half_px = viewport.xy * 0.5;
    vec2 sa = a.xy / a.w * half_px;
    vec2 sb = b.xy / b.w * half_px;
    vec2 d = sb - sa;
    d = length(d) > 1e-4 ? normalize(d) : vec2(1.0, 0.0);
    vec2 n = vec2(-d.y, d.x);
    a.xy += (n * offset.x + d * offset.y) * (0.5 * viewport.z) / half_px * a.w;
    vNormal = vec3(0.0, 0.0, 1.0);
    vValue = value;
    gl_Position = a;
}
