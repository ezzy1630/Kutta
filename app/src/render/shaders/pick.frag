#version 440
// GPU picking: one pixel, two float targets (world position + object, value).
layout(location = 0) in vec3 vNormal;
layout(location = 1) in float vValue;
layout(location = 2) in vec3 vWorld;
layout(location = 0) out vec4 outPosition;
layout(location = 1) out vec4 outValue;

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
    outPosition = vec4(vWorld, params.w + 1.0);
    outValue = vec4(vValue, 1.0, 0.0, 1.0);
}
