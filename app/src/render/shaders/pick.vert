#version 440
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in float value;

layout(location = 0) out vec3 vNormal;
layout(location = 1) out float vValue;
layout(location = 2) out vec3 vWorld;

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
    vNormal = normal;
    vValue = value;
    vWorld = position;
    gl_Position = mvp * vec4(position, 1.0);
}
