#version 440
// Mesh edges and outlines: one colour, drawn over surfaces they belong to.
layout(location = 0) in vec3 position;

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
    gl_Position = mvp * vec4(position, 1.0);
}
