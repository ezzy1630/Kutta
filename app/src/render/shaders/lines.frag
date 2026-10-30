#version 440
layout(location = 0) out vec4 fragColor;

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
    fragColor = vec4(solid.rgb * solid.a, solid.a);
}
