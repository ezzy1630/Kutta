#version 440
layout(location = 0) in vec3 vNormal;
layout(location = 1) in float vValue;
layout(location = 0) out vec4 fragColor;

layout(std140, binding = 0) uniform Draw {
    mat4 mvp;
    mat4 view;
    vec4 range;
    vec4 solid;
    vec4 params;
    vec4 viewport;
};
layout(binding = 1) uniform sampler2D colormap;

void main()
{
    vec3 base = solid.rgb;
    if (range.w > 0.5) {
        if (vValue > 1.0e37) {
            // Missing sample (no Fluent value here): never coloured as data.
            base = vec3(params.z);
        } else {
            float t = clamp((vValue - range.x) / max(range.y - range.x, 1e-30), 0.0, 1.0);
            float bands = range.z;
            if (bands > 0.5) {
                float b = min(floor(t * bands), bands - 1.0);
                base = texture(colormap, vec2((b + 0.5) / bands, 0.5)).rgb;
                // Thin contour line on each band edge, constant width in pixels.
                float f = t * bands;
                float d = abs(fract(f + 0.5) - 0.5);
                float w = max(fwidth(f), 1e-6);
                float line = 1.0 - smoothstep(0.35 * w, 1.1 * w, d);
                if (t > 0.0 && t < 1.0) base = mix(base, base * 0.45, line * params.y);
            } else {
                base = texture(colormap, vec2(t, 0.5)).rgb;
            }
        }
    }
    // Two-sided headlight; params.x = 0 renders flat (slices, ribbons).
    vec3 n = normalize(vNormal);
    float ndl = abs(n.z);
    float light = mix(1.0, 0.42 + 0.58 * ndl + 0.10 * pow(ndl, 24.0), params.x);
    fragColor = vec4(base * light * solid.a, solid.a);
}
