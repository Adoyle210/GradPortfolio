#version 120

uniform vec2 u_resolution;
uniform float u_time;

//this is from https://iquilezles.org/articles/palettes/
vec3 palette(float t) {
    vec3 a = vec3(0.5, 0.5, 0.5);
    vec3 b = vec3(0.5, 0.5, 0.5);
    vec3 c = vec3(2.0, 1.0, 0.0);
    vec3 d = vec3(0.50, 0.20, 0.25);
    return a + b*cos(6.28318*(c*t + d));
}

//red peach green
vec3 paletteTWO(float t) {
    vec3 a = vec3(0.8, 0.5, 0.4);
    vec3 b = vec3(0.2, 0.4, 0.2);
    vec3 c = vec3(2.0, 1.0, 1.0);
    vec3 d = vec3(0.00, 0.25, 0.25);
    return a + b*cos(6.28318*(c*t + d));
}

void main() {
    vec2 uv = (gl_FragCoord.xy * 2.0 - u_resolution.xy) / u_resolution.y;
    vec2 uv0 = uv;
    vec3 finalColor = vec3(0.0);
    
    for(float i = 0.0; i < 4.0; i++) {
        uv = fract(uv * 1.5) - 0.5;
        float d = length(uv) * exp(-length(uv0));
        vec3 col = palette(length(uv0) + i*0.4 + u_time*0.4);
        
        d = sin(d*8.0 + u_time) / 8.0;
        d = abs(d);
        d = pow(0.01 / d, 1.2);
        
        finalColor += col * d;
    }
    
    gl_FragColor = vec4(finalColor, 1.0);
}
