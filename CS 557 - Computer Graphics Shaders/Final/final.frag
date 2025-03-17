#version 120

// Input uniforms from OpenGL
uniform vec2 u_resolution;  // Screen resolution
uniform float u_time;       // Time in seconds
uniform float u_zoom;       // Zoom level - controlled by keyboard

// Shape distance functions for creating geometric shapes

// Creates an uneven capsule shape (pill-like shape with different end radii)
float sdUnevenCapsule( vec2 p, float r1, float r2, float h )
{
    p.x = abs(p.x);
    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(p,vec2(-b,a));
    if( k < 0.0 ) return length(p) - r1;
    if( k > a*h ) return length(p-vec2(0.0,h)) - r2;
    return dot(p, vec2(a,b) ) - r1;
}

// Creates a crescent moon shape using two circles
float sdMoon(vec2 p, float d, float ra, float rb) {
    p.y -= ra*0.25;  // Offset the moon position
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra-a*a,0.0));
    
    vec2 p2 = vec2(p.x, abs(p.y));
    float d1 = length(p2-vec2(a,0)) - rb;  // Inner circle
    float d2 = length(p2+vec2(a,0)) - ra;  // Outer circle
    return max(d1,-d2);  // Combine circles to create crescent
}

// Creates a five-pointed star shape
float sdStar(vec2 p, float r, float rf) {
    float an = 3.141593/5.0;  // Angle for pentagon
    float en = 3.141593/2.5;  // Angle for star points
    vec2 acs = vec2(cos(an),sin(an));
    vec2 ecs = vec2(cos(en),sin(en));
    float bn = mod(atan(p.x,p.y),2.0*an) - an;
    p = length(p)*vec2(cos(bn),abs(sin(bn)));
    p -= r*acs;
    p += ecs*clamp( -dot(p,ecs), 0.0, r*acs.y/ecs.y);
    return length(p)*sign(p.x);
}

// Creates a smooth color palette using cosine waves
// From Inigo Quilez's color palette formula
vec3 palette(float t) {
    vec3 a = vec3(0.5, 0.5, 0.5);    // Offset
    vec3 b = vec3(0.5, 0.5, 0.5);    // Amplitude
    vec3 c = vec3(2.0, 1.0, 0.0);    // Frequency
    vec3 d = vec3(0.50, 0.20, 0.25); // Phase
    return a + b*cos(6.28318*(c*t + d));
}

// Smooth minimum function for blending distances
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5*(b-a)/k, 0.0, 1.0);
    return mix(b, a, h) - k*h*(1.0-h);
}

void main() {
    // Convert pixel coordinates to normalized space (-0.5 to 0.5)
    // Scale based on zoom uniform
    vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution.xy) / min(u_resolution.x, u_resolution.y) * u_zoom;
    vec2 uv0 = uv;  // Store original coordinates for later use
    vec3 finalColor = vec3(0.0);
    
    // Calculate distances to celestial shapes
    float moonDist = sdMoon(uv0 - vec2(0.125, 0.075), 0.2, 0.15, 0.125);
    float starDist = sdStar(uv0 + vec2(0.075, -0.05), 0.05, 0.025);
    
    // Create rotating capsule
    vec2 capsulePos = uv0 + vec2(-0.05, 0.1);
    // Rotate the capsule over time
    capsulePos = mat2(cos(u_time), -sin(u_time), sin(u_time), cos(u_time)) * capsulePos;
    float capsuleDist = sdUnevenCapsule(capsulePos, 0.075, 0.025, 0.1);
    
    // Fractal iteration loop
    for(float i = 0.0; i < 4.0; i++) {
        // Create fractal pattern by scaling and repeating
        uv = fract(uv * 1.5) - 0.5;
        // Calculate base distance field with exponential falloff
        float d = length(uv) * exp(-length(uv0));
        
        // Blend all celestial shapes together
        float celestialDist = smin(moonDist, starDist, 0.3);
        celestialDist = smin(celestialDist, capsuleDist, 0.3);
        // Mix celestial shapes with base pattern
        d = mix(d, abs(celestialDist), 0.2);
        
        // Generate colors based on position and time
        vec3 col = palette(length(uv0) + i*0.4 + u_time*0.4);
        
        // Create final pattern by modulating distance field
        d = sin(d*8.0 + u_time)/8.0;   // Create wave pattern
        d = abs(d);                    // Make pattern symmetrical
        d = pow(0.01 / d, 1.2);       // Create sharp edges
        
        // Accumulate colors weighted by pattern
        finalColor += col * d;
    }
    
    // Output final color
    gl_FragColor = vec4(finalColor, 1.0);
}
