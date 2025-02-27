#version 120  // Add version declaration

// Input uniforms from C++ program
uniform float iTime;
uniform float uKa, uKd, uKs;        // lighting coefficients
uniform vec4 uSpecularColor;        // specular highlight color
uniform float uShininess;           // specular exponent

// Input varyings from vertex shader
varying vec2 vST;        // texture coords
varying vec3 vN;         // normal vector
varying vec3 vL;         // vector to light
varying vec3 vE;         // vector to eye

const float PI = 3.14159265;

float sdStar(vec2 p, float r, int n, float m)
{
    float an = 3.141593/float(n);
    float en = 3.141593/m;
    vec2 acs = vec2(cos(an),sin(an));
    vec2 ecs = vec2(cos(en),sin(en));
    float bn = mod(atan(p.x,p.y),2.0*an) - an;
    p = length(p)*vec2(cos(bn),abs(sin(bn)));
    p -= r*acs;
    p += ecs*clamp( -dot(p,ecs), 0.0, r*acs.y/ecs.y);
    return length(p)*sign(p.x);
}

vec2 snakePath(float t) {
    return vec2(
        sin(t * 0.1) * 0.5,
        cos(t * 0.2) * 0.3
    );
}

// Color palette function
vec3 palette(float t) {
    vec3 a = vec3(0.5, 0.5, 0.5);       //grey
    vec3 b = vec3(0.5, 0.5, 0.5);       //grey
    //vec3 b = vec3(1.0, 0.0, 0.0);       //red
    vec3 c = vec3(1.0, 1.0, 1.0);      //white
    vec3 d = vec3(0.263, 0.416, 0.557); //blue
    return a + b * cos(6.28318 * (c * t + d));
}

// Function to create a star pattern
float starPattern(vec2 st, float size, int points) 
{
    vec2 pos = st - vec2(0.5);
    float r = length(pos) * 2.0;
    float a = atan(pos.y, pos.x);
    float f = PI * 2.0 / float(points);
    float star = cos(a * float(points)) * size;
    return 1.0 - smoothstep(star, star + 0.01, r);
}

// Function to create a star shape
// uv: 2D coordinates centered at origin
// size: controls the size of the star
float starShape(vec2 uv, float size) {
    float a = atan(uv.y, uv.x);  // Angle from x-axis
    float r = length(uv);        // Distance from center
    
    // Create a 5-pointed star using cosine
    float stars = cos(a * 5.0) * 0.5 + 0.5;
    // Create sharp edges for the star shape
    float star = smoothstep(stars * size, stars * size + 0.01, r);
    
    return 1.0 - star;  // Invert so star is bright on dark background
}

void main()
{
    // Use texture coordinates directly, scaled to create more stars along the snake
    vec2 uv = vST * vec2(16.0, 4.0);  // More repetitions horizontally than vertically
    
    // Create repeating star pattern
    vec2 scaled_uv = fract(uv) - 0.5;  // Create repeating grid of stars
    float star = starShape(scaled_uv, 0.4);
    
    // Create color from palette with position-based offset
    float colorOffset = vST.x + vST.y + iTime * 0.5;  // Make colors flow along the snake
    vec3 starColor = palette(star + colorOffset);
    
    // Add some glow
    float glow = exp(-length(scaled_uv) * 4.0) * 0.5;
    vec3 finalColor = mix(starColor, vec3(1.0), glow);
    
    // Basic lighting
    vec3 Normal = normalize(vN);
    vec3 Light = normalize(vL);
    vec3 Eye = normalize(vE);
    
    // Calculate lighting components
    float diff = max(dot(Normal, Light), 0.0);
    vec3 ambient = uKa * finalColor;
    vec3 diffuse = uKd * diff * finalColor;
    
    // Specular component
    vec3 specular = vec3(0.0);
    if(diff > 0.0) {
        vec3 ref = normalize(reflect(-Light, Normal));
        float spec = pow(max(dot(Eye, ref), 0.0), uShininess);
        specular = uKs * spec * uSpecularColor.rgb;
    }
    
    // Combine lighting with pattern
    finalColor = ambient + diffuse + specular;
    finalColor = clamp(finalColor, 0.0, 1.0);
    
    gl_FragColor = vec4(finalColor, 1.0);
}
