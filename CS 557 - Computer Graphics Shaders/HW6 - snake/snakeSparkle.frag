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

// Color palette function
vec3 palette(float t) {
    vec3 a = vec3(0.5, 0.5, 0.5);       //grey
    vec3 b = vec3(0.5, 0.5, 0.5);       //grey
    //vec3 b = vec3(1.0, 0.0, 0.0);       //red
    vec3 c = vec3(1.0, 1.0, 1.0);      //white
    vec3 d = vec3(0.263, 0.416, 0.557); //blue
    return a + b * cos(6.28318 * (c * t + d));
}

// Function to create a star shape
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
    // Create UV coordinates from the normal
    vec3 pos = normalize(vN);
    vec2 uv = vec2(
        atan(pos.z, pos.x) / (2.0 * PI) + 0.5,
        acos(pos.y) / PI
    );
    
    // Create repeating star pattern
    vec2 scaled_uv = fract(uv * 8.0) - 0.5;  // Create 8x8 grid of stars
    float star = starShape(scaled_uv, 0.4);
    
    // Create color from palette
    vec3 starColor = palette(star + iTime * 0.5);
    
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
