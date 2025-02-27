#version 120  // Add version declaration

// These variables will be passed to the fragment shader
varying vec2 vST;        // Texture coordinates (s,t) for mapping patterns
varying vec3 vN;         // Normal vector for lighting calculations
varying vec3 vL;         // Vector from vertex to light source
varying vec3 vE;         // Vector from vertex to eye/camera

void main() 
{
    // Get texture coordinates from OpenGL
    vST = gl_MultiTexCoord0.st;
    
    // Transform vertex position from model space to eye space
    vec4 ECposition = gl_ModelViewMatrix * gl_Vertex;
    
    // Transform the normal vector to eye space and normalize it
    vN = normalize(gl_NormalMatrix * gl_Normal);
    
    // Define light position in eye coordinates
    vec4 lightPos = vec4(5.0, 5.0, 10.0, 1.0);
    vL = normalize(lightPos.xyz - ECposition.xyz);
    
    // Calculate vector from vertex to eye
    vE = normalize(-ECposition.xyz);
    
    // Transform vertex to clip space for final position
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
