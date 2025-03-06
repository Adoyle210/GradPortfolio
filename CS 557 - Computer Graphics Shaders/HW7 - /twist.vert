#version 120

uniform float uTwist;  // Twist amount from C++

varying vec2 vST;      // texture coords
varying vec3 vN;       // normal vector
varying vec3 vL;       // vector from point to light
varying vec3 vE;       // vector from point to eye
varying vec3 vMC;      // model coordinates for hatching
varying float vZ;      // eye coordinate Z for ChromaDepth

// Rotation function for Y-axis twist
vec3 RotateY(vec3 xyz, float radians)
{
    float c = cos(radians);
    float s = sin(radians);
    vec3 newxyz = xyz;
    newxyz.xz = vec2(
        dot(xyz.xz, vec2(c, s)),
        dot(xyz.xz, vec2(-s, c))
    );
    return newxyz;
}

void main()
{
    // Get the vertex position in model coordinates
    vMC = gl_Vertex.xyz;
    
    // Calculate twist angle based on Y position
    float twistAngle = uTwist * vMC.y;
    
    // Apply the twist
    vec4 twistedVertex = vec4(RotateY(vMC, twistAngle), 1.0);
    
    // Transform vertex to eye coordinates
    vec4 ECposition = gl_ModelViewMatrix * twistedVertex;
    
    // Save the eye coordinate Z for ChromaDepth
    vZ = -ECposition.z;
    
    // Transform the normal
    vec3 twistedNormal = RotateY(gl_Normal, twistAngle);
    vN = normalize(gl_NormalMatrix * twistedNormal);
    
    // Calculate vectors for lighting
    vec3 lightPos = vec3(5.0, 5.0, 10.0);
    vL = normalize(lightPos - ECposition.xyz);
    vE = normalize(-ECposition.xyz);
    
    // Pass through texture coordinates
    vST = gl_MultiTexCoord0.st;
    
    // Final position
    gl_Position = gl_ModelViewProjectionMatrix * twistedVertex;
}
