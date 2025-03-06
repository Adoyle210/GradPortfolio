#version 120

uniform float uKa, uKd, uKs;
uniform vec4 uSpecularColor;
uniform float uShininess;
uniform float uRedDepth, uBlueDepth;
uniform int uUseChromaDepth;

varying vec2 vST;      // texture coords
varying vec3 vN;       // normal vector
varying vec3 vL;       // vector from point to light
varying vec3 vE;       // vector from point to eye
varying vec3 vMC;      // model coordinates for hatching
varying float vZ;      // eye coordinate Z for ChromaDepth

// Smooth pulse function for hatching
float SmoothPulse(float left, float right, float value, float tol)
{
    return smoothstep(left-tol, left+tol, value) - 
           smoothstep(right-tol, right+tol, value);
}

// Rainbow function for ChromaDepth
vec3 Rainbow(float t)
{
    t = clamp(t, 0., 1.);         // 0.00 is red, 0.33 is green, 0.67 is blue

    float r = 1.;
    float g = 0.0;
    float b = 1. - 6. * (t - (5./6.));

    if(t <= (5./6.))
    {
        r = 6. * (t - (4./6.));
        g = 0.;
        b = 1.;
    }

    if(t <= (4./6.))
    {
        r = 0.;
        g = 1. - 6. * (t - (3./6.));
        b = 1.;
    }

    if(t <= (3./6.))
    {
        r = 0.;
        g = 1.;
        b = 6. * (t - (2./6.));
    }

    if(t <= (2./6.))
    {
        r = 1. - 6. * (t - (1./6.));
        g = 1.;
        b = 0.;
    }

    if(t <= (1./6.))
    {
        r = 1.;
        g = 6. * t;
    }

    return vec3(r, g, b);
}

void main()
{
    // Normalize vectors for lighting
    vec3 Normal = normalize(vN);
    vec3 Light = normalize(vL);
    vec3 Eye = normalize(vE);
    
    // Determine base color
    vec3 myColor;
    if(uUseChromaDepth != 0)
    {
        float depth = -vZ;
        
        // Debug visualization - uncomment this to see raw depth values
        myColor = vec3(depth/50.0);
        return;  // Add this line to see only the depth values
        
        // Normalize depth between 0 and 1
        float t = (depth - uRedDepth) / (uBlueDepth - uRedDepth);
        t = clamp(t, 0.0, 1.0);
        
        // Map to color range (red = close, blue = far)
        myColor = Rainbow(t);
    }
    else
    {
        // Original color/hatching logic
        float hatchFreq = 20.0;
        float hatchWidth = 0.3;
        float tolerance = 0.1;
        
        float st = vST.s * hatchFreq;
        float numLines = floor(st);
        float lineCenter = numLines + 0.5;
        float hatch = SmoothPulse(lineCenter - hatchWidth, 
                                 lineCenter + hatchWidth, 
                                 st, 
                                 tolerance);
        
        myColor = mix(vec3(0.2), vec3(0.7), hatch);
    }

    // Apply lighting
    float diff = max(dot(Normal, Light), 0.0);
    vec3 ambient = uKa * myColor;
    vec3 diffuse = uKd * diff * myColor;
    
    vec3 specular = vec3(0.0);
    if(diff > 0.0)
    {
        vec3 ref = normalize(reflect(-Light, Normal));
        float spec = pow(max(dot(Eye, ref), 0.0), uShininess);
        specular = uKs * spec * uSpecularColor.rgb;
    }
    
    vec3 finalColor = ambient + diffuse + specular;
    gl_FragColor = vec4(finalColor, 1.0);
}
