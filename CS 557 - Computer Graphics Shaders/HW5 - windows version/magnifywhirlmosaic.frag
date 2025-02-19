#version 330 compatibility

in vec2 vST;

uniform float uSc;
uniform float uTc;
uniform float uRad;
uniform float uMag;
uniform float uWhirl;
uniform float uMosaic;
uniform sampler2D uImageUnit;

void main()
{
    vec2 st = vST - vec2(uSc, uTc);  // make (0,0) the center of the circle
    float dist = length(st);

    if (dist > uRad)
    {
        vec3 rgb = texture(uImageUnit, vST).rgb;
        gl_FragColor = vec4(rgb, 1.);
    }
    else
    {
        // Magnifying
        float r = dist;
        float rPrime = r / uMag;

        // Whirling
        float theta = atan(st.t, st.s);
        float thetaPrime = theta - uWhirl * rPrime;

        // Restoring (s,t)
        st = rPrime * vec2(cos(thetaPrime), sin(thetaPrime));
        st += vec2(uSc, uTc);

        // Mosaic'ing
        int numins = int(1.0 / uMosaic);
        int numint = int(1.0 / uMosaic);
        float sc = floor(st.s * numins) / numins + 0.5 / numins;
        float tc = floor(st.t * numint) / numint + 0.5 / numint;
        st.s = sc;
        st.t = tc;

        vec3 rgb = texture(uImageUnit, st).rgb;
        gl_FragColor = vec4(rgb, 1.);
    }
}
