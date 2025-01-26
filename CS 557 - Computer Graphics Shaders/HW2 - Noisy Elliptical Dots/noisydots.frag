//#version 120 compatibility     //running on Mac
// you can set these 4 uniform variables dynamically or hardwire them:

uniform float	uKa, uKd, uKs;	// coefficients of each type of lighting
uniform float	uShininess;	// specular exponent

// in Project #1, these have to be set dynamically from keytime animations or by keyboard hits:
uniform float	uAd, uBd;
uniform float	uTol;

//added for noise (slide24)
uniform sampler3D Noise3;
uniform float uNoiseFreq, uNoiseAmp;
float ellipseEQ, nv, newDist, scale;


// interpolated from the vertex shader:
varying  vec2  vST;                  // texture coords
varying  vec3  vN;                   // normal vector
varying  vec3  vL;                   // vector from point to light
varying  vec3  vE;                   // vector from point to eye
varying  vec3  vMC;			         // model coordinates


const vec3 OBJECTCOLOR          = vec3( 1., 0., 1. );           // color to make the object
const vec3 ELLIPSECOLOR         = vec3( 0., 1., 1. );           // color to make the ellipse
const vec3 SPECULARCOLOR        = vec3( 1., 1., 1. );

void
main( )
{
    vec3 myColor = OBJECTCOLOR;
	vec2 st = vST;

    //adding from index noise form 3d model coords
    //vec4 nv = texture3D (Noise3, uNoiseFreq * vec3(vST,0.));
    vec4 nv = texture3D (Noise3, uNoiseFreq * vMC);
    float n = nv.r + nv.g + nv.b + nv.a;    //one to three 
    n = n - 2.;                              //positve one to negative one 
    n *= uNoiseAmp; 

	// blend OBJECTCOLOR and ELLIPSECOLOR by using the ellipse equation to decide how close
	// 	this fragment is to the ellipse border:

    //For ellipse equation: 
    float Ar = uAd / 2.;
    float Br = uBd / 2.;
    int numins = int( st.s / uAd );
	int numint = int( st.t / uBd );
    float sc = float(numins) * uAd + Ar;
    float tc = float(numint) * uBd + Br;

    //adding for ellipse noise:
    float ds = st.s - sc;                           //wrt elipse center 
    float dt = st.t - tc;                           //wrt elipse center 
    float oldDist = sqrt( ds*ds + dt*dt);
    float newDist = oldDist + n;                   
    float scale = newDist / oldDist;                 //this could be < 1., = 1. or > 1. 

    ds *= scale;         // scale by noise factor
    ds /= Ar;            //ellipse equation 
    dt *= scale;         //scale by noise factor 
    dt /= Br;             //ellipse equation
    ellipseEQ = ds*ds + dt*dt;

    //ellipse equation 
	//float ellipse = (((st.s - sc) * (st.s - sc)) / (Ar * Ar)) + (((st.t - tc) * (st.t - tc)) / (Br * Br));
	
    float t = smoothstep( 1.-uTol, 1.+uTol, ellipseEQ );
        myColor = mix( ELLIPSECOLOR, OBJECTCOLOR, t );
    gl_FragColor = vec4(myColor,1.);

	// now use myColor in the per-fragment lighting equations:

        vec3 Normal    = normalize(vN);
        vec3 Light     = normalize(vL);
        vec3 Eye       = normalize(vE);

        vec3 ambient = uKa * myColor;

        float d = max( dot(Normal,Light), 0. );       // only do diffuse if the light can see the point
        vec3 diffuse = uKd * d * myColor;

        float s = 0.;
        if( d > 0. )              // only do specular if the light can see the point
        {
                vec3 ref = normalize(  reflect( -Light, Normal )  );
                float cosphi = dot( Eye, ref );
                if( cosphi > 0. )
                        s = pow( max( cosphi, 0. ), uShininess );
        }
        vec3 specular = uKs * s * SPECULARCOLOR.rgb;
        gl_FragColor = vec4( ambient + diffuse + specular,  1. );

}