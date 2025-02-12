//#version 120 compatibility     //running on Mac

uniform float	uKa, uKd, uKs;	// coefficients of each type of lighting
uniform float	uShininess;	// specular exponent

//adding in some noise 
uniform sampler3D Noise3;
uniform float uNoiseFreq, uNoiseAmp;

//adding to global scope for pertubing the normals 
vec4 nvx, nvy;
float angx, angy;

// interpolated from the vertex shader:
varying  vec2  vST;                  // texture coords
varying  vec3  vN;                   // normal vector
varying  vec3  vL;                   // vector from point to light
varying  vec3  vE;                   // vector from point to eye
varying  vec3  vMC;			         // model coordinates

const vec3 OBJECTCOLOR          = vec3( 1., 0., 1. );           // color to make the object 
const vec3 SPECULARCOLOR        = vec3( 1., 1., 1. );


//added to perturb the normal  
vec3
PerturbNormal2( float angx, float angy, vec3 n )
{
        float cx = cos( angx );
        float sx = sin( angx );
        float cy = cos( angy );
        float sy = sin( angy );

        // rotate about x:
        float yp =  n.y*cx - n.z*sx;    // y'
        n.z      =  n.y*sx + n.z*cx;    // z'
        n.y      =  yp;
        // n.x      =  n.x;

        // rotate about y:
        float xp =  n.x*cy + n.z*sy;    // x'
        n.z      = -n.x*sy + n.z*cy;    // z'
        n.x      =  xp;
        // n.y      =  n.y;

        return normalize( n );
}

void
main( )
{
    vec3 myColor = OBJECTCOLOR;
	vec2 st = vST;

        //adding for noise lol its getting noisy in here 
        vec4 nvx = texture3D( Noise3, uNoiseFreq*vMC );
	float angx = nvx.r + nvx.g + nvx.b + nvx.a  -  2.;	// -1. to +1.
	angx *= uNoiseAmp;

        vec4 nvy = texture3D( Noise3, uNoiseFreq*vec3(vMC.xy,vMC.z+0.5) );
	float angy = nvy.r + nvy.g + nvy.b + nvy.a  -  2.;	// -1. to +1.
	angy *= uNoiseAmp;

	vec3 n = PerturbNormal2( angx, angy, vN );
        n = normalize(  gl_NormalMatrix * n  );

	// now use myColor in the per-fragment lighting equations:

        vec3 Normal    = n; //normalize(vN);
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