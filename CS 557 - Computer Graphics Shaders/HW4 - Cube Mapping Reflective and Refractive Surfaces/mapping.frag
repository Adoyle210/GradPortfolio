//#version 120 compatibility //using mac

uniform sampler3D	Noise3;
uniform float 		uNoiseAmp;
uniform float 		uNoiseFreq;
uniform float		uEta;
uniform float 		uMix;
uniform float 		uWhiteMix;
uniform samplerCube uReflectUnit;
uniform samplerCube uRefractUnit;

varying vec3	vNormal;
varying vec3	vEyeDir;
varying vec3	vMC;

const vec3  WHITE = vec3( 1.,1.,1. );

float angx, angy, angz; 
vec3 reflectColor;
vec4 nvx, nvy, nvz;

vec3
PerturbNormal3( float angx, float angy, float angz, vec3 n )
{
	float cx = cos( angx );
	float sx = sin( angx );
	float cy = cos( angy );
	float sy = sin( angy );
	float cz = cos( angz );
	float sz = sin( angz );
	
	// rotate about x:
	float yp =  n.y*cx - n.z*sx;	// y'
	n.z      =  n.y*sx + n.z*cx;	// z'
	n.y      =  yp;
	// n.x      =  n.x;

	// rotate about y:
	float xp =  n.x*cy + n.z*sy;	// x'
	n.z      = -n.x*sy + n.z*cy;	// z'
	n.x      =  xp;
	// n.y      =  n.y;

	// rotate about z:
	      xp =  n.x*cz - n.y*sz;	// x'
	n.y      =  n.x*sz + n.y*cz;	// y'
	n.x      = xp;
	// n.z      =  n.z;

	return normalize( n );
}


void
main( )
{
	vec3 Normal = normalize(vNormal);	// remember to unitize this
	vec3 Eye =    normalize(vEyeDir);	// remember to unitize this

	vec4 nvx = texture3D( Noise3, uNoiseFreq*vMC );
	vec4 nvy = texture3D( Noise3, uNoiseFreq*vec3(vMC.xy,vMC.z+0.33) );
	vec4 nvz = texture3D( Noise3, uNoiseFreq*vec3(vMC.xy,vMC.z+0.67) );

	float angx = nvx.r + nvx.g + nvx.b + nvx.a;	//  1. -> 3.
	angx = angx - 2.;				// -1. -> 1.
	angx *= uNoiseAmp;

	float angy = nvy.r + nvy.g + nvy.b + nvy.a;	//  1. -> 3.
	angy = angy - 2.;				// -1. -> 1.
	angy *= uNoiseAmp;

	float angz = nvz.r + nvz.g + nvz.b + nvz.a;	//  1. -> 3.
	angz = angz - 2.;				// -1. -> 1.
	angz *= uNoiseAmp;

	Normal = PerturbNormal3( angx, angy, angz, Normal );
	Normal = normalize( gl_NormalMatrix * Normal );

	vec3 reflectVector = reflect (Eye, Normal);
	vec3 reflectColor =  textureCube(uReflectUnit, reflectVector).rgb;

	vec3 refractVector = refract (Eye, Normal, uEta);

	vec3 refractColor;
	if( all( equal( refractVector, vec3(0.,0.,0.) ) ) )
	{
		refractColor = reflectColor;
	}
	else
	{
		refractColor = textureCube( uRefractUnit, refractVector ).rgb;
		refractColor = mix( refractColor, WHITE, uWhiteMix );
	}
	gl_FragColor = mix(vec4(refractColor, 1.0), vec4(reflectColor, 1.0), uMix);
}
