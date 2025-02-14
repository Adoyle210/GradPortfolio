//#version 120 compatibility // using mac 

varying vec3	vNormal;
varying vec3	vEyeDir;
varying vec3	vMC;


void
main( )
{    
	vMC = gl_Vertex.xyz;
	vec3 ECposition = ( gl_ModelViewMatrix * gl_Vertex ).xyz;
	vEyeDir = ECposition.xyz - vec3( 0., 0., 0. ) ; 
	       		// vector from the eye position to the point
	vNormal = normalize( gl_Normal );

	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
