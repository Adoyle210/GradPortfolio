// out variables to be interpolated in the rasterizer and sent to each fragment shader:

varying  vec3  vN;	  // normal vector
varying  vec3  vL;	  // vector from point to light
varying  vec3  vE;	  // vector from point to eye
varying  vec2  vST;	  // (s,t) texture coordinates
varying  vec3  vMC;   // model coordinates

//adding uniform varibles for curtin 
uniform float uA;  //Sine wave amplitude	
uniform float uP;  //Sine wave period

const float uY0 = 1.0;
const float PI = 3.1415926535;

//adding light:
uniform float uLightX, uLightY, uLightZ;
vec3 LIGHTPOS = vec3( uLightX, uLightY, uLightZ );

void
main( )
{
	vec4 vert = gl_Vertex;

	//adding equation
	//this is a sine wave where the amplitude increases as you go down in -Y
	float z = uA * (uY0 - vert.y) * sin( 2. * PI * vert.x / uP );

	vec4 curtin = vec4 (vert.x, vert.y, z, 1.0);
	vMC = curtin.xyz;

	//added tangent vectors
	float dzdx = uA * (uY0-vert.y) * (2. * PI / uP) * cos( 2. * PI * vert.x / uP);
	float dzdy = -uA * sin( 2. * PI * vert.x / uP );

	//added vec3 tangent vectors 
	vec3 Tx = vec3(1., 0., dzdx );
	vec3 Ty = vec3(0., 1., dzdy );
	vec3 normal = normalize( cross( Tx, Ty ) );

	vST = gl_MultiTexCoord0.st;
	vec4 ECposition = gl_ModelViewMatrix * curtin;

	vN = normalize(gl_NormalMatrix * normal);  // normal vector
	vL = (gl_ModelViewMatrix * vec4(LIGHTPOS, 1.0)).xyz - ECposition.xyz;   // vector from the point to the light position
	vE = -ECposition.xyz;       // vector from the point  to the eye position

	gl_Position = gl_ModelViewProjectionMatrix * curtin;
}
