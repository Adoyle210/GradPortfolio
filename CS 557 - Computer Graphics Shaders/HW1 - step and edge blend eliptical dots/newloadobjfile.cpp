#include "newloadobjfile.h"


GLuint
LoadObjFile( char *name )
{
	char *cmd;		// the command string
	char *str;		// argument string

	std::vector <struct Vertex> Vertices(10000);
	std::vector <struct Normal> Normals(10000);
	std::vector <struct TextureCoord> TextureCoords(10000);

	Vertices.clear();
	Normals.clear();
	TextureCoords.clear();

	struct Vertex sv;
	struct Normal sn;
	struct TextureCoord st;

	Mtls   myMaterials;
	char * currentMtlFileName = strdup( "" );
	char * currentMtlName     = strdup( "" );
	Mtl *  currentMtl;


	// open the input file:

	FILE *fp = fopen( name, "r" );
	if( fp == NULL )
	{
		fprintf( stderr, "Cannot open .obj file '%s'\n", name );
		return 1;
	}


	float xmin = 1.e+37f;
	float ymin = 1.e+37f;
	float zmin = 1.e+37f;
	float xmax = -xmin;
	float ymax = -ymin;
	float zmax = -zmin;

	std::vector<GLuint>	DisplayLists(50);
	DisplayLists.clear( );

	GLuint currentDL = glGenLists( 1 );
	glNewList( currentDL, GL_COMPILE );
	DisplayLists.push_back( currentDL );

	glBegin( GL_TRIANGLES );

	for( ; ; )
	{
		char *line = ReadRestOfLine( fp );
		if( line == NULL )
			break;


		// skip this line if it is a comment:

		if( line[0] == '#' )
			continue;


		// skip this line if it is something we don't feel like handling today:

		if( line[0] == 'g' )
			continue;

		if( line[0] == 's' )
			continue;

		// get the command string:

		cmd = strtok( line, OBJDELIMS );


		// skip this line if it is empty:

		if( cmd == NULL )
			continue;


		if( strcmp( cmd, "mtllib" )  ==  0 )
		{
			glEnd();
			glEndList();				// end the current display list

			str = strtok( NULL, OBJDELIMS );
			delete [ ] currentMtlFileName;
			currentMtlFileName = strdup( str );
			fprintf( stderr, "Material file name = '%s'\n", currentMtlFileName );
			if(  myMaterials.Open( currentMtlFileName ) != 0 )
			{
				fprintf( stderr, "Cannot open mtl filename '%s'\n", currentMtlFileName );
				return 1;
			}
			fprintf( stderr, "Opened mtl filename '%s'\n", currentMtlFileName );

			myMaterials.ReadMtlFile( );
			myMaterials.Close( );

			GLuint currentDL = glGenLists( 1 );	// grab a new dl name
			fprintf( stderr, "currentDL = %d\n", currentDL );
			glNewList( currentDL, GL_COMPILE );	// create that new dl
			glEnable( GL_TEXTURE_2D );
			glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

			glBegin(GL_TRIANGLES);
			DisplayLists.push_back( currentDL );	// save this list in the list of lists
			continue;
		}

		if( strcmp( cmd, "usemtl" )  ==  0 )
		{
			glEnd();
			glEndList();				// end the current display list

			str = strtok( NULL, OBJDELIMS );
			delete [ ] currentMtlName;
			currentMtlName = strdup( str );
			currentMtl = myMaterials.FindMtl( currentMtlName );	// a pointer to the mtl struct

			GLuint currentDL = glGenLists( 1 );	// grab a new dl name
			fprintf( stderr, "currentDL = %d\n", currentDL );
			glNewList( currentDL, GL_COMPILE );	// create that new dl
			currentMtl->SetOpenglMtlProperties( );	// set all the opengl stuff

			glBegin(GL_TRIANGLES);
			DisplayLists.push_back( currentDL );	// save this list in the list of lists
			continue;
		}

		if( strcmp( cmd, "v" )  ==  0 )
		{
			str = strtok( NULL, OBJDELIMS );
			sv.x = atof(str);

			str = strtok( NULL, OBJDELIMS );
			sv.y = atof(str);

			str = strtok( NULL, OBJDELIMS );
			sv.z = atof(str);

			Vertices.push_back( sv );

			if( sv.x < xmin )	xmin = sv.x;
			if( sv.x > xmax )	xmax = sv.x;
			if( sv.y < ymin )	ymin = sv.y;
			if( sv.y > ymax )	ymax = sv.y;
			if( sv.z < zmin )	zmin = sv.z;
			if( sv.z > zmax )	zmax = sv.z;
			continue;
		}


		if( strcmp( cmd, "vn" )  ==  0 )
		{
			str = strtok( NULL, OBJDELIMS );
			sn.nx = atof( str );

			str = strtok( NULL, OBJDELIMS );
			sn.ny = atof( str );

			str = strtok( NULL, OBJDELIMS );
			sn.nz = atof( str );

			Normals.push_back( sn );
			continue;
		}


		if( strcmp( cmd, "vt" )  ==  0 )
		{
			st.s = st.t = st.p = 0.;

			str = strtok( NULL, OBJDELIMS );
			st.s = atof( str );

			str = strtok( NULL, OBJDELIMS );
			if( str != NULL )
				st.t = atof( str );

			str = strtok( NULL, OBJDELIMS );
			if( str != NULL )
				st.p = atof( str );

			TextureCoords.push_back( st );
			continue;
		}


		if( strcmp( cmd, "f" )  ==  0 )
		{
			struct face vertices[10];
			for( int i = 0; i < 10; i++ )
			{
				vertices[i].v = 0;
				vertices[i].n = 0;
				vertices[i].t = 0;
			}

			int sizev = (int)Vertices.size();
			int sizen = (int)Normals.size();
			int sizet = (int)TextureCoords.size();

			int numVertices = 0;
			bool valid = true;
			int vtx = 0;
			char *str;
			while( ( str = strtok( NULL, OBJDELIMS ) )  !=  NULL )
			{
				int v=0, n=0, t=0;
				ReadObjVTN( str, &v, &t, &n );

				// if v, n, or t are negative, they are wrt the end of their respective list:

				if( v < 0 )
					v += ( sizev + 1 );

				if( n < 0 )
					n += ( sizen + 1 );

				if( t < 0 )
					t += ( sizet + 1 );


				// be sure we are not out-of-bounds (<vector> will abort):

				if( t > sizet )
				{
					if( t != 0 )
						fprintf( stderr, "Read texture coord %d, but only have %d so far\n", t, sizet );
					t = 0;
				}

				if( n > sizen )
				{
					if( n != 0 )
						fprintf( stderr, "Read normal %d, but only have %d so far\n", n, sizen );
					n = 0;
				}

				if( v > sizev )
				{
					if( v != 0 )
						fprintf( stderr, "Read vertex coord %d, but only have %d so far\n", v, sizev );
					v = 0;
					valid = false;
				}

				vertices[vtx].v = v;
				vertices[vtx].n = n;
				vertices[vtx].t = t;
				vtx++;

				if( vtx >= 10 )
					break;

				numVertices++;
			}


			// if vertices are invalid, don't draw anything this time:

			if( ! valid )
				continue;

			if( numVertices < 3 )
				continue;


			// list the vertices:

			int numTriangles = numVertices - 2;

			for( int it = 0; it < numTriangles; it++ )
			{
				int vv[3];
				vv[0] = 0;
				vv[1] = it + 1;
				vv[2] = it + 2;

				// get the planar normal, in case vertex normals are not defined:

				struct Vertex *v0 = &Vertices[ vertices[ vv[0] ].v - 1 ];
				struct Vertex *v1 = &Vertices[ vertices[ vv[1] ].v - 1 ];
				struct Vertex *v2 = &Vertices[ vertices[ vv[2] ].v - 1 ];

				float v01[3], v02[3], norm[3];
				v01[0] = v1->x - v0->x;
				v01[1] = v1->y - v0->y;
				v01[2] = v1->z - v0->z;
				v02[0] = v2->x - v0->x;
				v02[1] = v2->y - v0->y;
				v02[2] = v2->z - v0->z;
				Cross( v01, v02, norm );
				Unit( norm, norm );
				glNormal3fv( norm );

				for( int vtx = 0; vtx < 3 ; vtx++ )
				{
					if( vertices[ vv[vtx] ].t != 0 )
					{
						struct TextureCoord *tp = &TextureCoords[ vertices[ vv[vtx] ].t - 1 ];
						glTexCoord2f( tp->s, tp->t );
					}

					if( vertices[ vv[vtx] ].n != 0 )
					{
						struct Normal *np = &Normals[ vertices[ vv[vtx] ].n - 1 ];
						glNormal3f( np->nx, np->ny, np->nz );
					}

					struct Vertex *vp = &Vertices[ vertices[ vv[vtx] ].v - 1 ];
					glVertex3f( vp->x, vp->y, vp->z );
				}
			}
			continue;
		}


		if( strcmp( cmd, "s" )  ==  0 )
		{
			continue;
		}

	}
	fclose( fp );

	glEnd();
	glEndList( );				// end the current display list

	// create the master dl:

	GLuint masterDL = glGenLists( 1 );
	glNewList( masterDL, GL_COMPILE );
	for( int i = 0; i < (int)DisplayLists.size( ); i++ )
	{
		glCallList( DisplayLists[i] );
	}
	glEndList( );

	fprintf( stderr, "Obj file range: [%8.3f,%8.3f,%8.3f] -> [%8.3f,%8.3f,%8.3f]\n",
		xmin, ymin, zmin,  xmax, ymax, zmax );
	fprintf( stderr, "Obj file center = (%8.3f,%8.3f,%8.3f)\n",
		(xmin+xmax)/2., (ymin+ymax)/2., (zmin+zmax)/2. );
	fprintf( stderr, "Obj file  span = (%8.3f,%8.3f,%8.3f)\n",
		xmax-xmin, ymax-ymin, zmax-zmin );

	return masterDL;
}



void
Cross( float v1[3], float v2[3], float vout[3] )
{
	float tmp[3];

	tmp[0] = v1[1]*v2[2] - v2[1]*v1[2];
	tmp[1] = v2[0]*v1[2] - v1[0]*v2[2];
	tmp[2] = v1[0]*v2[1] - v2[0]*v1[1];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}



float
Unit( float v[3] )
{
	float dist;

	dist = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		v[0] /= dist;
		v[1] /= dist;
		v[2] /= dist;
	}

	return dist;
}



float
Unit( float vin[3], float vout[3] )
{
	float dist;

	dist = vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}


char *
ReadRestOfLine( FILE *fp )
{
	static char *line;
	std::vector<char> tmp(1000);
	tmp.clear();

	for( ; ; )
	{
		int c = getc( fp );

		if( c == EOF  &&  tmp.size() == 0 )
		{
			return NULL;
		}

		if( c == EOF  ||  c == '\n' )
		{
			delete [] line;
			line = new char [ tmp.size()+1 ];
			for( int i = 0; i < (int)tmp.size(); i++ )
			{
				line[i] = tmp[i];
			}
			line[ tmp.size() ] = '\0';	// terminating null
			return line;
		}
		else
		{
			tmp.push_back( c );
		}
	}

	return (char *) "";
}


void
ReadObjVTN( char *str, int *v, int *t, int *n )
{
	// can be one of v, v//n, v/t, v/t/n:

	if( strstr( str, "//") )				// v//n
	{
		*t = 0;
		sscanf( str, "%d//%d", v, n );
		return;
	}
	else if( sscanf( str, "%d/%d/%d", v, t, n ) == 3 )	// v/t/n
	{
		return;
	}
	else
	{
		*n = 0;
		if( sscanf( str, "%d/%d", v, t ) == 2 )		// v/t
		{
			return;
		}
		else						// v
		{
			*n = *t = 0;
			sscanf( str, "%d", v );
		}
	}
}


int
Readline( FILE *fp, char *line )
{
	char *cp = line;
	int c;
	int numChars;

	while( ( c = fgetc(fp) ) != '\n'  &&  c != EOF )
	{
		if( c == '\r' )
		{
			continue;
		}
		//fprintf( stderr, "c = '%c' (0x%02x)\n", c, c );
		*cp++ = c;
	}
	*cp = '\0';

	if( c == EOF )
		return EOF;

	// reject an empty line:

	numChars = 0;
	char cc;
	for( char *cp = line; ( cc = *cp ) != '\0'; cp++ )
	{
		if( cc != '\t'  &&  cc != ' ' )
			numChars++;
	}
	if( numChars == 0 )
	{
		line[0] = '\0';
	}

	return c;	// the character that caused the while-loop to terminate
}

float *
Array3( float a, float b, float c )
{
	static float array[4];
	array[0] = a;
	array[1] = b;
	array[2] = c;
	array[3] = 1.;
	return array;
}

float *
Array3( float * abc )
{
	return Array3( abc[0], abc[1], abc[2] );
}
	


#define TEST

#ifdef TEST

int
main( int argc, char *argv[ ] )
{
	char *mtlFilename;

	if( argc == 1 )
	{
		//fprintf( stderr, "Usage: %s mtlFilename.mtl\n", argv[0] );
		//return 1;
		mtlFilename = (char *)"Lowpoly_tree_sample.mtl";
	}
	else
	{
		mtlFilename = argv[1];
	}

	Mtls myMaterials;
	if(  myMaterials.Open( mtlFilename ) != 0 )
	{
		fprintf( stderr, "Could not read the mtl file\n" );
		return 1;
	}
	myMaterials.ReadMtlFile( );
	myMaterials.Close( );

	fprintf( stderr, "\n\nTesting 'Bark'\n" );
	Mtl * barkMtl = myMaterials.FindMtl( (char *)"bark3SG" );
	if( barkMtl != (Mtl *)NULLPTR )
	{
		fprintf( stderr, "Kd = %9.5f, %9.5f, %9.5f\n",
			barkMtl->Kd[0], barkMtl->Kd[1], barkMtl->Kd[2] );
	}


	return 0;
}
#endif
