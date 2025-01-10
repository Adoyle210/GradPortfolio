#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

const void  *	NULLPTR = (void *)0;
const char *	DELIMS = " \t";


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
	

class
Mtl
{
    public:
	char *		Name;
	float		Ka[3];
	float		Kd[3];
	float		Ke[3];
	float		Ks[3];
	float		D;
	float		Ni;
	float		Ns;
	float		Illum;
	char *		MapKd;
	Mtl *		Next;

	Mtl( )
	{
		Next = (Mtl *) NULLPTR;
	};

};


class
Mtls
{
    public:
	FILE *		Fp;
	Mtl *		First;
	Mtl *		Last;

	Mtls( )
	{
		Fp = (FILE *) NULLPTR;
		First = Last = (Mtl *) NULLPTR;
	};

	int			// 0 = success, !0 = failure
	Open( char * fileName )
	{
		Fp = fopen( fileName, "r" );
		if( Fp == NULLPTR )
		{
			fprintf( stderr, "Cannot open file '%s'\n", fileName );
			return 1;
		}

		return 0;
	};

	void
	ReadMtlFile( )
	{
		char line[256];
		Mtl * thisMtl = (Mtl *) NULLPTR;
		int c;
		while( ( c = Readline( Fp, line) ) != EOF )
		{
			if( c == EOF )
				break;

			if( line[0] == '\0' )
				continue;

			//fprintf( stderr, "line = '%s'\n", line );

			char * tok = strtok( line, DELIMS );			// first token is the variable name
			char cc;
			for( char *cp = tok; ( cc = *cp ) != '\0'; cp++ )
			{
				if( isupper(cc) )
					*cp = tolower( cc );
			}
			//fprintf( stderr, "Variable name = '%s'\n", tok );

			// process the rest of the line:

			while( tok != (char *) NULLPTR )
			{
				if( tok[0] == '#' )
					break;		// comment causes us to stop reading this line
				//fprintf( stderr, "tok = *%s*\n", tok );

				if( strcmp( tok, "newmtl" ) == 0 )
				{
					thisMtl = new Mtl;
					thisMtl->Next = (Mtl *) NULLPTR;
					if( First == (Mtl *) NULLPTR )
					{
						First = thisMtl;
					}
					if( Last != (Mtl *) NULLPTR )
					{
						Last->Next = thisMtl;
					}
					Last = thisMtl;
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Name = strdup(tok);
					fprintf( stderr, "Material name = '%s'\n", thisMtl->Name );
					break;
				}

				if( strcmp( tok, "ka" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Ka[0] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Ka[1] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Ka[2] = atof( tok );
					break;
				}

				if( strcmp( tok, "kd" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Kd[0] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Kd[1] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					if( tok != (char *)NULLPTR )	thisMtl->Kd[2] = atof( tok );
					break;
				}

				if( strcmp( tok, "ke" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ke[0] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ke[1] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ke[2] = atof( tok );
					break;
				}

				if( strcmp( tok, "ks" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ks[0] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ks[1] = atof( tok );
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ks[2] = atof( tok );
					break;
				}

				if( strcmp( tok, "d" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->D = atof( tok );
					break;
				}

				if( strcmp( tok, "ni" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ni = atof( tok );
					break;
				}

				if( strcmp( tok, "ns" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Ns = atof( tok );
					break;
				}

				if( strcmp( tok, "illum" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->Illum = atof( tok );
					break;
				}

				if( strcmp( tok, "map_kd" ) == 0 )
				{
					tok = strtok( (char *)NULLPTR, DELIMS );
					thisMtl->MapKd = strdup( tok );
					break;
				}

				fprintf( stderr, "Don't recognize Mtl file token '%s'\n", tok );
			}
		}
	};

	void
	Close( )
	{
		fclose( Fp );
		Fp = (FILE *)NULLPTR;
	};

	Mtl *
	FindMtl( char *mtlName )
	{
		for( Mtl * mtl = First; mtl != (Mtl *)NULLPTR; mtl = mtl->Next )
		{
			if( strcmp( mtl->Name, mtlName ) == 0 )
			{
				return mtl;
			}
		}
		return (Mtl *)NULLPTR;
	};

};


#define TEST_THIS

#ifdef  TEST_THIS

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
	Mtl * barkMtl = myMaterials.FindMtl( (char *)"Bark" );
	if( barkMtl != (Mtl *)NULLPTR )
	{
		fprintf( stderr, "Kd = %9.5f, %9.5f, %9.5f\n",
			barkMtl->Kd[0], barkMtl->Kd[1], barkMtl->Kd[2] );
	}

	fprintf( stderr, "\n\nTesting 'Tree'\n" );
	Mtl * treeMtl = myMaterials.FindMtl( (char *)"Tree" );
	if( treeMtl != (Mtl *)NULLPTR )
	{
		fprintf( stderr, "Kd = %9.5f, %9.5f, %9.5f\n",
			treeMtl->Kd[0], treeMtl->Kd[1], treeMtl->Kd[2] );
	}


	return 0;
}
#endif
