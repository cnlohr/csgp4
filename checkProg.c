#include <stdio.h>
#include "sattrack.h"
#include "os_generic.h"

#include "sgp4.h"

int main( int argc, char ** argv )
{
	if( argc != 2 )
	{
		fprintf( stderr, "Error: Usage: checkProg [.txt TLE file]\n" );
		return -6;
	}

	FILE * f = fopen( argv[1], "r" );
	if( !f || ferror( f ) )
	{
		fprintf( stderr, "Error: could not open %s\n", argv[1] );
		return -6;
	}

	struct TLEObject * obj = 0;
	int numObjects = 0;
	int r = ParseFile( f, &obj, &numObjects );
	if( r )
	{
		fprintf( stderr, "Error: Parsing failed\n" );
		return r;
	}
	printf( "Read %d objects.\n", numObjects );
	if( numObjects < 1 )
	{
		fprintf( stderr, "Did not read any objects\n" );
		return -4;
	}
	int i;
	for( i = 0; i < numObjects; i++ )
	{
		struct TLEObject * o = &obj[i];
		double diff = OGGetAbsoluteTime() - o->epoch;
		printf( "%24s %f %f %f\n", o->objectName, diff, OGGetAbsoluteTime(), o->epoch );
	}
}

