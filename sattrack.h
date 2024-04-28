#ifndef _SATTRACK_H
#define _SATTRACK_H

// Not good, but fast and tiny.
// MIT License.
// WARNING: This code is only intended for visualization purposes.
// It is NOT ACCURATE AND SHOULD NOT BE USED FOR ANY SERIOUS APPLICATION.

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "csgp4.h"

struct TLEObject
{
	int valid; // 0 for invalid, 1 for name, 2 for line 1, 4 for line 2.
	char objectName[28]; // Actually maxes at 25; [0] is the classification designator, but this pads
	uint32_t catalogNumber; 
	char internationalDesignator[12]; // Actually maxes at 8, but we pad.

	// Either-or
	double epoch; // In Unix Time.

	double epochDay;
	int epochYear;

	float meanMotion1;
	float meanMotion2; // I.e. 16003-4 = 0.000016003
	float dragTerm; // I.e. 34469-3 = 0.00034469
	int elementSetNumber;
	// Assume Ephemeris type = 0
	// Don't worry about Element Set Number
	// Check checksum.

	float inclination;
	float rightAscensionOfTheAscendingNode;
	float eccentricity; // 0115263 = 0.0115263
	float argumentOfPerigee;
	float meanAnomaly;
	double meanMotion;
	int revolutionNumberAtEpoch;
};

static float ParseFixedEponential( const char * le, int lineno, int * )
{
	int i;

	// Get rid of white space.
	while( le[0] == ' ' ) le++;

	int len = strlen( le );
	int eportion = 0;
	int iportion = 0;
	for( i = len-1; i >= 0; i-- )
	{
		char c = le[i];
		if( c == '-' || c == '+' )
		{
			eportion = atoi( le + i );
			break;
		}
	}
	int iportiondigits = i;
	if( le[0] == '-' ) iportiondigits--;
	if( le[0] == '+' || le[0] == ' ' ) { iportiondigits--; le++; }
	iportion = atoi( le );
	float ret = iportion;
	int shiftdown = iportiondigits - eportion;
	for( i = 0; i < shiftdown; i++ )
	{
		ret /= 10;
	}
	for( i = 0; i < -shiftdown; i++ )
	{
		ret *= 10;
	}
	return ret;
}


static double ConvertEpochYearAndDayToUnix( int epochYear, double epochDay )
{
	epochYear = ( epochYear > 56 ) ? ( epochYear + 1900 ) : ( epochYear + 2000 );
	int year;
	int day = 0;
	for( year = 1970; year < epochYear; year++ )
	{
		if (year % 400 == 0) {
			day += 366;
		} else if (year % 100 == 0) {
			day += 365;
		} else if (year % 4 == 0) {
			day += 366;
		} else {
			day += 365;
		}
	}

	epochDay--; // It's 1 indexed.

	return ( day + epochDay ) * 60 * 60 * 24;
}


static int ParseFileOrString( FILE * f, const char * sLineSet, struct TLEObject ** objects, int * numObjects )
{
	int i;
	ssize_t s;
	char line[256];
	char * lineptr = line;
	size_t n = sizeof( line ) - 1;
	int lineno = 0;

	int thisValid = 0;
	int aborted = 0;
	int ret = 0;
	struct TLEObject * thisObject;

	while( 1 )
	{
		if( f )
		{
			s = getline( &lineptr, &n, f );
			if( s < 0 )
				break;
		}
		else if( sLineSet )
		{
			char c;
			s = 0;
			while( ( c = *(sLineSet++) ) )			{
				line[s++] = c;
				if( c == '\n' ) break;
			}
			if( c == 0 ) break;
			line[s] = 0;
			n = s;
		}
		if( line[s-1] == '\r' || line[s-1] == '\n' ) s--;
		if( line[s-1] == '\r' || line[s-1] == '\n' ) s--;
		lineno++;
		if( thisValid == 0 )
		{
			int nObject = *numObjects;
			if( !aborted )
				++*numObjects;
			else
				ret = -5;

			aborted = 0;
			*objects = realloc( *objects, sizeof( **objects ) * *numObjects );
			thisObject = *objects + nObject;
			memset( thisObject, 0, sizeof( *thisObject ) );
			thisObject->objectName[0] = ' ';
			thisObject->objectName[1] = ' ';
			memcpy( thisObject->objectName + 2, line, 24 );

			int i;
			for( i = 24; i > 0; i-- ) 
				if( thisObject->objectName[i] == ' ' )
					thisObject->objectName[i] = 0;
				else
					break;

			thisValid = 1;
			continue;
		}

		if( line[1] != ' ' )
		{
			fprintf( stderr, "Parsing error on line %d - unexpected char at space 1\n", lineno );
			aborted = 1;
			continue; 
		}

		int checksum_check, checksum;

		switch( line[0] )
		{
		case '1':
			// 1 25544U 98067A   24108.06679608  .00019473  00000+0  34469-3 0  9999
			if( s < 68 )
			{
				fprintf( stderr, "Parsing error on line %d; too short\n", lineno );
				aborted = 1;
				continue;
			}
			checksum_check = 0;
			for( i = 0; i < s-1; i++ )
			{
				char c = line[i];
				if( c == '-' ) checksum_check++;
				if( c >= '0' && c <= '9' ) checksum_check += c - '0';
			}


			thisObject->objectName[0] = line[7];
			line[7] = 0;
			thisObject->catalogNumber = atoi( line + 2 );
			memcpy( thisObject->internationalDesignator, line + 9, 8 );

			checksum = atoi( line + 68 );
			line[68] = 0;
			if( checksum != checksum_check % 10 )
			{
				fprintf( stderr, "Checkum error on line %d.  %d != %d\n", lineno, checksum, checksum_check % 10 );
				aborted = 1;
				continue;
			}

			thisObject->elementSetNumber = atoi( line + 64 );
			line[61] = 0;
			thisObject->dragTerm = ParseFixedEponential( line + 53, lineno, &aborted );
			line[52] = 0;
			thisObject->meanMotion2 = ParseFixedEponential( line + 44, lineno, &aborted );
			line[43] = 0;
			thisObject->meanMotion1 = atof( line + 33 );
			line[32] = 0;
			// For some reason, -0.0000000005 matches the python SGP4 - but, I don't think that's right.  That said - without this delta, we sometimes seem to get squirely math on the data conversion functions.
			thisObject->epochDay = atof( line + 20 )-0.0000000005;
			line[20] = 0;
			thisObject->epochYear = atoi( line + 18 );
			line[17] = 0;

			thisObject->epoch = ConvertEpochYearAndDayToUnix( thisObject->epochYear, thisObject->epochDay );

			break;
		case '2':
			// 2 25544  51.6384 258.4693 0004986  70.1481  42.8477 15.50291806449066
			thisValid = 0; // Finalizing this record.

			if( s < 68 )
			{
				fprintf( stderr, "Parsing error on line %d; too short\n", lineno );
				aborted = 1;
				continue;
			}
			checksum_check = 0;
			for( i = 0; i < s-1; i++ )
			{
				char c = line[i];
				if( c == '-' ) checksum_check++;
				if( c >= '0' && c <= '9' ) checksum_check += c - '0';
			}

			int catalogNumber = atoi( line + 2 );
			if( thisObject->catalogNumber != catalogNumber )
			{
				fprintf( stderr, "Error: Mismatching catalog numbers at line %d\n", lineno );
				aborted = 1;
			}

			thisObject->objectName[0] = line[7];
			memcpy( thisObject->internationalDesignator, line + 9, 8 );

			checksum = atoi( line + 68 );
			line[68] = 0;
			if( checksum != checksum_check % 10 )
			{
				fprintf( stderr, "Checkum error on line %d.  %d != %d\n", lineno, checksum, checksum_check % 10 );
				aborted = 1;
				continue;
			}

			thisObject->revolutionNumberAtEpoch = atoi( line + 63 );
			line[63] = 0;
			thisObject->meanMotion = atof( line + 52 );
			line[51] = 0;
			thisObject->meanAnomaly = atof( line + 43 );
			line[42] = 0;
			thisObject->argumentOfPerigee = atof( line + 34 );
			line[33] = 0;

			thisObject->eccentricity = atof( line + 26 );
			if( line[26] == '-' ) thisObject->eccentricity *= 0.000001;
			else                  thisObject->eccentricity *= 0.0000001;
			line[25] = 0;

			thisObject->rightAscensionOfTheAscendingNode = atof( line + 17 );
			line[16] = 0;
			thisObject->inclination = atof( line + 8 );
			line[7] = 0;
			thisObject->valid = 1;

			break;
		default:
			fprintf( stderr, "Unknown record type %c on line %d\n", line[0], lineno );
		}
	}

	// If last element is complete, drop it.
	if( *numObjects && thisValid ) --*numObjects;

	return ret;
}



static int ConvertTLEToSGP4( struct elsetrec * satrec, struct TLEObject * obj )
{ 
	// Test from TestSGP4.cpp

	// Example
	//ISS (ZARYA)             
	//1 25544U 98067A   24108.06679608  .00019473  00000+0  34469-3 0  9999
	//2 25544  51.6384 258.4693 0004986  70.1481  42.8477 15.50291806449066

	const SGPF deg2rad  =   SGPPI / 180.0;         //   0.0174532925199433
	const SGPF xpdotp   =  1440.0 / (2.0 *SGPPI);  // 229.1831180523293

	if( !obj->valid ) return -1;

	// Options are: enum gravconsttype { wgs72old, wgs72, wgs84 }; // wgs72 is the standard and should be used with JSPOC TLEs
	enum gravconsttype whichconst = wgs72;
	char opsmode = 'a';
	snprintf( satrec->satnum, sizeof( satrec->satnum ), "%d", obj->catalogNumber );

	// ----------------------------------------------------------------
	// find sgp4epoch time of element set
	// remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
	// and minutes from the epoch (time)
	// ----------------------------------------------------------------

	// ---------------- temp fix for years from 1957-2056 -------------------
	// --------- correct fix will occur when year is 4-digit in tle ---------
	satrec->epochyr = obj->epochYear;
	satrec->epochdays = obj->epochDay;

	int mon, day, hr, minute, year;
	SGPF sec;
	if (satrec->epochyr < 57)
		year = satrec->epochyr + 2000;
	else
		year = satrec->epochyr + 1900;
	days2mdhms(year, satrec->epochdays, &mon, &day, &hr, &minute, &sec);
	jday(year, mon, day, hr, minute, sec, &satrec->jdsatepoch, &satrec->jdsatepochF);

	satrec->no_kozai = obj->meanMotion; //2.00491383;
	satrec->ecco = obj->eccentricity; //0.6877146;
	satrec->inclo = obj->inclination; //64.1586;
	satrec->nodeo = obj->rightAscensionOfTheAscendingNode; //279.0717;
	satrec->argpo = obj->argumentOfPerigee; //264.7651;
	satrec->mo = obj->meanAnomaly;          //20.2257;
	satrec->nddot = obj->meanMotion2;       //0.00000e0;
	satrec->bstar = obj->dragTerm;          //0.11873e-3;
	satrec->ndot = obj->meanMotion1;        //0.00000099;
	satrec->elnum = obj->elementSetNumber;  //813;
	satrec->revnum = obj->revolutionNumberAtEpoch; //22565;
	satrec->classification = obj->objectName[0];
	strncpy(satrec->intldesg, &obj->objectName[1], sizeof(satrec->intldesg) );

	satrec->ephtype = 0;

	// convert units and initialize
	satrec->no_kozai = satrec->no_kozai / xpdotp; //* rad/min
	satrec->ndot = satrec->ndot  / (xpdotp*1440.0);  //* ? * minperday
	satrec->nddot= satrec->nddot / (xpdotp*1440.0*1440);
	satrec->inclo = satrec->inclo  * deg2rad;
	satrec->nodeo = satrec->nodeo  * deg2rad;
	satrec->argpo = satrec->argpo  * deg2rad;
	satrec->mo    = satrec->mo     * deg2rad;

	// sgp4fix not needed here
	// satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
	// satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

	sgp4init( whichconst, opsmode, satrec->satnum, satrec->jdsatepoch-2433281.5 /* ???!!?? */, satrec->bstar,
		 satrec->ndot, satrec->nddot, satrec->ecco, satrec->argpo, satrec->inclo, satrec->mo, satrec->no_kozai,
		 satrec->nodeo, satrec );

	return 0;
}

#endif

