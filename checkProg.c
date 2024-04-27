#include <stdio.h>
#include "os_generic.h"
#include "sattrack.h"

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
		//double diff = OGGetAbsoluteTime() - o->epoch;
		printf( "%24s %f\n", o->objectName, o->epoch );//, diff, OGGetAbsoluteTime(), o->epoch );
	}


	struct TLEObject * o = &obj[0];

	// Let's just pick out ISS (Zarya)	
	struct elsetrec satrec;

	if( ConvertTLEToSGP4( &satrec, o ) )
	{
		printf( "failed to convert tle\n" );
		return -5;
	}

	// set start/stop times for propagation, in minutes.
	double startmfe = 1440.0;//(OGGetAbsoluteTime() - o->epoch)/60.0; // Convert to minutes.
	double stopmfe  = startmfe + 45.0;
	double deltamin = 1.0;

	double tsince = startmfe;
	while ((tsince < stopmfe) && (satrec.error == 0))
	{
		SGPF ro[3], vo[3];

		if(tsince > stopmfe)
			tsince = stopmfe;

		// .25us per call
		sgp4 (&satrec, tsince, ro,  vo);
#if 0
		double jd = satrec.jdsatepoch + satrec.jdsatepochF;
		double jdfrac = tsince/1440.0;
		int year, mon, day, hr, min;
		double sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f /", jd, jdfrac, year, mon, day, hr, min, sec );
#endif
		printf( "%16.8f %16.8f %16.8f %16.8f [%f] %12.9f %12.9f %12.9f\n",
			tsince,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2]);

		tsince = tsince + deltamin;
	}

#if 0
	// Test from TestSGP4.cpp

	// sgp4fix demonstrate method of running SGP4 directly from orbital element values
	//1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
	//2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656
	const double deg2rad  =   SGPPI / 180.0;         //   0.0174532925199433
	const double xpdotp   =  1440.0 / (2.0 *SGPPI);  // 229.1831180523293

	struct elsetrec satrec;

	enum gravconsttype whichconst = wgs72;
	char opsmode = 'a';
	strcpy( satrec.satnum, "8195" );
	satrec.jdsatepoch = 2453911.8321544402;
	satrec.no_kozai = 2.00491383;
	satrec.ecco = 0.6877146;
	satrec.inclo = 64.1586;
	satrec.nodeo = 279.0717;
	satrec.argpo = 264.7651;
	satrec.mo = 20.2257;
	satrec.nddot = 0.00000e0;
	satrec.bstar = 0.11873e-3;
	satrec.ndot = 0.00000099;
	satrec.elnum = 813;
	satrec.revnum = 22565;
	satrec.classification = 'U';
	strcpy(satrec.intldesg, "        ");
	satrec.ephtype = 0;

	// convert units and initialize
	satrec.no_kozai = satrec.no_kozai / xpdotp; //* rad/min
	satrec.ndot = satrec.ndot  / (xpdotp*1440.0);  //* ? * minperday
	satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);
	satrec.inclo = satrec.inclo  * deg2rad;
	satrec.nodeo = satrec.nodeo  * deg2rad;
	satrec.argpo = satrec.argpo  * deg2rad;
	satrec.mo    = satrec.mo     * deg2rad;

	// set start/stop times for propagation
	double startmfe =     0.0;
	double stopmfe  =  2880.0;
	double deltamin =   120.0;

	sgp4init( whichconst, opsmode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar,
		 satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
		 satrec.nodeo, &satrec);

	double tsince = startmfe;
	while ((tsince < stopmfe) && (satrec.error == 0))
	{
		double ro[3], vo[3];
		tsince = tsince + deltamin;

		if(tsince > stopmfe)
			tsince = stopmfe;

		sgp4 (&satrec,  tsince, ro,  vo);

		double jd = satrec.jdsatepoch + tsince/1440.0;
		//invjday( jd, year,mon,day,hr,min, sec );

		printf( " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
			tsince,ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);
	}
#endif

}

