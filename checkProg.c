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
	int r = ParseFileOrString( f, 0, &obj, &numObjects );
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




	if( 0 )
	{
		struct TLEObject * o = &obj[0];

		// Let's just pick out ISS (Zarya)	
		struct elsetrec satrec;

		if( ConvertTLEToSGP4( &satrec, o ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}

		// set start/stop times for propagation, in minutes.
		double startmfe = (OGGetAbsoluteTime() - o->epoch)/60.0; // Convert to minutes.
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
	}


	{
		SGPF ro[3], vo[3];
		struct TLEObject * ss = 0;
		struct elsetrec iss;
		int numSS = 0;
		r = ParseFileOrString( 0, ""
			"ISS (ZARYA)             \n"
			"1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995\n"
			"2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710\n",
			&ss, &numSS );
		if( r || numSS != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}
		if( ConvertTLEToSGP4( &iss, &ss[0] ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}
		double startmfe = (1714240800 - ss->epoch)/60.0;
		sgp4 (&iss, startmfe, ro,  vo);

		double jd = iss.jdsatepoch + iss.jdsatepochF;
		double jdfrac = startmfe/1440.0;
		int year, mon, day, hr, min;
		double sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f /", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[%16.8f] %16.8f %16.8f %16.8f [%f] %12.9f %12.9f %12.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2]);
		/*
		  These are from:
			>>> from sgp4.api import jday
			>>> from sgp4.api import Satrec
			>>> s = "1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995"
			>>> t = "2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710"
			>>> satellite = Satrec.twoline2rv(s, t)
			>>> jd, fr = jday(2024, 4, 27, 18, 0, 0)
			>>> e, r, v = satellite.sgp4(jd, fr)
			>>> r
			(-4582.664719456509, -4356.875102968861, 2475.1474001054107)
			>>> v
			(5.036095414394779, -2.2591278380385664, 5.3188560672302145)

		  For some reason, they disagree withh these coords by a few km.
			From https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.txt
			2024-04-27T18:00:00.000 -4587.597832963610 -4337.824413854870 2498.905090951730 5.05095945730461 -2.27273088144798 5.29906134757423
			1714240800 = Sat Apr 27 2024 18:00:00 GMT+0000
		*/
		double pysgp4[3] = { -4582.664719456509, -4356.875102968861, 2475.1474001054107 };
		double pysgp4v[3] = { 5.036095414394779, -2.2591278380385664, 5.3188560672302145 };
		double rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km\n",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.05 )
		{
			printf( "PASS\n" );
		}
		else
		{
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km\n",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( rmse < 0.0005 )
		{
			printf( "PASS\n" );
		}
		else
		{
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed.  Fail\n" );
			return -6;
		}

		int iter;
		double dStartSetup = OGGetAbsoluteTime();
		for( iter = 0; iter < 1000000; iter++ )
		{
			ConvertTLEToSGP4( &iss, &ss[0] );
		}
		double dEndSetup = OGGetAbsoluteTime();
		double dStartRun = OGGetAbsoluteTime();
		for( iter = 0; iter < 1000000; iter++ )
		{
			double startmfe = (1714240800 + iter - ss->epoch)/60.0;
			sgp4 (&iss, startmfe, ro,  vo);
		}
		double dEndRun = OGGetAbsoluteTime();
		printf( "Init: %f us/iteration\n", dEndSetup - dStartSetup );
		printf( "Run:  %f us/iteration\n", dEndRun - dStartRun );
	}


	return 0;
}

