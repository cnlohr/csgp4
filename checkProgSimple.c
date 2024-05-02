#include <stdio.h>
#include "os_generic.h"

#include "csgp4_simple.h"

#include "csgp4.h"

CSGP4_DECORATOR int ConvertTLEToSGP4GPU( struct TLEObject * obj, SGPF initial_time, SGPF * initial_r, SGPF* initial_v, SGPF * a_alta_altp )
{
	// can use 'a' or 'i' methods.
	// *	epoch	   - epoch time in days from jan 0, 1950. 0 hr
	// But, jdsatepoch is in days from 4713 bc
	sgp4init_simple( obj->jdsatepoch-2433281.5 + obj->jdsatepochF, obj->dragTerm,
		 obj->meanMotion1, obj->meanMotion2, obj->eccentricity, obj->argumentOfPerigee, obj->inclination, obj->meanAnomaly, obj->meanMotion,
		 obj->rightAscensionOfTheAscendingNode, initial_time, initial_r, initial_v, a_alta_altp );
	return 0;
}


int main( int argc, char ** argv )
{
	puts( "Simple Tests");

	double pysgp4[3] = { -4582.664719456509, -4356.875102968861, 2475.1474001054107 };
	double pysgp4v[3] = { 5.036095414394779, -2.2591278380385664, 5.3188560672302145 };
	SGPF ro[3], vo[3], a_alta_altp[3];
	double rmse;

	if( 1 )
	{
		struct TLEObject * ss = 0;
		int numSS = 0;
		int r = ParseFileOrString( 0, ""
			"ISS (ZARYA)             \n"
			"1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995\n"
			"2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710\n",
			&ss, &numSS );
		if( r || numSS != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}
		double startmfe = (1714240800 - ss->epoch)/60.0;
		if( ConvertTLEToSGP4GPU( &ss[0], startmfe, ro,  vo, a_alta_altp ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}
		puts( ss->objectName );

		SGPF jd = ss->jdsatepoch + floor(startmfe/1440.0);
		SGPF jdfrac = fmod(startmfe/1440.0, 1.0) + ss->jdsatepochF;
		int year, mon, day, hr, min;
		SGPF sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f\n", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[Δt%14.8f] %16.8f %16.8f %16.8f %16.8f \n                   %16.9f %16.9f %16.9f %16.9f\n                   %16.9f %16.9f %16.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2], sqrt(vo[0]*vo[0]+vo[1]*vo[1]+vo[2]*vo[2]),
			a_alta_altp[0], a_alta_altp[1], a_alta_altp[2]
		);

		/*
		  These are from the following python code:

from sgp4.api import jday
from sgp4.api import Satrec
s = "1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995"
t = "2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710"
satellite = Satrec.twoline2rv(s, t)
jd, fr = jday(2024, 4, 27, 18, 0, 0)
e, r, v = satellite.sgp4(jd, fr)
	>>> r
	(-4582.664719456509, -4356.875102968861, 2475.1474001054107)
	>>> v
	(5.036095414394779, -2.2591278380385664, 5.3188560672302145)

		  For some reason, they disagree withh these coords by a few km.
			From https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.txt
			2024-04-27T18:00:00.000 -4587.597832963610 -4337.824413854870 2498.905090951730 5.05095945730461 -2.27273088144798 5.29906134757423
			1714240800 = Sat Apr 27 2024 18:00:00 GMT+0000
		*/
		pysgp4[0] = -4582.664719456509;
		pysgp4[1] = -4356.875102968861;
		pysgp4[2] = 2475.1474001054107;
		pysgp4v[0] = 5.036095414394779;
		pysgp4v[1] = -2.2591278380385664;
		pysgp4v[2] = 5.3188560672302145;
		rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km ... ",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.005 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km/s ... ",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( vrmse < 0.00005 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed.  Fail\n" );
			return -6;
		}

/* And forwarded out a day... 
jd, fr = jday(2024, 4, 28, 18, 0, 0)
e, r, v = satellite.sgp4(jd, fr)
	>>> r
	(-5357.955394050937, -2614.4147285937256, 3236.6540120050795)
	>>> v
	(0.04217831766517118, -6.0043223779983865, -4.77009742984791)
*/
	}


	if( 1 )
	{
		struct TLEObject * ss = 0;
		int numSS = 0;
		int r = ParseFileOrString( 0, ""
			"ISS (ZARYA)             \n"
			"1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995\n"
			"2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710\n",
			&ss, &numSS );
		if( r || numSS != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}
		puts( ss->objectName );
		double startmfe = (1714327200 - ss->epoch)/60.0;
		if( ConvertTLEToSGP4GPU( &ss[0], startmfe, ro,  vo, 0 ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}

		double jd = ss->jdsatepoch + floor(startmfe/1440.0);
		double jdfrac = fmod(startmfe/1440.0, 1.0) + ss->jdsatepochF;
		int year, mon, day, hr, min;
		SGPF sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f\n", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[Δt%14.8f] %16.8f %16.8f %16.8f %16.8f \n                   %16.9f %16.9f %16.9f %16.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2], sqrt(vo[0]*vo[0]+vo[1]*vo[1]+vo[2]*vo[2]));

		pysgp4[0] = 4435.874337686209; 
		pysgp4[1] = 4191.117631955163;
		pysgp4[2] = -2991.3331751931737;
		pysgp4v[0] = -5.401088744185615;
		pysgp4v[1] = 2.177125902892223;
		pysgp4v[2] = -4.972896867609246;
		rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km ... ",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.005 + CSGP4_USE_FLOAT * .05 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km/s ... ",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( vrmse < 0.00005 + CSGP4_USE_FLOAT * .001 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed (%f).  Fail\n", vrmse );
			return -6;
		}
	}


// Forward a day, but, with Arase, a higher orbit satellite, ALSO TEST wth only init.
/*
from sgp4.api import jday
from sgp4.api import Satrec
s = "1 41896U 16080A   24120.85435253  .00000361  00000+0  77632-3 0  9994"
t = "2 41896  31.9360 317.6417 7003302 291.8988  10.4729  2.55468020 68280"
satellite = Satrec.twoline2rv(s, t)
jd, fr = jday(2024, 5, 2, 18, 0, 0)
e, r, v = satellite.sgp4(jd, fr)
	>>> r
	(13105.747459057653, 29909.425386662235, 19138.1545898884)
	>>> v
	(-1.3853033635228822, 1.227640540353242, -0.039183565331240844)
*/
	if( 1 )
	{
		struct TLEObject * ss_arase = 0;
		int numSS_arase = 0;
		int r = ParseFileOrString( 0, ""
			"ARASE (ERG)             \n"
			"1 41896U 16080A   24120.85435253  .00000361  00000+0  77632-3 0  9994\n"
			"2 41896  31.9360 317.6417 7003302 291.8988  10.4729  2.55468020 68280\n",
			&ss_arase, &numSS_arase );
		if( r || numSS_arase != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}

		double startmfe = (1714672800 - ss_arase->epoch)/60.0;
		if( ConvertTLEToSGP4GPU( &ss_arase[0], startmfe, ro,  vo, 0 ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}
		puts( ss_arase->objectName );

		double jd = ss_arase->jdsatepoch + floor(startmfe/1440.0);
		double jdfrac = fmod(startmfe/1440.0, 1.0) + ss_arase->jdsatepochF;
		int year, mon, day, hr, min;
		SGPF sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f\n", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[Δt%14.8f] %16.8f %16.8f %16.8f %16.8f \n                   %16.9f %16.9f %16.9f %16.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2], sqrt(vo[0]*vo[0]+vo[1]*vo[1]+vo[2]*vo[2]));

		pysgp4[0] = 13105.747459057653; 
		pysgp4[1] = 29909.425386662235;
		pysgp4[2] = 19138.1545898884;
		pysgp4v[0] = -1.3853033635228822;
		pysgp4v[1] = 1.227640540353242;
		pysgp4v[2] = -0.039183565331240844;
		rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km ... ",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.0005 + CSGP4_USE_FLOAT * .03 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km/s ... ",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( vrmse < 0.00005 + CSGP4_USE_FLOAT * .002 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed (%f).  Fail\n", vrmse );
			return -6;
		}
	}

// Forward a day, but, with THEMIS, a much higher orbit satellite
/*
from sgp4.api import jday
from sgp4.api import Satrec
s = "1 30798U 07004E   24121.56859868 -.00000520  00000+0  00000+0 0  9999"
t = "2 30798   7.9414 283.2653 8322119 302.3633 341.0652  0.87832507 59071"
satellite = Satrec.twoline2rv(s, t)
jd, fr = jday(2024, 5, 2, 18, 0, 0)
e, r, v = satellite.sgp4(jd, fr)
	>>> r
	(13105.747459057653, 29909.425386662235, 19138.1545898884)
	>>> v
	(-1.3853033635228822, 1.227640540353242, -0.039183565331240844)
*/
	{
		struct TLEObject * ss_THEMIS = 0;
		int numSS_THEMIS = 0;
		int r = ParseFileOrString( 0, ""
			"THEMIS E                \n"
			"1 30798U 07004E   24121.56859868 -.00000520  00000+0  00000+0 0  9999\n"
			"2 30798   7.9414 283.2653 8322119 302.3633 341.0652  0.87832507 59071\n",
			&ss_THEMIS, &numSS_THEMIS );
		if( r || numSS_THEMIS != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}
		puts( ss_THEMIS->objectName );
		double startmfe = (1714672800 - ss_THEMIS->epoch)/60.0;
		if( ConvertTLEToSGP4GPU( &ss_THEMIS[0], startmfe, ro,  vo, 0 ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}

		//Dumpelsetrec( &iss_THEMIS );

		double jd = ss_THEMIS->jdsatepoch + floor(startmfe/1440.0);
		double jdfrac = fmod(startmfe/1440.0, 1.0) + ss_THEMIS->jdsatepochF;
		int year, mon, day, hr, min;
		SGPF sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f\n", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[Δt%14.8f] %16.8f %16.8f %16.8f %16.8f \n                   %16.9f %16.9f %16.9f %16.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2], sqrt(vo[0]*vo[0]+vo[1]*vo[1]+vo[2]*vo[2]));

		pysgp4[0] = 12197.874919034643;
		pysgp4[1] = 48838.6479258657;
		pysgp4[2] = 3134.8786450225966;
		pysgp4v[0] = -1.96146751028316;
		pysgp4v[1] = -1.7885483032900817;
		pysgp4v[2] = -0.32311386956142407;
		rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km ... ",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.0005 + CSGP4_USE_FLOAT * .055 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km/s ... ",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( vrmse < 0.00005 + CSGP4_USE_FLOAT * .001 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed (%f).  Fail\n", vrmse );
			return -6;
		}
	}


/* And forwarded out a month... 
jd, fr = jday(2024, 5, 27, 18, 0, 0)
e, r, v = satellite.sgp4(jd, fr)
	>>> r
	(5099.551520031815, 1808.2683576301836, -4104.365753671076)
	>>> v
	(0.7417244009740825, 6.593295105250736, 3.825415802504736)
*/
	if( 1 )
	{
		struct TLEObject * ss = 0;
		int numSS = 0;
		int r = ParseFileOrString( 0, ""
			"ISS (ZARYA)             \n"
			"1 25544U 98067A   24118.69784154  .00029521  00000+0  51116-3 0  9995\n"
			"2 25544  51.6397 205.7509 0003603 115.2045 341.2688 15.50662375450710\n",
			&ss, &numSS );
		if( r || numSS != 1 )
		{
			fprintf( stderr, "Error: SS can't load.\n" );
			return -5;
		}

		double startmfe = (1716832800 - ss->epoch)/60.0;
		if( ConvertTLEToSGP4GPU( &ss[0], startmfe, ro,  vo, 0 ) )
		{
			printf( "failed to convert tle\n" );
			return -5;
		}
		puts( ss->objectName );
		double jd = ss->jdsatepoch + floor(startmfe/1440.0);
		double jdfrac = fmod(startmfe/1440.0, 1.0) + ss->jdsatepochF;
		int year, mon, day, hr, min;
		SGPF sec;
		invjday( jd, jdfrac, &year, &mon, &day, &hr, &min, &sec );
		printf( "%f %f %04d %02d %02d %02d:%02d:%05.02f\n", jd, jdfrac, year, mon, day, hr, min, sec );
		printf( "[Δt%14.8f] %16.8f %16.8f %16.8f %16.8f \n                   %16.9f %16.9f %16.9f %16.9f\n",
			startmfe,ro[0],ro[1],ro[2], sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]),vo[0],vo[1],vo[2], sqrt(vo[0]*vo[0]+vo[1]*vo[1]+vo[2]*vo[2]));

		pysgp4[0] = 5099.551520031815;
		pysgp4[1] = 1808.2683576301836;
		pysgp4[2] = -4104.365753671076;
		pysgp4v[0] = 0.7417244009740825;
		pysgp4v[1] = 6.593295105250736;
		pysgp4v[2] = 3.825415802504736;
		rmse = sqrt( (ro[0] - pysgp4[0])*(ro[0] - pysgp4[0]) + (ro[1] - pysgp4[1]) * (ro[1] - pysgp4[1]) + ( ro[2] - pysgp4[2] ) * (ro[2] - pysgp4[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km ... ",
			 ro[0] - pysgp4[0], ro[1] - pysgp4[1], ro[2] - pysgp4[2],
			rmse );
		if( rmse < 0.005 + CSGP4_USE_FLOAT * 5 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in position.  Fail\n" );
			return -5;
		}

		double vrmse = sqrt( (vo[0] - pysgp4v[0])*(vo[0] - pysgp4v[0]) + (vo[1] - pysgp4v[1]) * (vo[1] - pysgp4v[1]) + ( vo[2] - pysgp4v[2] ) * (vo[2] - pysgp4v[2]) );
		printf( "Python / Our SGP4 Disagreement: %f %f %f RMS: %f km/s ... ",
			 vo[0] - pysgp4v[0], vo[1] - pysgp4v[1], vo[2] - pysgp4v[2],
			vrmse );
		if( vrmse < 0.00005 + CSGP4_USE_FLOAT * .006 )
		{
			printf( "PASS\n" );
		}
		else
		{
			printf( "FAIL\n" );
			fprintf( stderr, "Error: SGP Algorithm disagrees in speed (%f).  Fail\n", vrmse );
			return -6;
		}

	}










	// Perf test
	{
		struct TLEObject * ss = 0;
		int numss = 0;
		int iter = 0;
		numss = 0;
		int r = ParseFileOrString( 0, ""
			"THEMIS E                \n"
			"1 30798U 07004E   24121.56859868 -.00000520  00000+0  00000+0 0  9999\n"
			"2 30798   7.9414 283.2653 8322119 302.3633 341.0652  0.87832507 59071\n",
			&ss, &numss );
		if( r || numss != 1 )
		{
			fprintf( stderr, "Error: THEMIS can't load.\n" );
			return -5;
		}
		double dStartFullInit = OGGetAbsoluteTime();
		for( iter = 0; iter < 1000000; iter++ )
		{
			double startmfe = (1714240800 + iter - ss->epoch)/60.0;
			ConvertTLEToSGP4GPU( &ss[0], startmfe, ro, vo, 0 );
		}
		double dEndFullInit = OGGetAbsoluteTime();
		printf( "Deep Space Full At Init: %.4f us/iteration\n", dEndFullInit - dStartFullInit );



	}
}
