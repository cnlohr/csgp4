// A shader-compilable, C header transliterated from David Vallado's SGP4Lib.cs, version: "SGP4 Version 2020-03-12";
// https://celestrak.org/software/vallado-sw.php
// I could not find a license, but assume whatever license the original code is under.
// I (Charles Lohr) just ported the code back to C from C#. This code is much less featureful than the original code
// AND ABSOLUTELY DOES NOT MAINTAIN THE SAME LEVEL OF ACCURACY.
//
// USE AT YOUR OWN RISK
//
// This file gets rid of basically all optional features of this system, and only does supports 'a' computation mode.
//
// Please note: There are still branches based off static initialization, please be sensitive to organize your
// satellites in such a way you don't accidentally parallelize multiple satellites together.

#ifndef _CSGP4_SIMPLE_H
#define _CSGP4_SIMPLE_H

#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CSGP4_OUT * restrict
#define CSGP4_REF(x) (*x)

#define SGPPI 3.1415926535897932384626433

// Default to double calculations.
#ifndef CSGP4_USE_FLOAT
#define CSGP4_USE_FLOAT 0
#endif

#ifndef CSGP4_INIT
#define CSGP4_INIT 1
#endif


#ifndef CSGP4_DECORATOR
#define CSGP4_DECORATOR static
#endif




#if CSGP4_USE_FLOAT
// Use single precision.
#define SGPF float
#define FLOOR floorf
#define FABS  fabsf
#define COS   cosf
#define SIN   sinf
#define ATAN2 atan2f
#define SQRT  sqrtf
#define POW   powf
#define FMOD  fmodf
#else
// Use double precision.
#define SGPF double
#define FLOOR floor
#define FABS  fabs
#define COS   cos
#define SIN   sin
#define ATAN2 atan2
#define SQRT  sqrt
#define POW   pow
#define FMOD fmod
#endif


#define BOOL int


struct elsetrec_simple
{
	int error;
	bool method; // true for 'd' false for 'n' (Deep / Near Eart)

	/* Near Earth */
	int isimp;

	SGPF aycof, con41, cc1, cc4, cc5, d2, d3, d4,
		   delmo, eta, argpdot, omgcof, sinmao, t, t2cof, t3cof,
		   t4cof, t5cof, x1mth2, x7thm1, mdot, nodedot, xlcof, xmcof, nodecf;


	SGPF  bstar,  inclo, nodeo, ecco, argpo, mo;

	/* Deep Space */
	int irez;
	SGPF d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232,
		   d5421, d5433, dedt, del1, del2, del3, didt, dmdt,
		   dnodt, domdt, e3, ee2, peo, pgho, pho, pinco,
		   plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3,
		   si2, si3, sl2, sl3, sl4, gsto, xfact, xgh2,
		   xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3,
		   xl4, xlamo, zmol, zmos, atime, xli, xni;
	// sgp4fix add unkozai'd variable
	SGPF no_unkozai;
	// sgp4fix add singly averaged variables
	SGPF am, em, im, Om, om, mm, nm;
	// sgp4fix add constant parameters to eliminate mutliple calls during execution
	SGPF tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;


	SGPF a /* Not actually used in algo but fun to look at */;
	SGPF no_kozai,rcse;
	SGPF ndot /* Not actually used in algo*/;
	SGPF nddot /* Not actually used in algo*/;
	SGPF alta /* Not used in algo, but cool anyway to look at */;
	SGPF altp /* Not used in algo, but cool anyway to look at */;
};

CSGP4_DECORATOR void dspace_simple
	 (
	   int irez,
	   SGPF d2201, SGPF d2211, SGPF d3210, SGPF d3222, SGPF d4410,
	   SGPF d4422, SGPF d5220, SGPF d5232, SGPF d5421, SGPF d5433,
	   SGPF dedt, SGPF del1, SGPF del2, SGPF del3, SGPF didt,
	   SGPF dmdt, SGPF dnodt, SGPF domdt, SGPF argpo, SGPF argpdot,
	   SGPF t, SGPF tc, SGPF gsto, SGPF xfact, SGPF xlamo,
	   SGPF no,
	   SGPF CSGP4_OUT atime, SGPF CSGP4_OUT em, SGPF CSGP4_OUT argpm, SGPF CSGP4_OUT inclm, SGPF CSGP4_OUT xli,
	   SGPF CSGP4_OUT mm, SGPF CSGP4_OUT xni, SGPF CSGP4_OUT nodem, SGPF CSGP4_OUT dndt, SGPF CSGP4_OUT nm
	 )
{
	const SGPF twopi = 2.0 * SGPPI;
	int iretn;  //, iret;
	SGPF delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, g22, g32,
		 g44, g52, g54, fasx2, fasx4, fasx6, rptim, step2, stepn, stepp;

	fasx2 = 0.13130908;
	fasx4 = 2.8843198;
	fasx6 = 0.37448087;
	g22 = 5.7686396;
	g32 = 0.95240898;
	g44 = 1.8014998;
	g52 = 1.0508330;
	g54 = 4.4108898;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	stepp = 720.0;
	stepn = -720.0;
	step2 = 259200.0;

	/* ----------- calculate deep space resonance effects ----------- */
	CSGP4_REF(dndt) = 0.0;
	theta = fmod( (gsto + tc * rptim), twopi );
	CSGP4_REF(em) = CSGP4_REF(em) + dedt * t;

	CSGP4_REF(inclm) = CSGP4_REF(inclm) + didt * t;
	CSGP4_REF(argpm) = CSGP4_REF(argpm) + domdt * t;
	CSGP4_REF(nodem) = CSGP4_REF(nodem) + dnodt * t;
	CSGP4_REF(mm) = CSGP4_REF(mm) + dmdt * t;

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//	inclm = -inclm;
	//	argpm = argpm - Math.PI;
	//	nodem = nodem + Math.PI;
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted
	// sgp4fix set initial values for c#
	xndt = 0.0;
	xnddt = 0.0;
	xldot = 0.0;

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft = 0.0;
	if (irez != 0)
	{
		// sgp4fix streamline check
		if ((CSGP4_REF(atime) == 0.0) || (t * CSGP4_REF(atime) <= 0.0) || (FABS(t) < FABS(CSGP4_REF(atime))))
		{
			CSGP4_REF(atime) = 0.0;
			CSGP4_REF(xni) = no;
			CSGP4_REF(xli) = xlamo;
		}
		SGPF xli_cache = CSGP4_REF(xli);
		// sgp4fix move check outside loop
		if (t > 0.0)
			delt = stepp;
		else
			delt = stepn;

		//iretn = 381; // added for do loop
		//iret = 0; // added for loop
		//while (iretn == 381)
		bool iretn = true;
		do
		{
			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if (irez != 2)
			{
				xndt = del1 * SIN(xli_cache - fasx2) + del2 * SIN(2.0 * (xli_cache - fasx4)) +
						del3 * SIN(3.0 * (xli_cache - fasx6));
				xldot = *xni + xfact;
				xnddt = del1 * COS(xli_cache - fasx2) +
						2.0 * del2 * COS(2.0 * (xli_cache - fasx4)) +
						3.0 * del3 * COS(3.0 * (xli_cache - fasx6));
				xnddt = xnddt * xldot;
			}
			else
			{
				/* --------- near - half-day resonance terms -------- */
				xomi = argpo + argpdot * *atime;
				x2omi = xomi + xomi;
				x2li = xli_cache + xli_cache;
				xndt = d2201 * SIN(x2omi + xli_cache - g22) + d2211 * SIN(xli_cache - g22) +
					  d3210 * SIN(xomi + xli_cache - g32) + d3222 * SIN(-xomi + xli_cache - g32) +
					  d4410 * SIN(x2omi + x2li - g44) + d4422 * SIN(x2li - g44) +
					  d5220 * SIN(xomi + xli_cache- g52) + d5232 * SIN(-xomi + *xli - g52) +
					  d5421 * SIN(xomi + x2li - g54) + d5433 * SIN(-xomi + x2li - g54);
				xldot = *xni + xfact;
				xnddt = d2201 * COS(x2omi + xli_cache - g22) + d2211 * COS(xli_cache - g22) +
					  d3210 * COS(xomi + xli_cache - g32) + d3222 * COS(-xomi + xli_cache - g32) +
					  d5220 * COS(xomi + xli_cache - g52) + d5232 * COS(-xomi + xli_cache - g52) +
					  2.0 * (d4410 * COS(x2omi + x2li - g44) +
					  d4422 * COS(x2li - g44) + d5421 * COS(xomi + x2li - g54) +
					  d5433 * COS(-xomi + x2li - g54));
				xnddt = xnddt * xldot;
			}

			/* ----------------------- integrator ------------------- */
			// sgp4fix move end checks to end of routine
			if (FABS(t - CSGP4_REF(atime)) >= stepp)
			{
				//iret = 0;
				iretn = true;
			}
			else // exit here
			{
				ft = t - CSGP4_REF(atime);
				iretn = false;
			}

			if (iretn)
			{
				xli_cache = xli_cache + xldot * delt + xndt * step2;
				CSGP4_REF(xni) = CSGP4_REF(xni) + xndt * delt + xnddt * step2;
				CSGP4_REF(atime) = CSGP4_REF(atime) + delt;
			}
		} while( iretn );

		CSGP4_REF(nm) = CSGP4_REF(xni) + xndt * ft + xnddt * ft * ft * 0.5;
		xl = xli_cache + xldot * ft + xndt * ft * ft * 0.5;
		CSGP4_REF(xli) = xli_cache;
		if (irez != 1)
		{
			CSGP4_REF(mm) = xl - 2.0 * CSGP4_REF(nodem) + 2.0 * theta;
			CSGP4_REF(dndt) = CSGP4_REF(nm) - no;
		}
		else
		{
			CSGP4_REF(mm) = xl - CSGP4_REF(nodem) - CSGP4_REF(argpm) + theta;
			CSGP4_REF(dndt) = CSGP4_REF(nm) - no;
		}
		CSGP4_REF(nm) = no + CSGP4_REF(dndt);
	}

	//#include "debug4.cpp"
}  // end dsspace

CSGP4_DECORATOR void dpper_simple
	 (
	   SGPF e3, SGPF ee2, SGPF peo, SGPF pgho, SGPF pho,
	   SGPF pinco, SGPF plo, SGPF se2, SGPF se3, SGPF sgh2,
	   SGPF sgh3, SGPF sgh4, SGPF sh2, SGPF sh3, SGPF si2,
	   SGPF si3, SGPF sl2, SGPF sl3, SGPF sl4, SGPF t,
	   SGPF xgh2, SGPF xgh3, SGPF xgh4, SGPF xh2, SGPF xh3,
	   SGPF xi2, SGPF xi3, SGPF xl2, SGPF xl3, SGPF xl4,
	   SGPF zmol, SGPF zmos, SGPF inclo,
	   BOOL init,
	   SGPF CSGP4_OUT ep, SGPF CSGP4_OUT inclp, SGPF CSGP4_OUT nodep, SGPF CSGP4_OUT argpp, SGPF CSGP4_OUT mp
	 )
{
	/* --------------------- local variables ------------------------ */
	const SGPF twopi = 2.0 * SGPPI;
	SGPF alfdp, betdp, cosip, cosop, dalf, dbet, dls,
		 f2, f3, pe, pgh, ph, pinc, pl,
		 sel, ses, sghl, sghs, shll, shs, sil,
		 sinip, sinop, sinzf, sis, sll, sls, xls,
		 xnoh, zf, zm, zel, zes, znl, zns;

	/* ---------------------- constants ----------------------------- */
	zns = 1.19459e-5;
	zes = 0.01675;
	znl = 1.5835218e-4;
	zel = 0.05490;

	/* --------------- calculate time varying periodics ----------- */
	zm = zmos + zns * t;
	// be sure that the initial call has time set to zero
	if (init)
		zm = zmos;
	zf = zm + 2.0 * zes * SIN(zm);
	sinzf = SIN(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * COS(zf);
	ses = se2 * f2 + se3 * f3;
	sis = si2 * f2 + si3 * f3;
	sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
	sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
	shs = sh2 * f2 + sh3 * f3;
	zm = zmol + znl * t;
	if (init)
		zm = zmol;
	zf = zm + 2.0 * zel * SIN(zm);
	sinzf = SIN(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * COS(zf);
	sel = ee2 * f2 + e3 * f3;
	sil = xi2 * f2 + xi3 * f3;
	sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
	sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
	shll = xh2 * f2 + xh3 * f3;
	pe = ses + sel;
	pinc = sis + sil;
	pl = sls + sll;
	pgh = sghs + sghl;
	ph = shs + shll;

	if (!init)
	{
		pe = pe - peo;
		pinc = pinc - pinco;
		pl = pl - plo;
		pgh = pgh - pgho;
		ph = ph - pho;
		CSGP4_REF(inclp) = CSGP4_REF(inclp) + pinc;
		CSGP4_REF(ep) = CSGP4_REF(ep) + pe;
		sinip = SIN(CSGP4_REF(inclp));
		cosip = COS(CSGP4_REF(inclp));

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if (*inclp >= 0.2)
		{
			ph = ph / sinip;
			pgh = pgh - cosip * ph;
			CSGP4_REF(argpp) = CSGP4_REF(argpp) + pgh;
			CSGP4_REF(nodep) = CSGP4_REF(nodep) + ph;
			CSGP4_REF(mp) = CSGP4_REF(mp) + pl;
		}
		else
		{
			/* ---- apply periodics with lyddane modification ---- */
			sinop = SIN(CSGP4_REF(nodep));
			cosop = COS(CSGP4_REF(nodep));
			alfdp = sinip * sinop;
			betdp = sinip * cosop;
			dalf = ph * cosop + pinc * cosip * sinop;
			dbet = -ph * sinop + pinc * cosip * cosop;
			alfdp = alfdp + dalf;
			betdp = betdp + dbet;
			CSGP4_REF(nodep) = fmod( CSGP4_REF(nodep), twopi );
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if ((CSGP4_REF(nodep) < 0.0)) // && (opsmode == 'a'))
				CSGP4_REF(nodep) = CSGP4_REF(nodep) + twopi;
			xls = CSGP4_REF(mp) + CSGP4_REF(argpp) + cosip * CSGP4_REF(nodep);
			dls = pl + pgh - pinc * CSGP4_REF(nodep) * sinip;
			xls = xls + dls;
			xnoh = CSGP4_REF(nodep);
			CSGP4_REF(nodep) = ATAN2(alfdp, betdp);
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if ((CSGP4_REF(nodep) < 0.0)) // && (opsmode == 'a'))
				CSGP4_REF(nodep) = CSGP4_REF(nodep) + twopi;
			if (FABS(xnoh - CSGP4_REF(nodep)) > SGPPI)
				if (CSGP4_REF(nodep) < xnoh)
					CSGP4_REF(nodep) = CSGP4_REF(nodep) + twopi;
				else
					CSGP4_REF(nodep) = CSGP4_REF(nodep) - twopi;
			CSGP4_REF(mp) = CSGP4_REF(mp) + pl;
			CSGP4_REF(argpp) = xls - CSGP4_REF(mp) - cosip * CSGP4_REF(nodep);
		}
	}   // if init == 'n'
	//#include "debug1.cpp"
}  // end dpper


CSGP4_DECORATOR void sgp4_simple
	 (
	   struct elsetrec_simple * satrec, SGPF tsince,
	   SGPF r[3], SGPF v[3]
	 )
{
	SGPF am, axnl, aynl, betal, cosim, cnod,
		cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu,
		delm, delomg, em, emsq, ecose, el2, eo1,
		ep, esine, argpm, argpp, argpdf, pl, mrt = 0.0,
		mvt, rdotl, rl, rvdot, rvdotl, sinim,
		sin2u, sineo1, sini, sinip, sinsu, sinu,
		snod, su, t2, t3, t4, tem5, temp,
		temp1, temp2, tempa, tempe, templ, u, ux,
		uy, uz, vx, vy, vz, inclm, mm,
		nm, nodem, xinc, xincp, xl, xlm, mp,
		xmdf, xmx, xmy, nodedf, xnode, nodep, tc, dndt,
		twopi, x2o3,  //, j2, j3, tumin, j4, xke, j3oj2, radiusearthkm,
		vkmpersec, delmtemp;   // mu, 
	int ktr;

	// assign initial values
	r[0] = 0.0;
	r[1] = 0.0;
	r[2] = 0.0;
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;

	/* ------------------ set mathematical constants --------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(Math.PI-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	const SGPF temp4 = 1.5e-12;
	twopi = 2.0 * SGPPI;
	x2o3 = 2.0 / 3.0;
	// sgp4fix identify constants and allow alternate values
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	vkmpersec = satrec->radiusearthkm * satrec->xke / 60.0;

	/* --------------------- clear sgp4 error flag ----------------- */
	satrec->t = tsince;
	satrec->error = 0;

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf = satrec->mo + satrec->mdot * satrec->t;
	argpdf = satrec->argpo + satrec->argpdot * satrec->t;
	nodedf = satrec->nodeo + satrec->nodedot * satrec->t;
	argpm = argpdf;
	mm = xmdf;
	t2 = satrec->t * satrec->t;
	nodem = nodedf + satrec->nodecf * t2;
	tempa = 1.0 - satrec->cc1 * satrec->t;
	tempe = satrec->bstar * satrec->cc4 * satrec->t;
	templ = satrec->t2cof * t2;

	if (satrec->isimp != 1)
	{
		delomg = satrec->omgcof * satrec->t;
		// sgp4fix use mutliply for speed instead of Math.POW
		delmtemp = 1.0 + satrec->eta * COS(xmdf);
		delm = satrec->xmcof *
				 (delmtemp * delmtemp * delmtemp -
				 satrec->delmo);
		temp = delomg + delm;
		mm = xmdf + temp;
		argpm = argpdf - temp;
		t3 = t2 * satrec->t;
		t4 = t3 * satrec->t;
		tempa = tempa - satrec->d2 * t2 - satrec->d3 * t3 -
						 satrec->d4 * t4;
		tempe = tempe + satrec->bstar * satrec->cc5 * (SIN(mm) -
						 satrec->sinmao);
		templ = templ + satrec->t3cof * t3 + t4 * (satrec->t4cof +
						 satrec->t * satrec->t5cof);
	}

	nm = satrec->no_unkozai;
	em = satrec->ecco;
	inclm = satrec->inclo;

	if (satrec->method)
	{
		tc = satrec->t;
		dspace_simple
			(
			  satrec->irez,
			  satrec->d2201, satrec->d2211, satrec->d3210,
			  satrec->d3222, satrec->d4410, satrec->d4422,
			  satrec->d5220, satrec->d5232, satrec->d5421,
			  satrec->d5433, satrec->dedt, satrec->del1,
			  satrec->del2, satrec->del3, satrec->didt,
			  satrec->dmdt, satrec->dnodt, satrec->domdt,
			  satrec->argpo, satrec->argpdot, satrec->t, tc,
			  satrec->gsto, satrec->xfact, satrec->xlamo,
			  satrec->no_unkozai, &satrec->atime,
			  &em, &argpm, &inclm, &satrec->xli, &mm, &satrec->xni,
			  &nodem, &dndt, &nm
			);
	} // if method = true

	if (nm <= 0.0)
	{
		//		 printf("# error nm %f\n", nm);
		satrec->error = 2;
		// sgp4fix add return
		//return false;
	}

	am = POW((satrec->xke / nm), x2o3) * tempa * tempa;
	nm = satrec->xke / POW(am, 1.5);
	em = em - tempe;

	// fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/ )
	{
		//		 printf("# error em %f\n", em);
		satrec->error = 1;
		// sgp4fix to return if there is an error in eccentricity
		//return false;
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if (em < 1.0e-6)
		em = 1.0e-6;
	mm = mm + satrec->no_unkozai * templ;
	xlm = mm + argpm + nodem;
	emsq = em * em;
	temp = 1.0 - emsq;

	nodem = fmod((nodem), twopi );
	argpm = fmod((argpm), twopi );
	xlm = fmod((xlm), twopi );
	mm = fmod((xlm - argpm - nodem), twopi);

	// sgp4fix recover singly averaged mean elements
	satrec->am = am;
	satrec->em = em;
	satrec->im = inclm;
	satrec->Om = nodem;
	satrec->om = argpm;
	satrec->mm = mm;
	satrec->nm = nm;

	/* ----------------- compute extra mean quantities ------------- */
	sinim = SIN(inclm);
	cosim = COS(inclm);

	/* -------------------- add lunar-solar periodics -------------- */
	ep = em;
	xincp = inclm;
	argpp = argpm;
	nodep = nodem;
	mp = mm;
	sinip = sinim;
	cosip = cosim;
	if (satrec->method)
	{
		dpper_simple
			(
			  satrec->e3, satrec->ee2, satrec->peo,
			  satrec->pgho, satrec->pho, satrec->pinco,
			  satrec->plo, satrec->se2, satrec->se3,
			  satrec->sgh2, satrec->sgh3, satrec->sgh4,
			  satrec->sh2, satrec->sh3, satrec->si2,
			  satrec->si3, satrec->sl2, satrec->sl3,
			  satrec->sl4, satrec->t, satrec->xgh2,
			  satrec->xgh3, satrec->xgh4, satrec->xh2,
			  satrec->xh3, satrec->xi2, satrec->xi3,
			  satrec->xl2, satrec->xl3, satrec->xl4,
			  satrec->zmol, satrec->zmos, satrec->inclo,
			  false, &ep, &xincp, &nodep, &argpp, &mp
			);
		if (xincp < 0.0)
		{
			xincp = -xincp;
			nodep = nodep + SGPPI;
			argpp = argpp - SGPPI;
		}
		if ((ep < 0.0) || (ep > 1.0))
		{
			//			printf("# error ep %f\n", ep);
			satrec->error = 3;
			// sgp4fix add return
			//return false;
		}
	} // if method = true

	/* -------------------- long period periodics ------------------ */
	if (satrec->method)
	{
		sinip = SIN(xincp);
		cosip = COS(xincp);
		satrec->aycof = -0.5 * satrec->j3oj2 * sinip;
		// sgp4fix for divide by zero for xincp = 180 deg
		if (FABS(cosip + 1.0) > 1.5e-12)
			satrec->xlcof = -0.25 * satrec->j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
		else
			satrec->xlcof = -0.25 * satrec->j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
	}
	axnl = ep * COS(argpp);
	temp = 1.0 / (am * (1.0 - ep * ep));
	aynl = ep * SIN(argpp) + temp * satrec->aycof;
	xl = mp + argpp + nodep + temp * satrec->xlcof * axnl;

	/* --------------------- solve kepler's equation --------------- */
	u = fmod( (xl - nodep), twopi );
	eo1 = u;
	tem5 = 9999.9;
	ktr = 1;
	// sgp4fix for c# intiialize
	coseo1 = 0.0;
	sineo1 = 0.0;
	//   sgp4fix for kepler iteration
	//   the following iteration needs better limits on corrections
	while ((FABS(tem5) >= 1.0e-12) && (ktr <= 10))
	{
		sineo1 = SIN(eo1);
		coseo1 = COS(eo1);
		tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
		tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
		if (FABS(tem5) >= 0.95)
			tem5 = tem5 > 0.0 ? 0.95 : -0.95;
		eo1 = eo1 + tem5;
		ktr = ktr + 1;
	}

	/* ------------- short period preliminary quantities ----------- */
	ecose = axnl * coseo1 + aynl * sineo1;
	esine = axnl * sineo1 - aynl * coseo1;
	el2 = axnl * axnl + aynl * aynl;
	pl = am * (1.0 - el2);
	if (pl < 0.0)
	{
		//		 printf("# error pl %f\n", pl);
		satrec->error = 4;
		// sgp4fix add return
		//return false;
	}
	else
	{
		rl = am * (1.0 - ecose);
		rdotl = SQRT(am) * esine / rl;
		rvdotl = SQRT(pl) / rl;
		betal = SQRT(1.0 - el2);
		temp = esine / (1.0 + betal);
		sinu = am / rl * (sineo1 - aynl - axnl * temp);
		cosu = am / rl * (coseo1 - axnl + aynl * temp);
		su = ATAN2(sinu, cosu);
		sin2u = (cosu + cosu) * sinu;
		cos2u = 1.0 - 2.0 * sinu * sinu;
		temp = 1.0 / pl;
		temp1 = 0.5 * satrec->j2 * temp;
		temp2 = temp1 * temp;

		/* -------------- update for short period periodics ------------ */
		if (satrec->method  ) // == 'd')
		{
			cosisq = cosip * cosip;
			satrec->con41 = 3.0 * cosisq - 1.0;
			satrec->x1mth2 = 1.0 - cosisq;
			satrec->x7thm1 = 7.0 * cosisq - 1.0;
		}
		mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec->con41) +
				0.5 * temp1 * satrec->x1mth2 * cos2u;
		su = su - 0.25 * temp2 * satrec->x7thm1 * sin2u;
		xnode = nodep + 1.5 * temp2 * cosip * sin2u;
		xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
		mvt = rdotl - nm * temp1 * satrec->x1mth2 * sin2u / satrec->xke;
		rvdot = rvdotl + nm * temp1 * (satrec->x1mth2 * cos2u +
				1.5 * satrec->con41) / satrec->xke;

		/* --------------------- orientation vectors ------------------- */
		sinsu = SIN(su);
		cossu = COS(su);
		snod = SIN(xnode);
		cnod = COS(xnode);
		sini = SIN(xinc);
		cosi = COS(xinc);
		xmx = -snod * cosi;
		xmy = cnod * cosi;
		ux = xmx * sinsu + cnod * cossu;
		uy = xmy * sinsu + snod * cossu;
		uz = sini * sinsu;
		vx = xmx * cossu - cnod * sinsu;
		vy = xmy * cossu - snod * sinsu;
		vz = sini * cossu;

		/* --------- position and velocity (in km and km/sec) ---------- */
		r[0] = (mrt * ux) * satrec->radiusearthkm;
		r[1] = (mrt * uy) * satrec->radiusearthkm;
		r[2] = (mrt * uz) * satrec->radiusearthkm;
		v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
		v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
		v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
	}  // if pl > 0

	// sgp4fix for decaying satellites
	if (mrt < 1.0)
	{
		//		 printf("# decay condition %11.6f \n",mrt);
		satrec->error = 6;
		//return false;
	}

	//#include "debug7.cpp"
	//return true;
}  // end sgp4


CSGP4_DECORATOR void initl_simple
	 (
	// sgp4fix satn not needed. include in satrec in case needed later  
	// int satn,	  
	// sgp4fix just pass in xke and j2
	// gravconsttype whichconst, 
	   SGPF xke, SGPF j2,
	   SGPF ecco, SGPF epoch, SGPF inclo, SGPF no_kozai, bool CSGP4_OUT method,

	// Output params:
	   SGPF CSGP4_OUT ainv, SGPF CSGP4_OUT ao, SGPF CSGP4_OUT con41, SGPF CSGP4_OUT con42, SGPF CSGP4_OUT cosio,
	   SGPF CSGP4_OUT cosio2, SGPF CSGP4_OUT eccsq, SGPF CSGP4_OUT omeosq, SGPF CSGP4_OUT posq,
	   SGPF CSGP4_OUT rp, SGPF CSGP4_OUT rteosq, SGPF CSGP4_OUT sinio, SGPF CSGP4_OUT gsto,
	   SGPF CSGP4_OUT no_unkozai
	 )
{
	/* --------------------- local variables ------------------------ */
	SGPF ak, d1, del, adel, po, x2o3;

	// sgp4fix use old way of finding gst
	SGPF ds70;
	SGPF ts70, tfrac, c1, thgr70, fk5r, c1p2p;
	const SGPF twopi = 2.0 * SGPPI;

	/* ----------------------- earth constants ---------------------- */
	// sgp4fix identify constants and allow alternate values
	// only xke and j2 are used here so pass them in directly
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	x2o3 = 2.0 / 3.0;

	/* ------------- calculate auxillary epoch quantities ---------- */
	*eccsq = ecco * ecco;
	*omeosq = 1.0 - *eccsq;
	*rteosq = SQRT(*omeosq);
	*cosio = COS(inclo);
	*cosio2 = *cosio * *cosio;

	/* ------------------ un-kozai the mean motion ----------------- */
	ak = POW(xke / no_kozai, x2o3);
	d1 = 0.75 * j2 * (3.0 * *cosio2 - 1.0) / (*rteosq * *omeosq);
	del = d1 / (ak * ak);
	adel = ak * (1.0 - del * del - del *
			(1.0 / 3.0 + 134.0 * del * del / 81.0));
	del = d1 / (adel * adel);
	CSGP4_REF(no_unkozai) = no_kozai / (1.0 + del);

	CSGP4_REF(ao) = POW(xke / (CSGP4_REF(no_unkozai)), x2o3);
	CSGP4_REF(sinio) = SIN(inclo);
	po = CSGP4_REF(ao) * CSGP4_REF(omeosq);
	CSGP4_REF(con42) = 1.0 - 5.0 * CSGP4_REF(cosio2);
	CSGP4_REF(con41) = -CSGP4_REF(con42) - CSGP4_REF(cosio2) - CSGP4_REF(cosio2);
	CSGP4_REF(ainv) = 1.0 / CSGP4_REF(ao);
	CSGP4_REF(posq) = po * po;
	CSGP4_REF(rp) = CSGP4_REF(ao) * (1.0 - ecco);
	CSGP4_REF(method) = false;
	//*method = 'n';

	// sgp4fix modern approach to finding sidereal time
	// ALWAYS USE A OPS MODE in this form.
	//if (opsmode == 'a')
	{
		// sgp4fix use old way of finding gst
		// count integer number of days from 0 jan 1970
		ts70 = epoch - 7305.0;
		ds70 = FLOOR(ts70 + 1.0e-8);
		tfrac = ts70 - ds70;
		// find greenwich location at epoch
		c1 = 1.72027916940703639e-2;
		thgr70 = 1.7321343856509374;
		fk5r = 5.07551419432269442e-15;
		c1p2p = c1 + twopi;
		CSGP4_REF(gsto) = fmod(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
		if (CSGP4_REF(gsto) < 0.0)
			CSGP4_REF(gsto) = CSGP4_REF(gsto) + twopi;
	}
	//else
	//	*gsto = gstime(epoch + 2433281.5);

	//#include "debug5.cpp"
}  // end initl


CSGP4_DECORATOR void dscom_simple
	 (
	   SGPF epoch, SGPF ep, SGPF argpp, SGPF tc, SGPF inclp,
	   SGPF nodep, SGPF np,
	   SGPF CSGP4_OUT snodm, SGPF CSGP4_OUT cnodm, SGPF CSGP4_OUT sinim, SGPF CSGP4_OUT cosim, SGPF CSGP4_OUT sinomm,
	   SGPF CSGP4_OUT cosomm, SGPF CSGP4_OUT day, SGPF CSGP4_OUT e3, SGPF CSGP4_OUT ee2, SGPF CSGP4_OUT em,
	   SGPF CSGP4_OUT emsq, SGPF CSGP4_OUT gam, SGPF CSGP4_OUT peo, SGPF CSGP4_OUT pgho, SGPF CSGP4_OUT pho,
	   SGPF CSGP4_OUT pinco, SGPF CSGP4_OUT plo, SGPF CSGP4_OUT rtemsq, SGPF CSGP4_OUT se2, SGPF CSGP4_OUT se3,
	   SGPF CSGP4_OUT sgh2, SGPF CSGP4_OUT sgh3, SGPF CSGP4_OUT sgh4, SGPF CSGP4_OUT sh2, SGPF CSGP4_OUT sh3,
	   SGPF CSGP4_OUT si2, SGPF CSGP4_OUT si3, SGPF CSGP4_OUT sl2, SGPF CSGP4_OUT sl3, SGPF CSGP4_OUT sl4,
	   SGPF CSGP4_OUT s1, SGPF CSGP4_OUT s2, SGPF CSGP4_OUT s3, SGPF CSGP4_OUT s4, SGPF CSGP4_OUT s5,
	   SGPF CSGP4_OUT s6, SGPF CSGP4_OUT s7, SGPF CSGP4_OUT ss1, SGPF CSGP4_OUT ss2, SGPF CSGP4_OUT ss3,
	   SGPF CSGP4_OUT ss4, SGPF CSGP4_OUT ss5, SGPF CSGP4_OUT ss6, SGPF CSGP4_OUT ss7, SGPF CSGP4_OUT sz1,
	   SGPF CSGP4_OUT sz2, SGPF CSGP4_OUT sz3, SGPF CSGP4_OUT sz11, SGPF CSGP4_OUT sz12, SGPF CSGP4_OUT sz13,
	   SGPF CSGP4_OUT sz21, SGPF CSGP4_OUT sz22, SGPF CSGP4_OUT sz23, SGPF CSGP4_OUT sz31, SGPF CSGP4_OUT sz32,
	   SGPF CSGP4_OUT sz33, SGPF CSGP4_OUT xgh2, SGPF CSGP4_OUT xgh3, SGPF CSGP4_OUT xgh4, SGPF CSGP4_OUT xh2,
	   SGPF CSGP4_OUT xh3, SGPF CSGP4_OUT xi2, SGPF CSGP4_OUT xi3, SGPF CSGP4_OUT xl2, SGPF CSGP4_OUT xl3,
	   SGPF CSGP4_OUT xl4, SGPF CSGP4_OUT nm, SGPF CSGP4_OUT z1, SGPF CSGP4_OUT z2, SGPF CSGP4_OUT z3,
	   SGPF CSGP4_OUT z11, SGPF CSGP4_OUT z12, SGPF CSGP4_OUT z13, SGPF CSGP4_OUT z21, SGPF CSGP4_OUT z22,
	   SGPF CSGP4_OUT z23, SGPF CSGP4_OUT z31, SGPF CSGP4_OUT z32, SGPF CSGP4_OUT z33, SGPF CSGP4_OUT zmol,
	   SGPF CSGP4_OUT zmos
	 )
{
	/* -------------------------- constants ------------------------- */
	const SGPF zes = 0.01675;
	const SGPF zel = 0.05490;
	const SGPF c1ss = 2.9864797e-6;
	const SGPF c1l = 4.7968065e-7;
	const SGPF zsinis = 0.39785416;
	const SGPF zcosis = 0.91744867;
	const SGPF zcosgs = 0.1945905;
	const SGPF zsings = -0.98088458;
	const SGPF twopi = 2.0 * SGPPI;

	/* --------------------- local variables ------------------------ */
	int lsflg;
	SGPF a1, a2, a3, a4, a5, a6, a7,
	   a8, a9, a10, betasq, cc, ctem, stem,
	   x1, x2, x3, x4, x5, x6, x7,
	   x8, xnodce, xnoi, zcosg, zcosgl, zcosh, zcoshl,
	   zcosi, zcosil, zsing, zsingl, zsinh, zsinhl, zsini,
	   zsinil, zx, zy;

	// sgp4fix - initialize the parameters for c#
	CSGP4_REF(ss1) = 0.0;
	CSGP4_REF(ss2) = 0.0;
	CSGP4_REF(ss3) = 0.0;
	CSGP4_REF(ss4) = 0.0;
	CSGP4_REF(ss5) = 0.0;
	CSGP4_REF(ss6) = 0.0;
	CSGP4_REF(ss7) = 0.0;
	CSGP4_REF(s1) = 0.0;
	CSGP4_REF(s2) = 0.0;
	CSGP4_REF(s3) = 0.0;
	CSGP4_REF(s4) = 0.0;
	CSGP4_REF(s5) = 0.0;
	CSGP4_REF(s6) = 0.0;
	CSGP4_REF(s7) = 0.0;
	CSGP4_REF(sz11) = 0.0;
	CSGP4_REF(sz12) = 0.0;
	CSGP4_REF(sz13) = 0.0;
	CSGP4_REF(sz1) = 0.0;
	CSGP4_REF(sz2) = 0.0;
	CSGP4_REF(sz3) = 0.0;
	CSGP4_REF(sz21) = 0.0;
	CSGP4_REF(sz22) = 0.0;
	CSGP4_REF(sz23) = 0.0;
	CSGP4_REF(sz31) = 0.0;
	CSGP4_REF(sz32) = 0.0;
	CSGP4_REF(sz33) = 0.0;
	CSGP4_REF(z13) = 0.0;
	CSGP4_REF(z21) = 0.0;
	CSGP4_REF(z1) = 0.0;
	CSGP4_REF(z2) = 0.0;
	CSGP4_REF(z3) = 0.0;
	CSGP4_REF(z11) = 0.0;
	CSGP4_REF(z12) = 0.0;
	CSGP4_REF(z31) = 0.0;
	CSGP4_REF(z21) = 0.0;
	CSGP4_REF(z22) = 0.0;
	CSGP4_REF(z23) = 0.0;
	CSGP4_REF(z32) = 0.0;
	CSGP4_REF(z33) = 0.0;

	CSGP4_REF(nm) = np;
	CSGP4_REF(em) = ep;
	CSGP4_REF(snodm) = SIN(nodep);
	CSGP4_REF(cnodm) = COS(nodep);
	CSGP4_REF(sinomm) = SIN(argpp);
	CSGP4_REF(cosomm) = COS(argpp);
	CSGP4_REF(sinim) = SIN(inclp);
	CSGP4_REF(cosim) = COS(inclp);
	CSGP4_REF(emsq) = CSGP4_REF(em) * CSGP4_REF(em);
	betasq = 1.0 - CSGP4_REF(emsq);
	CSGP4_REF(rtemsq) = SQRT(betasq);

	/* ----------------- initialize lunar solar terms --------------- */
	CSGP4_REF(peo) = 0.0;
	CSGP4_REF(pinco) = 0.0;
	CSGP4_REF(plo) = 0.0;
	CSGP4_REF(pgho) = 0.0;
	CSGP4_REF(pho) = 0.0;
	CSGP4_REF(day) = epoch + 18261.5 + tc / 1440.0;
	xnodce = FMOD(4.5236020 - 9.2422029e-4 * CSGP4_REF(day), twopi);
	stem = SIN(xnodce);
	ctem = COS(xnodce);
	zcosil = 0.91375164 - 0.03568096 * ctem;
	zsinil = SQRT(1.0 - zcosil * zcosil);
	zsinhl = 0.089683511 * stem / zsinil;
	zcoshl = SQRT(1.0 - zsinhl * zsinhl);
	*gam = 5.8351514 + 0.0019443680 * *day;
	zx = 0.39785416 * stem / zsinil;
	zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
	zx = ATAN2(zx, zy);
	zx = *gam + zx - xnodce;
	zcosgl = COS(zx);
	zsingl = SIN(zx);

	/* ------------------------- do solar terms --------------------- */
	zcosg = zcosgs;
	zsing = zsings;
	zcosi = zcosis;
	zsini = zsinis;
	zcosh = *cnodm;
	zsinh = *snodm;
	cc = c1ss;
	xnoi = 1.0 / *nm;

	for (lsflg = 1; lsflg <= 2; lsflg++)
	{
		a1 = zcosg * zcosh + zsing * zcosi * zsinh;
		a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
		a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
		a8 = zsing * zsini;
		a9 = zsing * zsinh + zcosg * zcosi * zcosh;
		a10 = zcosg * zsini;
		a2 = *cosim * a7 + *sinim * a8;
		a4 = *cosim * a9 + *sinim * a10;
		a5 = -*sinim * a7 + *cosim * a8;
		a6 = -*sinim * a9 + *cosim * a10;

		x1 = a1 * *cosomm + a2 * *sinomm;
		x2 = a3 * *cosomm + a4 * *sinomm;
		x3 = -a1 * *sinomm + a2 * *cosomm;
		x4 = -a3 * *sinomm + a4 * *cosomm;
		x5 = a5 * *sinomm;
		x6 = a6 * *sinomm;
		x7 = a5 * *cosomm;
		x8 = a6 * *cosomm;

		*z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
		*z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
		*z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
		*z1 = 3.0 * (a1 * a1 + a2 * a2) + *z31 * *emsq;
		*z2 = 6.0 * (a1 * a3 + a2 * a4) + *z32 * *emsq;
		*z3 = 3.0 * (a3 * a3 + a4 * a4) + *z33 * *emsq;
		*z11 = -6.0 * a1 * a5 + *emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
		*z12 = -6.0 * (a1 * a6 + a3 * a5) + *emsq *
			   (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
		*z13 = -6.0 * a3 * a6 + *emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
		*z21 = 6.0 * a2 * a5 + *emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
		*z22 = 6.0 * (a4 * a5 + a2 * a6) + *emsq *
			   (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
		*z23 = 6.0 * a4 * a6 + *emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
		*z1 = *z1 + *z1 + betasq * *z31;
		*z2 = *z2 + *z2 + betasq * *z32;
		*z3 = *z3 + *z3 + betasq * *z33;
		*s3 = cc * xnoi;
		*s2 = -0.5 * *s3 / *rtemsq;
		*s4 = *s3 * *rtemsq;
		*s1 = -15.0 * *em * *s4;
		*s5 = x1 * x3 + x2 * x4;
		*s6 = x2 * x3 + x1 * x4;
		*s7 = x2 * x4 - x1 * x3;

		/* ----------------------- do lunar terms ------------------- */
		if (lsflg == 1)
		{
			*ss1 = *s1;
			*ss2 = *s2;
			*ss3 = *s3;
			*ss4 = *s4;
			*ss5 = *s5;
			*ss6 = *s6;
			*ss7 = *s7;
			*sz1 = *z1;
			*sz2 = *z2;
			*sz3 = *z3;
			*sz11 = *z11;
			*sz12 = *z12;
			*sz13 = *z13;
			*sz21 = *z21;
			*sz22 = *z22;
			*sz23 = *z23;
			*sz31 = *z31;
			*sz32 = *z32;
			*sz33 = *z33;
			zcosg = zcosgl;
			zsing = zsingl;
			zcosi = zcosil;
			zsini = zsinil;
			zcosh = zcoshl * *cnodm + zsinhl * *snodm;
			zsinh = *snodm * zcoshl - *cnodm * zsinhl;
			cc = c1l;
		}
	}

	*zmol = fmod( (4.7199672 + 0.22997150 * *day - *gam), twopi );
	*zmos = fmod( (6.2565837 + 0.017201977 * *day), twopi );

	/* ------------------------ do solar terms ---------------------- */
	*se2 = 2.0 * *ss1 * *ss6;
	*se3 = 2.0 * *ss1 * *ss7;
	*si2 = 2.0 * *ss2 * *sz12;
	*si3 = 2.0 * *ss2 * (*sz13 - *sz11);
	*sl2 = -2.0 * *ss3 * *sz2;
	*sl3 = -2.0 * *ss3 * (*sz3 - *sz1);
	*sl4 = -2.0 * *ss3 * (-21.0 - 9.0 * *emsq) * zes;
	*sgh2 = 2.0 * *ss4 * *sz32;
	*sgh3 = 2.0 * *ss4 * (*sz33 - *sz31);
	*sgh4 = -18.0 * *ss4 * zes;
	*sh2 = -2.0 * *ss2 * *sz22;
	*sh3 = -2.0 * *ss2 * (*sz23 - *sz21);

	/* ------------------------ do lunar terms ---------------------- */
	*ee2 = 2.0 * *s1 * *s6;
	*e3 = 2.0 * *s1 * *s7;
	*xi2 = 2.0 * *s2 * *z12;
	*xi3 = 2.0 * *s2 * (*z13 - *z11);
	*xl2 = -2.0 * *s3 * *z2;
	*xl3 = -2.0 * *s3 * (*z3 - *z1);
	*xl4 = -2.0 * *s3 * (-21.0 - 9.0 * *emsq) * zel;
	*xgh2 = 2.0 * *s4 * *z32;
	*xgh3 = 2.0 * *s4 * (*z33 - *z31);
	*xgh4 = -18.0 * *s4 * zel;
	*xh2 = -2.0 * *s2 * *z22;
	*xh3 = -2.0 * *s2 * (*z23 - *z21);

	//#include "debug2.cpp"
}  // end dscom


CSGP4_DECORATOR void dsinit_simple
	 (
	// sgp4fix just send in xke as a constant and eliminate getgravconst call
	// gravconsttype whichconst, 
	   SGPF xke,
	   SGPF cosim, SGPF emsq, SGPF argpo, SGPF s1, SGPF s2,
	   SGPF s3, SGPF s4, SGPF s5, SGPF sinim, SGPF ss1,
	   SGPF ss2, SGPF ss3, SGPF ss4, SGPF ss5, SGPF sz1,
	   SGPF sz3, SGPF sz11, SGPF sz13, SGPF sz21, SGPF sz23,
	   SGPF sz31, SGPF sz33, SGPF t, SGPF tc, SGPF gsto,
	   SGPF mo, SGPF mdot, SGPF no, SGPF nodeo, SGPF nodedot,
	   SGPF xpidot, SGPF z1, SGPF z3, SGPF z11, SGPF z13,
	   SGPF z21, SGPF z23, SGPF z31, SGPF z33, SGPF ecco,
	   SGPF eccsq, SGPF CSGP4_OUT em, SGPF CSGP4_OUT argpm, SGPF CSGP4_OUT inclm, SGPF CSGP4_OUT mm,
	   SGPF CSGP4_OUT nm, SGPF CSGP4_OUT nodem,
	   int * irez,
	   SGPF CSGP4_OUT atime, SGPF CSGP4_OUT d2201, SGPF CSGP4_OUT d2211, SGPF CSGP4_OUT d3210, SGPF CSGP4_OUT d3222,
	   SGPF CSGP4_OUT d4410, SGPF CSGP4_OUT d4422, SGPF CSGP4_OUT d5220, SGPF CSGP4_OUT d5232, SGPF CSGP4_OUT d5421,
	   SGPF CSGP4_OUT d5433, SGPF CSGP4_OUT dedt, SGPF CSGP4_OUT didt, SGPF CSGP4_OUT dmdt, SGPF CSGP4_OUT dndt,
	   SGPF CSGP4_OUT dnodt, SGPF CSGP4_OUT domdt, SGPF CSGP4_OUT del1, SGPF CSGP4_OUT del2, SGPF CSGP4_OUT del3,
	   SGPF CSGP4_OUT xfact, SGPF CSGP4_OUT xlamo, SGPF CSGP4_OUT xli, SGPF CSGP4_OUT xni
	 )
{
	/* --------------------- local variables ------------------------ */
	const SGPF twopi = 2.0 * SGPPI;

	SGPF ainv2, aonv = 0.0, cosisq, eoc, f220, f221, f311,
		 f321, f322, f330, f441, f442, f522, f523,
		 f542, f543, g200, g201, g211, g300, g310,
		 g322, g410, g422, g520, g521, g532, g533,
		 ses, sgs, sghl, sghs, shs, shll, sis,
		 sini2, sls, temp, temp1, theta, xno2, q22,
		 q31, q33, root22, root44, root54, rptim, root32,
		 root52, x2o3, znl, emo, zns, emsqo;

	q22 = 1.7891679e-6;
	q31 = 2.1460748e-6;
	q33 = 2.2123015e-7;
	root22 = 1.7891679e-6;
	root44 = 7.3636953e-9;
	root54 = 2.1765803e-9;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	root32 = 3.7393792e-7;
	root52 = 1.1428639e-7;
	x2o3 = 2.0 / 3.0;
	znl = 1.5835218e-4;
	zns = 1.19459e-5;

	// sgp4fix identify constants and allow alternate values
	// just xke is used here so pass it in rather than have multiple calls
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

	/* -------------------- deep space initialization ------------ */
	*irez = 0;
	if ((*nm < 0.0052359877) && (*nm > 0.0034906585))
		*irez = 1;
	if ((*nm >= 8.26e-3) && (*nm <= 9.24e-3) && (*em >= 0.5))
		*irez = 2;

	/* ------------------------ do solar terms ------------------- */
	ses = ss1 * zns * ss5;
	sis = ss2 * zns * (sz11 + sz13);
	sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
	sghs = ss4 * zns * (sz31 + sz33 - 6.0);
	shs = -zns * ss2 * (sz21 + sz23);
	// sgp4fix for 180 deg incl
	if ((*inclm < 5.2359877e-2) || (*inclm > SGPPI - 5.2359877e-2))
		shs = 0.0;
	if (sinim != 0.0)
		shs = shs / sinim;
	sgs = sghs - cosim * shs;

	/* ------------------------- do lunar terms ------------------ */
	*dedt = ses + s1 * znl * s5;
	*didt = sis + s2 * znl * (z11 + z13);
	*dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
	sghl = s4 * znl * (z31 + z33 - 6.0);
	shll = -znl * s2 * (z21 + z23);
	// sgp4fix for 180 deg incl
	if ((*inclm < 5.2359877e-2) || (*inclm > SGPPI - 5.2359877e-2))
		shll = 0.0;
	*domdt = sgs + sghl;
	*dnodt = shs;
	if (sinim != 0.0)
	{
		*domdt = *domdt - cosim / sinim * shll;
		*dnodt = *dnodt + shll / sinim;
	}

	/* ----------- calculate deep space resonance effects -------- */
	*dndt = 0.0;
	theta = fmod((gsto + tc * rptim), twopi );
	*em = *em + *dedt * t;
	*inclm = *inclm + *didt * t;
	*argpm = *argpm + *domdt * t;
	*nodem = *nodem + *dnodt * t;
	*mm = *mm + *dmdt * t;
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//if (inclm < 0.0)
	//  {
	//	inclm  = -inclm;
	//	argpm  = argpm - Math.PI;
	//	nodem = nodem + Math.PI;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if (*irez != 0)
	{
		aonv = POW(*nm / xke, x2o3);

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if (*irez == 2)
		{
			cosisq = cosim * cosim;
			emo = *em;
			*em = ecco;
			emsqo = emsq;
			emsq = eccsq;
			eoc = *em * emsq;
			g201 = -0.306 - (*em - 0.64) * 0.440;

			if (*em <= 0.65)
			{
				g211 = 3.616 - 13.2470 * *em + 16.2900 * emsq;
				g310 = -19.302 + 117.3900 * *em - 228.4190 * emsq + 156.5910 * eoc;
				g322 = -18.9068 + 109.7927 * *em - 214.6334 * emsq + 146.5816 * eoc;
				g410 = -41.122 + 242.6940 * *em - 471.0940 * emsq + 313.9530 * eoc;
				g422 = -146.407 + 841.8800 * *em - 1629.014 * emsq + 1083.4350 * eoc;
				g520 = -532.114 + 3017.977 * *em - 5740.032 * emsq + 3708.2760 * eoc;
			}
			else
			{
				g211 = -72.099 + 331.819 * *em - 508.738 * emsq + 266.724 * eoc;
				g310 = -346.844 + 1582.851 * *em - 2415.925 * emsq + 1246.113 * eoc;
				g322 = -342.585 + 1554.908 * *em - 2366.899 * emsq + 1215.972 * eoc;
				g410 = -1052.797 + 4758.686 * *em - 7193.992 * emsq + 3651.957 * eoc;
				g422 = -3581.690 + 16178.110 * *em - 24462.770 * emsq + 12422.520 * eoc;
				if (*em > 0.715)
					g520 = -5149.66 + 29936.92 * *em - 54087.36 * emsq + 31324.56 * eoc;
				else
					g520 = 1464.74 - 4664.75 * *em + 3763.64 * emsq;
			}
			if (*em < 0.7)
			{
				g533 = -919.22770 + 4988.6100 * *em - 9064.7700 * emsq + 5542.21 * eoc;
				g521 = -822.71072 + 4568.6173 * *em - 8491.4146 * emsq + 5337.524 * eoc;
				g532 = -853.66600 + 4690.2500 * *em - 8624.7700 * emsq + 5341.4 * eoc;
			}
			else
			{
				g533 = -37995.780 + 161616.52 * *em - 229838.20 * emsq + 109377.94 * eoc;
				g521 = -51752.104 + 218913.95 * *em - 309468.16 * emsq + 146349.42 * eoc;
				g532 = -40023.880 + 170470.89 * *em - 242699.48 * emsq + 115605.82 * eoc;
			}

			sini2 = sinim * sinim;
			f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
			f221 = 1.5 * sini2;
			f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
			f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
			f441 = 35.0 * sini2 * f220;
			f442 = 39.3750 * sini2 * sini2;
			f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
					0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
			f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
				   10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
			f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq *
				   (-12.0 + 8.0 * cosim + 10.0 * cosisq));
			f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq *
				   (12.0 + 8.0 * cosim - 10.0 * cosisq));
			xno2 = *nm * *nm;
			ainv2 = aonv * aonv;
			temp1 = 3.0 * xno2 * ainv2;
			temp = temp1 * root22;
			*d2201 = temp * f220 * g201;
			*d2211 = temp * f221 * g211;
			temp1 = temp1 * aonv;
			temp = temp1 * root32;
			*d3210 = temp * f321 * g310;
			*d3222 = temp * f322 * g322;
			temp1 = temp1 * aonv;
			temp = 2.0 * temp1 * root44;
			*d4410 = temp * f441 * g410;
			*d4422 = temp * f442 * g422;
			temp1 = temp1 * aonv;
			temp = temp1 * root52;
			*d5220 = temp * f522 * g520;
			*d5232 = temp * f523 * g532;
			temp = 2.0 * temp1 * root54;
			*d5421 = temp * f542 * g521;
			*d5433 = temp * f543 * g533;
			*xlamo = fmod( (mo + nodeo + nodeo - theta - theta), twopi );
			*xfact = mdot + *dmdt + 2.0 * (nodedot + *dnodt - rptim) - no;
			*em = emo;
			emsq = emsqo;
		}

		/* ---------------- synchronous resonance terms -------------- */
		if (*irez == 1)
		{
			g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
			g310 = 1.0 + 2.0 * emsq;
			g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
			f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
			f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
			f330 = 1.0 + cosim;
			f330 = 1.875 * f330 * f330 * f330;
			*del1 = 3.0 * *nm * *nm * aonv * aonv;
			*del2 = 2.0 * *del1 * f220 * g200 * q22;
			*del3 = 3.0 * *del1 * f330 * g300 * q33 * aonv;
			*del1 = *del1 * f311 * g310 * q31 * aonv;
			*xlamo = fmod(mo + nodeo + argpo - theta, twopi );
			*xfact = mdot + xpidot - rptim + *dmdt + *domdt + *dnodt - no;
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		*xli = *xlamo;
		*xni = no;
		*atime = 0.0;
		*nm = no + *dndt;
	}

	//#include "debug3.cpp"
}  // end dsinit


CSGP4_DECORATOR void getgravconst_simple
	 (
	  SGPF CSGP4_OUT tumin,
	  SGPF CSGP4_OUT mu,
	  SGPF CSGP4_OUT radiusearthkm,
	  SGPF CSGP4_OUT xke,
	  SGPF CSGP4_OUT j2,
	  SGPF CSGP4_OUT j3,
	  SGPF CSGP4_OUT j4,
	  SGPF CSGP4_OUT j3oj2
	 )
{
	//Hard-coded WGS72
			*mu = 398600.8;			// in km3 / s2
			*radiusearthkm = 6378.135;	 // km
			*xke = 60.0 / SQRT(*radiusearthkm * *radiusearthkm * *radiusearthkm / *mu);
			*tumin = 1.0 / *xke;
			*j2 = 0.001082616;
			*j3 = -0.00000253881;
			*j4 = -0.00000165597;
			*j3oj2 = *j3 / *j2;

}   // end getgravconst


CSGP4_DECORATOR void sgp4init_simple
	 (
	   /*enum gravconsttype whichconst,*/ SGPF epoch,
	   SGPF xbstar, SGPF xndot, SGPF xnddot, SGPF xecco, SGPF xargpo,
	   SGPF xinclo, SGPF xmo, SGPF xno_kozai,
	   SGPF xnodeo, SGPF initial_time, SGPF * initial_r, SGPF * initial_v,
	   SGPF * a_alta_altp
	 )
{
	struct elsetrec_simple g = { 0 };
	struct elsetrec_simple * satrec = &g;

	/* --------------------- local variables ------------------------ */
	SGPF ao, ainv, con42, cosio, sinio, cosio2, eccsq,
		 omeosq, posq, rp, rteosq,
		 cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq,
		 cc2, cc3, coef, coef1, cosio4, day, dndt,
		 em, emsq, eeta, etasq, gam, argpm, nodem,
		 inclm, mm, nm, perige, pinvsq, psisq, qzms24,
		 rtemsq, s1, s2, s3, s4, s5, s6,
		 s7, sfour, ss1, ss2, ss3, ss4, ss5,
		 ss6, ss7, sz1, sz2, sz3, sz11, sz12,
		 sz13, sz21, sz22, sz23, sz31, sz32, sz33,
		 tc, temp, temp1, temp2, temp3, tsi, xpidot,
		 xhdot1, z1, z2, z3, z11, z12, z13,
		 z21, z22, z23, z31, z32, z33,
		 qzms2t, ss, x2o3,
		 delmotemp, qzms2ttemp, qzms24temp;

	/* ------------------------ initialization --------------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(Math.PI-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	const SGPF temp4 = 1.5e-12;

	/* ----------- set all near earth variables to zero ------------ */
	/*satrec->isimp = 0; */ satrec->aycof = 0.0;
	satrec->con41 = 0.0; satrec->cc1 = 0.0; satrec->cc4 = 0.0;
	satrec->cc5 = 0.0; satrec->d2 = 0.0; satrec->d3 = 0.0;
	satrec->d4 = 0.0; satrec->delmo = 0.0; satrec->eta = 0.0;
	satrec->argpdot = 0.0; satrec->omgcof = 0.0; satrec->sinmao = 0.0;
	satrec->t = 0.0; satrec->t2cof = 0.0; satrec->t3cof = 0.0;
	satrec->t4cof = 0.0; satrec->t5cof = 0.0; satrec->x1mth2 = 0.0;
	satrec->x7thm1 = 0.0; satrec->mdot = 0.0; satrec->nodedot = 0.0;
	satrec->xlcof = 0.0; satrec->xmcof = 0.0; satrec->nodecf = 0.0;

	/* ----------- set all deep space variables to zero ------------ */
	satrec->irez = 0; satrec->d2201 = 0.0; satrec->d2211 = 0.0;
	satrec->d3210 = 0.0; satrec->d3222 = 0.0; satrec->d4410 = 0.0;
	satrec->d4422 = 0.0; satrec->d5220 = 0.0; satrec->d5232 = 0.0;
	satrec->d5421 = 0.0; satrec->d5433 = 0.0; satrec->dedt = 0.0;
	satrec->del1 = 0.0; satrec->del2 = 0.0; satrec->del3 = 0.0;
	satrec->didt = 0.0; satrec->dmdt = 0.0; satrec->dnodt = 0.0;
	satrec->domdt = 0.0; satrec->e3 = 0.0; satrec->ee2 = 0.0;
	satrec->peo = 0.0; satrec->pgho = 0.0; satrec->pho = 0.0;
	satrec->pinco = 0.0; satrec->plo = 0.0; satrec->se2 = 0.0;
	satrec->se3 = 0.0; satrec->sgh2 = 0.0; satrec->sgh3 = 0.0;
	satrec->sgh4 = 0.0; satrec->sh2 = 0.0; satrec->sh3 = 0.0;
	satrec->si2 = 0.0; satrec->si3 = 0.0; satrec->sl2 = 0.0;
	satrec->sl3 = 0.0; satrec->sl4 = 0.0; satrec->gsto = 0.0;
	satrec->xfact = 0.0; satrec->xgh2 = 0.0; satrec->xgh3 = 0.0;
	satrec->xgh4 = 0.0; satrec->xh2 = 0.0; satrec->xh3 = 0.0;
	satrec->xi2 = 0.0; satrec->xi3 = 0.0; satrec->xl2 = 0.0;
	satrec->xl3 = 0.0; satrec->xl4 = 0.0; satrec->xlamo = 0.0;
	satrec->zmol = 0.0; satrec->zmos = 0.0; satrec->atime = 0.0;
	satrec->xli = 0.0; satrec->xni = 0.0;

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow aernate values
	// this is now the only call for the constants
	getgravconst_simple(/*whichconst,*/ &satrec->tumin, &satrec->mu, &satrec->radiusearthkm, &satrec->xke,
				  &satrec->j2, &satrec->j3, &satrec->j4, &satrec->j3oj2);
	//-------------------------------------------------------------------------
	satrec->error = 0;
	//satrec->operationmode = opsmode;
	//strcpy( satrec->satnum, satn );

	// sgp4fix - note the following variables are also passed directly via satrec->
	// it is possible to streamline the sgp4init call by deleting the 'x'
	// variables, but the user would need to set the satrec* values first. we
	// include the additional assignments in case twoline2rv is not used.
	satrec->bstar = xbstar;
	// sgp4fix allow additional parameters in the struct
	satrec->ndot = xndot;
	satrec->nddot = xnddot;
	satrec->ecco = xecco;
	satrec->argpo = xargpo;
	satrec->inclo = xinclo;
	satrec->mo = xmo;
	// sgp4fix rename variables to clarify which mean motion is intended
	satrec->no_kozai = xno_kozai;
	satrec->nodeo = xnodeo;

	// single averaged mean elements
	satrec->am = satrec->em = satrec->im = satrec->Om = satrec->mm = satrec->nm = 0.0;

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow alternate values no longer needed
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	ss = 78.0 / satrec->radiusearthkm + 1.0;
	// sgp4fix use multiply for speed instead of Math.POW
	qzms2ttemp = (120.0 - 78.0) / satrec->radiusearthkm;
	qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
	x2o3 = 2.0 / 3.0;

	//satrec->init = 'y';
	satrec->t = 0.0;

	// sgp4fix remove satn as it is not needed in initl
	initl_simple
		(satrec->xke, satrec->j2, satrec->ecco, epoch, satrec->inclo, satrec->no_kozai, &satrec->method,
		  &ainv, &ao, &satrec->con41, &con42, &cosio, &cosio2, &eccsq, &omeosq,
		  &posq, &rp, &rteosq, &sinio, &satrec->gsto, &satrec->no_unkozai);
	satrec->a = POW(satrec->no_unkozai * satrec->tumin, (-2.0 / 3.0));
	satrec->alta = satrec->a * (1.0 + satrec->ecco) - 1.0;
	satrec->altp = satrec->a * (1.0 - satrec->ecco) - 1.0;
	satrec->error = 0;

	if( a_alta_altp )
	{
		a_alta_altp[0] = satrec->a;
		a_alta_altp[1] = satrec->alta;
		a_alta_altp[2] = satrec->altp;
	}

	// sgp4fix remove this check as it is unnecessary
	// the mrt check in sgp4 handles decaying satellite cases even if the starting
	// condition is below the surface of te earth
	//	 if (rp < 1.0)
	//	   {
	//		 printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
	//		 satrec->error = 5;
	//	   }

	if ((omeosq >= 0.0) || (satrec->no_unkozai >= 0.0))
	{
		satrec->isimp = 0;
		if (rp < (220.0 / satrec->radiusearthkm + 1.0))
			satrec->isimp = 1;
		sfour = ss;
		qzms24 = qzms2t;
		perige = (rp - 1.0) * satrec->radiusearthkm;

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if (perige < 156.0)
		{
			sfour = perige - 78.0;
			if (perige < 98.0)
				sfour = 20.0;
			// sgp4fix use multiply for speed instead of Math.POW
			qzms24temp = (120.0 - sfour) / satrec->radiusearthkm;
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
			sfour = sfour / satrec->radiusearthkm + 1.0;
		}
		pinvsq = 1.0 / posq;

		tsi = 1.0 / (ao - sfour);
		satrec->eta = ao * satrec->ecco * tsi;
		etasq = satrec->eta * satrec->eta;
		eeta = satrec->ecco * satrec->eta;
		psisq = FABS(1.0 - etasq);
		coef = qzms24 * POW(tsi, 4.0);
		coef1 = coef / POW(psisq, 3.5);
		cc2 = coef1 * satrec->no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
					   (4.0 + etasq)) + 0.375 * satrec->j2 * tsi / psisq * satrec->con41 *
					   (8.0 + 3.0 * etasq * (8.0 + etasq)));
		satrec->cc1 = satrec->bstar * cc2;
		cc3 = 0.0;
		if (satrec->ecco > 1.0e-4)
			cc3 = -2.0 * coef * tsi * satrec->j3oj2 * satrec->no_unkozai * sinio / satrec->ecco;
		satrec->x1mth2 = 1.0 - cosio2;
		satrec->cc4 = 2.0 * satrec->no_unkozai * coef1 * ao * omeosq *
						  (satrec->eta * (2.0 + 0.5 * etasq) + satrec->ecco *
						  (0.5 + 2.0 * etasq) - satrec->j2 * tsi / (ao * psisq) *
						  (-3.0 * satrec->con41 * (1.0 - 2.0 * eeta + etasq *
						  (1.5 - 0.5 * eeta)) + 0.75 * satrec->x1mth2 *
						  (2.0 * etasq - eeta * (1.0 + etasq)) * COS(2.0 * satrec->argpo)));
		satrec->cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
					   (etasq + eeta) + eeta * etasq);
		cosio4 = cosio2 * cosio2;
		temp1 = 1.5 * satrec->j2 * pinvsq * satrec->no_unkozai;
		temp2 = 0.5 * temp1 * satrec->j2 * pinvsq;
		temp3 = -0.46875 * satrec->j4 * pinvsq * pinvsq * satrec->no_unkozai;
		satrec->mdot = satrec->no_unkozai + 0.5 * temp1 * rteosq * satrec->con41 + 0.0625 *
						   temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
		satrec->argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
							(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
							temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
		xhdot1 = -temp1 * cosio;
		satrec->nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
							 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
		xpidot = satrec->argpdot + satrec->nodedot;
		satrec->omgcof = satrec->bstar * cc3 * COS(satrec->argpo);
		satrec->xmcof = 0.0;
		if (satrec->ecco > 1.0e-4)
			satrec->xmcof = -x2o3 * coef * satrec->bstar / eeta;
		satrec->nodecf = 3.5 * omeosq * xhdot1 * satrec->cc1;
		satrec->t2cof = 1.5 * satrec->cc1;
		// sgp4fix for divide by zero with xinco = 180 deg
		if (FABS(cosio + 1.0) > 1.5e-12)
			satrec->xlcof = -0.25 * satrec->j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
		else
			satrec->xlcof = -0.25 * satrec->j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
		satrec->aycof = -0.5 * satrec->j3oj2 * sinio;
		// sgp4fix use multiply for speed instead of Math.POW
		delmotemp = 1.0 + satrec->eta * COS(satrec->mo);
		satrec->delmo = delmotemp * delmotemp * delmotemp;
		satrec->sinmao = SIN(satrec->mo);
		satrec->x7thm1 = 7.0 * cosio2 - 1.0;

		/* --------------- deep space initialization ------------- */
		if ((2 * SGPPI / satrec->no_unkozai) >= 225.0)
		{
			satrec->method = true;
			satrec->isimp = 1;
			tc = 0.0;
			inclm = satrec->inclo;

			dscom_simple
				(
				  epoch, satrec->ecco, satrec->argpo, tc, satrec->inclo, satrec->nodeo,
				  satrec->no_unkozai, &snodm, &cnodm, &sinim, &cosim, &sinomm, &cosomm,
				  &day, &satrec->e3, &satrec->ee2, &em, &emsq, &gam,
				  &satrec->peo, &satrec->pgho, &satrec->pho, &satrec->pinco,
				  &satrec->plo, &rtemsq, &satrec->se2, &satrec->se3,
				  &satrec->sgh2, &satrec->sgh3, &satrec->sgh4,
				  &satrec->sh2, & satrec->sh3, &satrec->si2, &satrec->si3,
				  &satrec->sl2, &satrec->sl3, &satrec->sl4, &s1, &s2, &s3, &s4, &s5,
				  &s6, &s7, &ss1, &ss2, &ss3, &ss4, &ss5, &ss6, &ss7, &sz1, &sz2, &sz3,
				  &sz11, &sz12, &sz13, &sz21, &sz22, &sz23, &sz31, &sz32, &sz33,
				  &satrec->xgh2, &satrec->xgh3, &satrec->xgh4, &satrec->xh2,
				  &satrec->xh3, &satrec->xi2, &satrec->xi3, &satrec->xl2,
				  &satrec->xl3, &satrec->xl4, &nm, &z1, &z2, &z3, &z11,
				  &z12, &z13, &z21, &z22, &z23, &z31, &z32, &z33,
				  &satrec->zmol, &satrec->zmos
				);
			dpper_simple
				(
				  satrec->e3, satrec->ee2, satrec->peo, satrec->pgho,
				  satrec->pho, satrec->pinco, satrec->plo, satrec->se2,
				  satrec->se3, satrec->sgh2, satrec->sgh3, satrec->sgh4,
				  satrec->sh2, satrec->sh3, satrec->si2, satrec->si3,
				  satrec->sl2, satrec->sl3, satrec->sl4, satrec->t,
				  satrec->xgh2, satrec->xgh3, satrec->xgh4, satrec->xh2,
				  satrec->xh3, satrec->xi2, satrec->xi3, satrec->xl2,
				  satrec->xl3, satrec->xl4, satrec->zmol, satrec->zmos, inclm, true,
				  &satrec->ecco, &satrec->inclo, &satrec->nodeo, &satrec->argpo, &satrec->mo
				);

			argpm = 0.0;
			nodem = 0.0;
			mm = 0.0;

			dsinit_simple
				(
				  satrec->xke,
				  cosim, emsq, satrec->argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
				  ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec->t, tc,
				  satrec->gsto, satrec->mo, satrec->mdot, satrec->no_unkozai, satrec->nodeo,
				  satrec->nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
				  satrec->ecco, eccsq, &em, &argpm, &inclm, &mm, &nm, &nodem,
				  &satrec->irez, &satrec->atime,
				  &satrec->d2201, &satrec->d2211, &satrec->d3210, &satrec->d3222,
				  &satrec->d4410, &satrec->d4422, &satrec->d5220, &satrec->d5232,
				  &satrec->d5421, &satrec->d5433, &satrec->dedt, &satrec->didt,
				  &satrec->dmdt, &dndt, &satrec->dnodt, &satrec->domdt,
				  &satrec->del1, &satrec->del2, &satrec->del3, &satrec->xfact,
				  &satrec->xlamo, &satrec->xli, &satrec->xni
				);
		}

		/* ----------- set variables if not deep space ----------- */
		if (satrec->isimp != 1)
		{
			cc1sq = satrec->cc1 * satrec->cc1;
			satrec->d2 = 4.0 * ao * tsi * cc1sq;
			temp = satrec->d2 * tsi * satrec->cc1 / 3.0;
			satrec->d3 = (17.0 * ao + sfour) * temp;
			satrec->d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
							 satrec->cc1;
			satrec->t3cof = satrec->d2 + 2.0 * cc1sq;
			satrec->t4cof = 0.25 * (3.0 * satrec->d3 + satrec->cc1 *
							 (12.0 * satrec->d2 + 10.0 * cc1sq));
			satrec->t5cof = 0.2 * (3.0 * satrec->d4 +
							 12.0 * satrec->cc1 * satrec->d3 +
							 6.0 * satrec->d2 * satrec->d2 +
							 15.0 * cc1sq * (2.0 * satrec->d2 + cc1sq));
		}
	} // if omeosq = 0 ...


	/* finally propagate to zero epoch to initialize all others. */
	// sgp4fix take out check to let satellites process until they are actually below earth surface
	//	   if(satrec->error == 0)
	SGPF r[3];
	SGPF v[3];
	if( !initial_r ) initial_r = r;
	if( !initial_v ) initial_v = v;
	sgp4_simple(satrec, initial_time, initial_r, initial_v);

	//satrec->init = 'n';

	//#include "debug6.cpp"
	//sgp4fix return boolean. satrec->error contains any error codes
	//return true;
}  // end sgp4init

#endif

