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
#include <fenv.h>

#ifndef CSGP4_OUT
#define CSGP4_OUT * restrict
#endif

#ifndef CSGP4_DEREF
#define CSGP4_DEREF(x) (*x)
#endif

#ifndef CSGP4_REF
#define CSGP4_REF(x) (&x)
#endif

#define SGPPI 3.1415926535897932384626433

// Default to double calculations.
#ifndef CSGP4_USE_FLOAT
#define CSGP4_USE_FLOAT 0
#endif

#ifndef CSGP4_INIT
#define CSGP4_INIT 1
#endif


#ifndef CSGP4_DECORATOR
#define CSGP4_DECORATOR static inline
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
#define FMOD  fmod
#endif


#define BOOL int


struct elsetrec_simple
{
	int error;
	bool method; // true for 'd' false for 'n' (Deep / Near Eart)

	/* Near Earth */
	bool isimp;

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

#if 0
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
		CSGP4_DEREF(inclp) = CSGP4_DEREF(inclp) + pinc;
		CSGP4_DEREF(ep) = CSGP4_DEREF(ep) + pe;
		sinip = SIN(CSGP4_DEREF(inclp));
		cosip = COS(CSGP4_DEREF(inclp));

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if (CSGP4_DEREF(inclp) >= 0.2)
		{
			ph = ph / sinip;
			pgh = pgh - cosip * ph;
			CSGP4_DEREF(argpp) = CSGP4_DEREF(argpp) + pgh;
			CSGP4_DEREF(nodep) = CSGP4_DEREF(nodep) + ph;
			CSGP4_DEREF(mp) = CSGP4_DEREF(mp) + pl;
		}
		else
		{
			/* ---- apply periodics with lyddane modification ---- */
			sinop = SIN(CSGP4_DEREF(nodep));
			cosop = COS(CSGP4_DEREF(nodep));
			alfdp = sinip * sinop;
			betdp = sinip * cosop;
			dalf = ph * cosop + pinc * cosip * sinop;
			dbet = -ph * sinop + pinc * cosip * cosop;
			alfdp = alfdp + dalf;
			betdp = betdp + dbet;
			CSGP4_DEREF(nodep) = FMOD( CSGP4_DEREF(nodep), twopi );
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if ((CSGP4_DEREF(nodep) < 0.0)) // && (opsmode == 'a'))
				CSGP4_DEREF(nodep) = CSGP4_DEREF(nodep) + twopi;
			xls = CSGP4_DEREF(mp) + CSGP4_DEREF(argpp) + cosip * CSGP4_DEREF(nodep);
			dls = pl + pgh - pinc * CSGP4_DEREF(nodep) * sinip;
			xls = xls + dls;
			xnoh = CSGP4_DEREF(nodep);
			CSGP4_DEREF(nodep) = ATAN2(alfdp, betdp);
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if ((CSGP4_DEREF(nodep) < 0.0)) // && (opsmode == 'a'))
				CSGP4_DEREF(nodep) = CSGP4_DEREF(nodep) + twopi;
			if (FABS(xnoh - CSGP4_DEREF(nodep)) > SGPPI)
			{
				if (CSGP4_DEREF(nodep) < xnoh)
					CSGP4_DEREF(nodep) = CSGP4_DEREF(nodep) + twopi;
				else
					CSGP4_DEREF(nodep) = CSGP4_DEREF(nodep) - twopi;
			}
			CSGP4_DEREF(mp) = CSGP4_DEREF(mp) + pl;
			CSGP4_DEREF(argpp) = xls - CSGP4_DEREF(mp) - cosip * CSGP4_DEREF(nodep);
		}
	}   // if init == 'n'
	//#include "debug1.cpp"
}  // end dpper
#endif

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
	CSGP4_DEREF(eccsq) = ecco * ecco;
	CSGP4_DEREF(omeosq) = 1.0 - CSGP4_DEREF(eccsq);
	CSGP4_DEREF(rteosq) = SQRT(CSGP4_DEREF(omeosq));
	CSGP4_DEREF(cosio) = COS(inclo);
	CSGP4_DEREF(cosio2) = CSGP4_DEREF(cosio) * CSGP4_DEREF(cosio);

	/* ------------------ un-kozai the mean motion ----------------- */
	ak = POW(xke / no_kozai, x2o3);
	d1 = 0.75 * j2 * (3.0 * CSGP4_DEREF(cosio2) - 1.0) / (CSGP4_DEREF(rteosq) * CSGP4_DEREF(omeosq));
	del = d1 / (ak * ak);
	adel = ak * (1.0 - del * del - del *
			(1.0 / 3.0 + 134.0 * del * del / 81.0));
	del = d1 / (adel * adel);
	CSGP4_DEREF(no_unkozai) = no_kozai / (1.0 + del);

	CSGP4_DEREF(ao) = POW(xke / (CSGP4_DEREF(no_unkozai)), x2o3);
	CSGP4_DEREF(sinio) = SIN(inclo);
	po = CSGP4_DEREF(ao) * CSGP4_DEREF(omeosq);
	CSGP4_DEREF(con42) = 1.0 - 5.0 * CSGP4_DEREF(cosio2);
	CSGP4_DEREF(con41) = -CSGP4_DEREF(con42) - CSGP4_DEREF(cosio2) - CSGP4_DEREF(cosio2);
	CSGP4_DEREF(ainv) = 1.0 / CSGP4_DEREF(ao);
	CSGP4_DEREF(posq) = po * po;
	CSGP4_DEREF(rp) = CSGP4_DEREF(ao) * (1.0 - ecco);
	CSGP4_DEREF(method) = false;
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
		CSGP4_DEREF(gsto) = FMOD(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
		if (CSGP4_DEREF(gsto) < 0.0)
			CSGP4_DEREF(gsto) = CSGP4_DEREF(gsto) + twopi;
	}
	//else
	//	*gsto = gstime(epoch + 2433281.5);

	//#include "debug5.cpp"
}  // end initl


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
	   int CSGP4_OUT irez,
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
	CSGP4_DEREF(irez) = 0;
	if ((CSGP4_DEREF(nm) < 0.0052359877) && (CSGP4_DEREF(nm) > 0.0034906585))
		CSGP4_DEREF(irez) = 1;
	if ((CSGP4_DEREF(nm) >= 8.26e-3) && (CSGP4_DEREF(nm) <= 9.24e-3) && (CSGP4_DEREF(em) >= 0.5))
		CSGP4_DEREF(irez) = 2;

	/* ------------------------ do solar terms ------------------- */
	ses = ss1 * zns * ss5;
	sis = ss2 * zns * (sz11 + sz13);
	sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
	sghs = ss4 * zns * (sz31 + sz33 - 6.0);
	shs = -zns * ss2 * (sz21 + sz23);
	// sgp4fix for 180 deg incl
	if ((CSGP4_DEREF(inclm) < 5.2359877e-2) || (CSGP4_DEREF(inclm) > SGPPI - 5.2359877e-2))
		shs = 0.0;
	if (sinim != 0.0)
		shs = shs / sinim;
	sgs = sghs - cosim * shs;

	/* ------------------------- do lunar terms ------------------ */
	CSGP4_DEREF(dedt) = ses + s1 * znl * s5;
	CSGP4_DEREF(didt) = sis + s2 * znl * (z11 + z13);
	CSGP4_DEREF(dmdt) = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
	sghl = s4 * znl * (z31 + z33 - 6.0);
	shll = -znl * s2 * (z21 + z23);
	// sgp4fix for 180 deg incl
	if ((CSGP4_DEREF(inclm) < 5.2359877e-2) || (CSGP4_DEREF(inclm) > SGPPI - 5.2359877e-2))
		shll = 0.0;
	CSGP4_DEREF(domdt) = sgs + sghl;
	CSGP4_DEREF(dnodt) = shs;
	if (sinim != 0.0)
	{
		CSGP4_DEREF(domdt) = CSGP4_DEREF(domdt) - cosim / sinim * shll;
		CSGP4_DEREF(dnodt) = CSGP4_DEREF(dnodt) + shll / sinim;
	}

	/* ----------- calculate deep space resonance effects -------- */
	CSGP4_DEREF(dndt) = 0.0;
	theta = FMOD((gsto + tc * rptim), twopi );
	CSGP4_DEREF(em) = CSGP4_DEREF(em) + CSGP4_DEREF(dedt) * t;
	CSGP4_DEREF(inclm) = CSGP4_DEREF(inclm) + CSGP4_DEREF(didt) * t;
	CSGP4_DEREF(argpm) = CSGP4_DEREF(argpm) + CSGP4_DEREF(domdt) * t;
	CSGP4_DEREF(nodem) = CSGP4_DEREF(nodem) + CSGP4_DEREF(dnodt) * t;
	CSGP4_DEREF(mm) = CSGP4_DEREF(mm) + CSGP4_DEREF(dmdt) * t;
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//if (inclm < 0.0)
	//  {
	//	inclm  = -inclm;
	//	argpm  = argpm - Math.PI;
	//	nodem = nodem + Math.PI;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if (CSGP4_DEREF(irez) != 0)
	{
		aonv = POW(CSGP4_DEREF(nm) / xke, x2o3);

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if (CSGP4_DEREF(irez) == 2)
		{
			cosisq = cosim * cosim;
			emo = CSGP4_DEREF(em);
			CSGP4_DEREF(em) = ecco;
			emsqo = emsq;
			emsq = eccsq;
			eoc = CSGP4_DEREF(em) * emsq;
			g201 = -0.306 - (CSGP4_DEREF(em) - 0.64) * 0.440;

			if (CSGP4_DEREF(em) <= 0.65)
			{
				g211 = 3.616 - 13.2470 * CSGP4_DEREF(em) + 16.2900 * emsq;
				g310 = -19.302 + 117.3900 * CSGP4_DEREF(em) - 228.4190 * emsq + 156.5910 * eoc;
				g322 = -18.9068 + 109.7927 * CSGP4_DEREF(em) - 214.6334 * emsq + 146.5816 * eoc;
				g410 = -41.122 + 242.6940 * CSGP4_DEREF(em) - 471.0940 * emsq + 313.9530 * eoc;
				g422 = -146.407 + 841.8800 * CSGP4_DEREF(em) - 1629.014 * emsq + 1083.4350 * eoc;
				g520 = -532.114 + 3017.977 * CSGP4_DEREF(em) - 5740.032 * emsq + 3708.2760 * eoc;
			}
			else
			{
				g211 = -72.099 + 331.819 * CSGP4_DEREF(em) - 508.738 * emsq + 266.724 * eoc;
				g310 = -346.844 + 1582.851 * CSGP4_DEREF(em) - 2415.925 * emsq + 1246.113 * eoc;
				g322 = -342.585 + 1554.908 * CSGP4_DEREF(em) - 2366.899 * emsq + 1215.972 * eoc;
				g410 = -1052.797 + 4758.686 * CSGP4_DEREF(em) - 7193.992 * emsq + 3651.957 * eoc;
				g422 = -3581.690 + 16178.110 * CSGP4_DEREF(em) - 24462.770 * emsq + 12422.520 * eoc;
				if (*em > 0.715)
					g520 = -5149.66 + 29936.92 * CSGP4_DEREF(em) - 54087.36 * emsq + 31324.56 * eoc;
				else
					g520 = 1464.74 - 4664.75 * CSGP4_DEREF(em) + 3763.64 * emsq;
			}
			if (*em < 0.7)
			{
				g533 = -919.22770 + 4988.6100 * CSGP4_DEREF(em) - 9064.7700 * emsq + 5542.21 * eoc;
				g521 = -822.71072 + 4568.6173 * CSGP4_DEREF(em) - 8491.4146 * emsq + 5337.524 * eoc;
				g532 = -853.66600 + 4690.2500 * CSGP4_DEREF(em) - 8624.7700 * emsq + 5341.4 * eoc;
			}
			else
			{
				g533 = -37995.780 + 161616.52 * CSGP4_DEREF(em) - 229838.20 * emsq + 109377.94 * eoc;
				g521 = -51752.104 + 218913.95 * CSGP4_DEREF(em) - 309468.16 * emsq + 146349.42 * eoc;
				g532 = -40023.880 + 170470.89 * CSGP4_DEREF(em) - 242699.48 * emsq + 115605.82 * eoc;
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
			xno2 = CSGP4_DEREF(nm) * CSGP4_DEREF(nm);
			ainv2 = aonv * aonv;
			temp1 = 3.0 * xno2 * ainv2;
			temp = temp1 * root22;
			CSGP4_DEREF(d2201) = temp * f220 * g201;
			CSGP4_DEREF(d2211) = temp * f221 * g211;
			temp1 = temp1 * aonv;
			temp = temp1 * root32;
			CSGP4_DEREF(d3210) = temp * f321 * g310;
			CSGP4_DEREF(d3222) = temp * f322 * g322;
			temp1 = temp1 * aonv;
			temp = 2.0 * temp1 * root44;
			CSGP4_DEREF(d4410) = temp * f441 * g410;
			CSGP4_DEREF(d4422) = temp * f442 * g422;
			temp1 = temp1 * aonv;
			temp = temp1 * root52;
			CSGP4_DEREF(d5220) = temp * f522 * g520;
			CSGP4_DEREF(d5232) = temp * f523 * g532;
			temp = 2.0 * temp1 * root54;
			CSGP4_DEREF(d5421) = temp * f542 * g521;
			CSGP4_DEREF(d5433) = temp * f543 * g533;
			CSGP4_DEREF(xlamo) = FMOD( (mo + nodeo + nodeo - theta - theta), twopi );
			CSGP4_DEREF(xfact) = mdot + CSGP4_DEREF(dmdt) + 2.0 * (nodedot + CSGP4_DEREF(dnodt) - rptim) - no;
			CSGP4_DEREF(em) = emo;
			emsq = emsqo;
		}

		/* ---------------- synchronous resonance terms -------------- */
		if (CSGP4_DEREF(irez) == 1)
		{
			g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
			g310 = 1.0 + 2.0 * emsq;
			g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
			f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
			f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
			f330 = 1.0 + cosim;
			f330 = 1.875 * f330 * f330 * f330;
			CSGP4_DEREF(del1) = 3.0 * CSGP4_DEREF(nm) * CSGP4_DEREF(nm) * aonv * aonv;
			CSGP4_DEREF(del2) = 2.0 * CSGP4_DEREF(del1) * f220 * g200 * q22;
			CSGP4_DEREF(del3) = 3.0 * CSGP4_DEREF(del1) * f330 * g300 * q33 * aonv;
			CSGP4_DEREF(del1) = CSGP4_DEREF(del1) * f311 * g310 * q31 * aonv;
			CSGP4_DEREF(xlamo) = FMOD(mo + nodeo + argpo - theta, twopi );
			CSGP4_DEREF(xfact) = mdot + xpidot - rptim + CSGP4_DEREF(dmdt) + CSGP4_DEREF(domdt) + CSGP4_DEREF(dnodt) - no;
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		CSGP4_DEREF(xli) = CSGP4_DEREF(xlamo);
		CSGP4_DEREF(xni) = no;
		CSGP4_DEREF(atime) = 0.0;
		CSGP4_DEREF(nm) = no + CSGP4_DEREF(dndt);
	}

	//#include "debug3.cpp"
}  // end dsinit


CSGP4_DECORATOR int sgp4init_simple
	 (
	   /*enum gravconsttype whichconst,*/ SGPF epoch,
	   SGPF xbstar, SGPF xndot, SGPF xnddot, SGPF xecco, SGPF xargpo,
	   SGPF xinclo, SGPF xmo, SGPF xno_kozai,
	   SGPF xnodeo, SGPF initial_time, SGPF * r, SGPF * v,
	   SGPF * a_alta_altp
	 )
{
	struct elsetrec_simple gsr = { 0 };

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

#if 0
	/* ----------- set all near earth variables to zero ------------ */
	gsr.isimp = false; gsr.aycof = 0.0;
	gsr.con41 = 0.0; gsr.cc1 = 0.0; gsr.cc4 = 0.0;
	gsr.cc5 = 0.0; gsr.d2 = 0.0; gsr.d3 = 0.0;
	gsr.d4 = 0.0; gsr.delmo = 0.0; gsr.eta = 0.0;
	gsr.argpdot = 0.0; gsr.omgcof = 0.0; gsr.sinmao = 0.0;
	gsr.t = 0.0; gsr.t2cof = 0.0; gsr.t3cof = 0.0;
	gsr.t4cof = 0.0; gsr.t5cof = 0.0; gsr.x1mth2 = 0.0;
	gsr.x7thm1 = 0.0; gsr.mdot = 0.0; gsr.nodedot = 0.0;
	gsr.xlcof = 0.0; gsr.xmcof = 0.0; gsr.nodecf = 0.0;

	/* ----------- set all deep space variables to zero ------------ */
	gsr.irez = 0; gsr.d2201 = 0.0; gsr.d2211 = 0.0;
	gsr.d3210 = 0.0; gsr.d3222 = 0.0; gsr.d4410 = 0.0;
	gsr.d4422 = 0.0; gsr.d5220 = 0.0; gsr.d5232 = 0.0;
	gsr.d5421 = 0.0; gsr.d5433 = 0.0; gsr.dedt = 0.0;
	gsr.del1 = 0.0; gsr.del2 = 0.0; gsr.del3 = 0.0;
	gsr.didt = 0.0; gsr.dmdt = 0.0; gsr.dnodt = 0.0;
	gsr.domdt = 0.0; gsr.e3 = 0.0; gsr.ee2 = 0.0;
	gsr.peo = 0.0; gsr.pgho = 0.0; gsr.pho = 0.0;
	gsr.pinco = 0.0; gsr.plo = 0.0; gsr.se2 = 0.0;
	gsr.se3 = 0.0; gsr.sgh2 = 0.0; gsr.sgh3 = 0.0;
	gsr.sgh4 = 0.0; gsr.sh2 = 0.0; gsr.sh3 = 0.0;
	gsr.si2 = 0.0; gsr.si3 = 0.0; gsr.sl2 = 0.0;
	gsr.sl3 = 0.0; gsr.sl4 = 0.0; gsr.gsto = 0.0;
	gsr.xfact = 0.0; gsr.xgh2 = 0.0; gsr.xgh3 = 0.0;
	gsr.xgh4 = 0.0; gsr.xh2 = 0.0; gsr.xh3 = 0.0;
	gsr.xi2 = 0.0; gsr.xi3 = 0.0; gsr.xl2 = 0.0;
	gsr.xl3 = 0.0; gsr.xl4 = 0.0; gsr.xlamo = 0.0;
	gsr.zmol = 0.0; gsr.zmos = 0.0; gsr.atime = 0.0;
	gsr.xli = 0.0; gsr.xni = 0.0;
#endif

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow aernate values
	// this is now the only call for the constants
//	getgravconst_simple(/*whichconst,*/ CSGP4_REF(gsr.tumin), CSGP4_REF(gsr.mu), CSGP4_REF(gsr.radiusearthkm), CSGP4_REF(gsr.xke),
//				  CSGP4_REF(gsr.j2), CSGP4_REF(gsr.j3), CSGP4_REF(gsr.j4), CSGP4_REF(gsr.j3oj2));
	//Hard-coded WGS72
	(gsr.mu) = 398600.8;			// in km3 / s2
	(gsr.radiusearthkm) = 6378.135;	 // km
	(gsr.xke) = 60.0 / SQRT((gsr.radiusearthkm) * (gsr.radiusearthkm) * (gsr.radiusearthkm) / (gsr.mu));
	(gsr.tumin) = 1.0 / (gsr.xke);
	(gsr.j2) = 0.001082616;
	(gsr.j3) = -0.00000253881;
	(gsr.j4) = -0.00000165597;
	(gsr.j3oj2) = (gsr.j3) / (gsr.j2);


	//-------------------------------------------------------------------------
	gsr.error = 0;
	//gsr.operationmode = opsmode;
	//strcpy( gsr.satnum, satn );

	// sgp4fix - note the following variables are also passed directly via gsr.
	// it is possible to streamline the sgp4init call by deleting the 'x'
	// variables, but the user would need to set the satrec* values first. we
	// include the additional assignments in case twoline2rv is not used.
	gsr.bstar = xbstar;
	// sgp4fix allow additional parameters in the struct
	gsr.ndot = xndot;
	gsr.nddot = xnddot;
	gsr.ecco = xecco;
	gsr.argpo = xargpo;
	gsr.inclo = xinclo;
	gsr.mo = xmo;
	// sgp4fix rename variables to clarify which mean motion is intended
	gsr.no_kozai = xno_kozai;
	gsr.nodeo = xnodeo;

	// single averaged mean elements
	gsr.am = gsr.em = gsr.im = gsr.Om = gsr.mm = gsr.nm = 0.0;

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow alternate values no longer needed
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	ss = 78.0 / gsr.radiusearthkm + 1.0;
	// sgp4fix use multiply for speed instead of Math.POW
	qzms2ttemp = (120.0 - 78.0) / gsr.radiusearthkm;
	qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
	x2o3 = 2.0 / 3.0;

	//gsr.init = 'y';
	gsr.t = 0.0;

	// sgp4fix remove satn as it is not needed in initl
	initl_simple
		(gsr.xke, gsr.j2, gsr.ecco, epoch, gsr.inclo, gsr.no_kozai, CSGP4_REF(gsr.method),
		  CSGP4_REF(ainv), CSGP4_REF(ao), CSGP4_REF(gsr.con41), CSGP4_REF(con42), CSGP4_REF(cosio), CSGP4_REF(cosio2), CSGP4_REF(eccsq), CSGP4_REF(omeosq),
		  CSGP4_REF(posq), CSGP4_REF(rp), CSGP4_REF(rteosq), CSGP4_REF(sinio), CSGP4_REF(gsr.gsto), CSGP4_REF(gsr.no_unkozai) );
	gsr.a = POW(gsr.no_unkozai * gsr.tumin, (-2.0 / 3.0));
	gsr.alta = gsr.a * (1.0 + gsr.ecco) - 1.0;
	gsr.altp = gsr.a * (1.0 - gsr.ecco) - 1.0;
	gsr.error = 0;

	if( a_alta_altp )
	{
		a_alta_altp[0] = gsr.a;
		a_alta_altp[1] = gsr.alta;
		a_alta_altp[2] = gsr.altp;
	}

	// sgp4fix remove this check as it is unnecessary
	// the mrt check in sgp4 handles decaying satellite cases even if the starting
	// condition is below the surface of te earth
	//	 if (rp < 1.0)
	//	   {
	//		 printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
	//		 gsr.error = 5;
	//	   }

	if ((omeosq >= 0.0) || (gsr.no_unkozai >= 0.0))
	{
		gsr.isimp = false;
		if (rp < (220.0 / gsr.radiusearthkm + 1.0))
			gsr.isimp = true;
		sfour = ss;
		qzms24 = qzms2t;
		perige = (rp - 1.0) * gsr.radiusearthkm;

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if (perige < 156.0)
		{
			sfour = perige - 78.0;
			if (perige < 98.0)
				sfour = 20.0;
			// sgp4fix use multiply for speed instead of Math.POW
			qzms24temp = (120.0 - sfour) / gsr.radiusearthkm;
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
			sfour = sfour / gsr.radiusearthkm + 1.0;
		}
		pinvsq = 1.0 / posq;

		tsi = 1.0 / (ao - sfour);
		gsr.eta = ao * gsr.ecco * tsi;
		etasq = gsr.eta * gsr.eta;
		eeta = gsr.ecco * gsr.eta;
		psisq = FABS(1.0 - etasq);
		coef = qzms24 * POW(tsi, 4.0);
		coef1 = coef / POW(psisq, 3.5);
		cc2 = coef1 * gsr.no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
					   (4.0 + etasq)) + 0.375 * gsr.j2 * tsi / psisq * gsr.con41 *
					   (8.0 + 3.0 * etasq * (8.0 + etasq)));
		gsr.cc1 = gsr.bstar * cc2;
		cc3 = 0.0;
		if (gsr.ecco > 1.0e-4)
			cc3 = -2.0 * coef * tsi * gsr.j3oj2 * gsr.no_unkozai * sinio / gsr.ecco;
		gsr.x1mth2 = 1.0 - cosio2;
		gsr.cc4 = 2.0 * gsr.no_unkozai * coef1 * ao * omeosq *
						  (gsr.eta * (2.0 + 0.5 * etasq) + gsr.ecco *
						  (0.5 + 2.0 * etasq) - gsr.j2 * tsi / (ao * psisq) *
						  (-3.0 * gsr.con41 * (1.0 - 2.0 * eeta + etasq *
						  (1.5 - 0.5 * eeta)) + 0.75 * gsr.x1mth2 *
						  (2.0 * etasq - eeta * (1.0 + etasq)) * COS(2.0 * gsr.argpo)));
		gsr.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
					   (etasq + eeta) + eeta * etasq);
		cosio4 = cosio2 * cosio2;
		temp1 = 1.5 * gsr.j2 * pinvsq * gsr.no_unkozai;
		temp2 = 0.5 * temp1 * gsr.j2 * pinvsq;
		temp3 = -0.46875 * gsr.j4 * pinvsq * pinvsq * gsr.no_unkozai;
		gsr.mdot = gsr.no_unkozai + 0.5 * temp1 * rteosq * gsr.con41 + 0.0625 *
						   temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
		gsr.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
							(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
							temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
		xhdot1 = -temp1 * cosio;
		gsr.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
							 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
		xpidot = gsr.argpdot + gsr.nodedot;
		gsr.omgcof = gsr.bstar * cc3 * COS(gsr.argpo);
		gsr.xmcof = 0.0;
		if (gsr.ecco > 1.0e-4)
			gsr.xmcof = -x2o3 * coef * gsr.bstar / eeta;
		gsr.nodecf = 3.5 * omeosq * xhdot1 * gsr.cc1;
		gsr.t2cof = 1.5 * gsr.cc1;
		// sgp4fix for divide by zero with xinco = 180 deg
		if (FABS(cosio + 1.0) > 1.5e-12)
			gsr.xlcof = -0.25 * gsr.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
		else
			gsr.xlcof = -0.25 * gsr.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
		gsr.aycof = -0.5 * gsr.j3oj2 * sinio;
		// sgp4fix use multiply for speed instead of Math.POW
		delmotemp = 1.0 + gsr.eta * COS(gsr.mo);
		gsr.delmo = delmotemp * delmotemp * delmotemp;
		gsr.sinmao = SIN(gsr.mo);
		gsr.x7thm1 = 7.0 * cosio2 - 1.0;

		/* --------------- deep space initialization ------------- */
		if ((2 * SGPPI / gsr.no_unkozai) >= 225.0)
		{
			gsr.method = true;
			gsr.isimp = true;
			tc = 0.0;
			inclm = gsr.inclo;

			//dscom_simple begin.

				SGPF nodep = gsr.nodeo;
				SGPF argpp = gsr.argpo;
				SGPF inclp = gsr.inclo;
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
				(ss1) = 0.0;
				(ss2) = 0.0;
				(ss3) = 0.0;
				(ss4) = 0.0;
				(ss5) = 0.0;
				(ss6) = 0.0;
				(ss7) = 0.0;
				(s1) = 0.0;
				(s2) = 0.0;
				(s3) = 0.0;
				(s4) = 0.0;
				(s5) = 0.0;
				(s6) = 0.0;
				(s7) = 0.0;
				(sz11) = 0.0;
				(sz12) = 0.0;
				(sz13) = 0.0;
				(sz1) = 0.0;
				(sz2) = 0.0;
				(sz3) = 0.0;
				(sz21) = 0.0;
				(sz22) = 0.0;
				(sz23) = 0.0;
				(sz31) = 0.0;
				(sz32) = 0.0;
				(sz33) = 0.0;
				(z13) = 0.0;
				(z21) = 0.0;
				(z1) = 0.0;
				(z2) = 0.0;
				(z3) = 0.0;
				(z11) = 0.0;
				(z12) = 0.0;
				(z31) = 0.0;
				(z21) = 0.0;
				(z22) = 0.0;
				(z23) = 0.0;
				(z32) = 0.0;
				(z33) = 0.0;

				(nm) = gsr.no_unkozai;
				(em) = gsr.ecco;
				(snodm) = SIN(nodep);
				(cnodm) = COS(nodep);
				(sinomm) = SIN(argpp);
				(cosomm) = COS(argpp);
				(sinim) = SIN(inclp);
				(cosim) = COS(inclp);
				(emsq) = (em) * (em);
				betasq = 1.0 - (emsq);
				(rtemsq) = SQRT(betasq);

				/* ----------------- initialize lunar solar terms --------------- */
				(gsr.peo) = 0.0;
				(gsr.pinco) = 0.0;
				(gsr.plo) = 0.0;
				(gsr.pgho) = 0.0;
				(gsr.pho) = 0.0;
				(day) = epoch + 18261.5 + tc / 1440.0;
				xnodce = FMOD(4.5236020 - 9.2422029e-4 * (day), twopi);
				stem = SIN(xnodce);
				ctem = COS(xnodce);
				zcosil = 0.91375164 - 0.03568096 * ctem;
				zsinil = SQRT(1.0 - zcosil * zcosil);
				zsinhl = 0.089683511 * stem / zsinil;
				zcoshl = SQRT(1.0 - zsinhl * zsinhl);
				(gam) = 5.8351514 + 0.0019443680 * (day);
				zx = 0.39785416 * stem / zsinil;
				zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
				zx = ATAN2(zx, zy);
				zx = (gam) + zx - xnodce;
				zcosgl = COS(zx);
				zsingl = SIN(zx);

				/* ------------------------- do solar terms --------------------- */
				zcosg = zcosgs;
				zsing = zsings;
				zcosi = zcosis;
				zsini = zsinis;
				zcosh = (cnodm);
				zsinh = (snodm);
				cc = c1ss;
				xnoi = 1.0 / nm;

				for (lsflg = 1; lsflg <= 2; lsflg++)
				{
					a1 = zcosg * zcosh + zsing * zcosi * zsinh;
					a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
					a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
					a8 = zsing * zsini;
					a9 = zsing * zsinh + zcosg * zcosi * zcosh;
					a10 = zcosg * zsini;
					a2 = (cosim) * a7 + (sinim) * a8;
					a4 = (cosim) * a9 + (sinim) * a10;
					a5 = -(sinim) * a7 + (cosim) * a8;
					a6 = -(sinim) * a9 + (cosim) * a10;

					x1 = a1 * (cosomm) + a2 * (sinomm);
					x2 = a3 * (cosomm) + a4 * (sinomm);
					x3 = -a1 * (sinomm) + a2 * (cosomm);
					x4 = -a3 * (sinomm) + a4 * (cosomm);
					x5 = a5 * (sinomm);
					x6 = a6 * (sinomm);
					x7 = a5 * (cosomm);
					x8 = a6 * (cosomm);

					(z31) = 12.0 * x1 * x1 - 3.0 * x3 * x3;
					(z32) = 24.0 * x1 * x2 - 6.0 * x3 * x4;
					(z33) = 12.0 * x2 * x2 - 3.0 * x4 * x4;
					(z1) = 3.0 * (a1 * a1 + a2 * a2) + z31 * (emsq);
					(z2) = 6.0 * (a1 * a3 + a2 * a4) + z32 * (emsq);
					(z3) = 3.0 * (a3 * a3 + a4 * a4) + z33 * (emsq);
					(z11) = -6.0 * a1 * a5 + (emsq) * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
					(z12) = -6.0 * (a1 * a6 + a3 * a5) + (emsq) *
						   (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
					(z13) = -6.0 * a3 * a6 + (emsq) * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
					(z21) = 6.0 * a2 * a5 + (emsq) * (24.0 * x1 * x5 - 6.0 * x3 * x7);
					(z22) = 6.0 * (a4 * a5 + a2 * a6) + (emsq) *
						   (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
					(z23) = 6.0 * a4 * a6 + (emsq) * (24.0 * x2 * x6 - 6.0 * x4 * x8);
					(z1) = (z1) + (z1) + betasq * (z31);
					(z2) = (z2) + (z2) + betasq * (z32);
					(z3) = (z3) + (z3) + betasq * (z33);
					(s3) = cc * xnoi;
					(s2) = -0.5 * (s3) / (rtemsq);
					(s4) = (s3) * (rtemsq);
					(s1) = -15.0 * (em) * (s4);
					(s5) = x1 * x3 + x2 * x4;
					(s6) = x2 * x3 + x1 * x4;
					(s7) = x2 * x4 - x1 * x3;

					/* ----------------------- do lunar terms ------------------- */
					if (lsflg == 1)
					{
						(ss1) = (s1);
						(ss2) = (s2);
						(ss3) = (s3);
						(ss4) = (s4);
						(ss5) = (s5);
						(ss6) = (s6);
						(ss7) = (s7);
						(sz1) = (z1);
						(sz2) = (z2);
						(sz3) = (z3);
						(sz11) = (z11);
						(sz12) = (z12);
						(sz13) = (z13);
						(sz21) = (z21);
						(sz22) = (z22);
						(sz23) = (z23);
						(sz31) = (z31);
						(sz32) = (z32);
						(sz33) = (z33);
						zcosg = zcosgl;
						zsing = zsingl;
						zcosi = zcosil;
						zsini = zsinil;
						zcosh = zcoshl * (cnodm) + zsinhl * (snodm);
						zsinh = (snodm) * zcoshl - (cnodm) * zsinhl;
						cc = c1l;
					}
				}

				(gsr.zmol) = FMOD( (4.7199672 + 0.22997150 * (day) - (gam)), twopi );
				(gsr.zmos) = FMOD( (6.2565837 + 0.017201977 * (day)), twopi );

				/* ------------------------ do solar terms ---------------------- */
				(gsr.se2) = 2.0 * (ss1) * (ss6);
				(gsr.se3) = 2.0 * (ss1) * (ss7);
				(gsr.si2) = 2.0 * (ss2) * (sz12);
				(gsr.si3) = 2.0 * (ss2) * ((sz13) - (sz11));
				(gsr.sl2) = -2.0 * (ss3) * (sz2);
				(gsr.sl3) = -2.0 * (ss3) * ((sz3) - (sz1));
				(gsr.sl4) = -2.0 * (ss3) * (-21.0 - 9.0 * (emsq)) * zes;
				(gsr.sgh2) = 2.0 * (ss4) * (sz32);
				(gsr.sgh3) = 2.0 * (ss4) * ((sz33) - (sz31));
				(gsr.sgh4) = -18.0 * (ss4) * zes;
				(gsr.sh2) = -2.0 * (ss2) * (sz22);
				(gsr.sh3) = -2.0 * (ss2) * ((sz23) - (sz21));

				/* ------------------------ do lunar terms ---------------------- */
				(gsr.ee2) = 2.0 * (s1) * (s6);
				(gsr.e3) = 2.0 * (s1) * (s7);
				(gsr.xi2) = 2.0 * (s2) * (z12);
				(gsr.xi3) = 2.0 * (s2) * ((z13) - (z11));
				(gsr.xl2) = -2.0 * (s3) * (z2);
				(gsr.xl3) = -2.0 * (s3) * ((z3) - (z1));
				(gsr.xl4) = -2.0 * (s3) * (-21.0 - 9.0 * (emsq)) * zel;
				(gsr.xgh2) = 2.0 * (s4) * (z32);
				(gsr.xgh3) = 2.0 * (s4) * ((z33) - (z31));
				(gsr.xgh4) = -18.0 * (s4) * zel;
				(gsr.xh2) = -2.0 * (s2) * (z22);
				(gsr.xh3) = -2.0 * (s2) * ((z23) - (z21));


			//dscom_simple end.

			// Original code ran dpper_simple( ... true ... ) here. But we don't actually need to.

			argpm = 0.0;
			nodem = 0.0;
			mm = 0.0;

			dsinit_simple
				(
				  gsr.xke,
				  cosim, emsq, gsr.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
				  ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, gsr.t, tc,
				  gsr.gsto, gsr.mo, gsr.mdot, gsr.no_unkozai, gsr.nodeo,
				  gsr.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
				  gsr.ecco, eccsq, CSGP4_REF(em), CSGP4_REF(argpm), CSGP4_REF(inclm), CSGP4_REF(mm), CSGP4_REF(nm), CSGP4_REF(nodem),
				  CSGP4_REF(gsr.irez), CSGP4_REF(gsr.atime),
				  CSGP4_REF(gsr.d2201), CSGP4_REF(gsr.d2211), CSGP4_REF(gsr.d3210), CSGP4_REF(gsr.d3222),
				  CSGP4_REF(gsr.d4410), CSGP4_REF(gsr.d4422), CSGP4_REF(gsr.d5220), CSGP4_REF(gsr.d5232),
				  CSGP4_REF(gsr.d5421), CSGP4_REF(gsr.d5433), CSGP4_REF(gsr.dedt), CSGP4_REF(gsr.didt),
				  CSGP4_REF(gsr.dmdt), CSGP4_REF(dndt), CSGP4_REF(gsr.dnodt), CSGP4_REF(gsr.domdt),
				  CSGP4_REF(gsr.del1), CSGP4_REF(gsr.del2), CSGP4_REF(gsr.del3), CSGP4_REF(gsr.xfact),
				  CSGP4_REF(gsr.xlamo), CSGP4_REF(gsr.xli), CSGP4_REF(gsr.xni)
				);
		}

		/* ----------- set variables if not deep space ----------- */
		if (gsr.isimp != true)
		{
			cc1sq = gsr.cc1 * gsr.cc1;
			gsr.d2 = 4.0 * ao * tsi * cc1sq;
			temp = gsr.d2 * tsi * gsr.cc1 / 3.0;
			gsr.d3 = (17.0 * ao + sfour) * temp;
			gsr.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
							 gsr.cc1;
			gsr.t3cof = gsr.d2 + 2.0 * cc1sq;
			gsr.t4cof = 0.25 * (3.0 * gsr.d3 + gsr.cc1 *
							 (12.0 * gsr.d2 + 10.0 * cc1sq));
			gsr.t5cof = 0.2 * (3.0 * gsr.d4 +
							 12.0 * gsr.cc1 * gsr.d3 +
							 6.0 * gsr.d2 * gsr.d2 +
							 15.0 * cc1sq * (2.0 * gsr.d2 + cc1sq));
		}
	} // if omeosq = 0 ...


	/* finally propagate to zero epoch to initialize all others. */
	// sgp4fix take out check to let satellites process until they are actually below earth surface
	//	   if(gsr.error == 0)
	SGPF r_default[3];
	SGPF v_default[3];
	if( !r ) r = r_default;
	if( !v ) v = v_default;

	//sgp4_simple(satrec, initial_time, initial_r, initial_v);
	// SGP4_SIMPLE

	SGPF tsince = initial_time;

	SGPF am, axnl, aynl, betal, cnod,
		cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu,
		delm, delomg, ecose, el2, eo1,
		ep, esine, argpdf, pl, mrt = 0.0,
		mvt, rdotl, rl, rvdot, rvdotl,
		sin2u, sineo1, sini, sinip, sinsu, sinu,
		snod, su, t2, t3, t4, tem5,
		tempa, tempe, templ, u, ux,
		uy, uz, vx, vy, vz, argpp,
		xinc, xincp, xl, xlm, mp,
		xmdf, xmx, xmy, nodedf, xnode, nodep,
		twopi,  //, j2, j3, tumin, j4, xke, j3oj2, radiusearthkm,
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
	//temp4 = 1.5e-12;
	twopi = 2.0 * SGPPI;
	x2o3 = 2.0 / 3.0;
	// sgp4fix identify constants and allow alternate values
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	vkmpersec = gsr.radiusearthkm * gsr.xke / 60.0;

	/* --------------------- clear sgp4 error flag ----------------- */
	gsr.t = tsince;
	gsr.error = 0;

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf = gsr.mo + gsr.mdot * gsr.t;
	argpdf = gsr.argpo + gsr.argpdot * gsr.t;
	nodedf = gsr.nodeo + gsr.nodedot * gsr.t;
	argpm = argpdf;
	mm = xmdf;
	t2 = gsr.t * gsr.t;
	nodem = nodedf + gsr.nodecf * t2;
	tempa = 1.0 - gsr.cc1 * gsr.t;
	tempe = gsr.bstar * gsr.cc4 * gsr.t;
	templ = gsr.t2cof * t2;

	if (gsr.isimp != true)
	{
		delomg = gsr.omgcof * gsr.t;
		// sgp4fix use mutliply for speed instead of Math.POW
		delmtemp = 1.0 + gsr.eta * COS(xmdf);
		delm = gsr.xmcof *
				 (delmtemp * delmtemp * delmtemp -
				 gsr.delmo);
		temp = delomg + delm;
		mm = xmdf + temp;
		argpm = argpdf - temp;
		t3 = t2 * gsr.t;
		t4 = t3 * gsr.t;
		tempa = tempa - gsr.d2 * t2 - gsr.d3 * t3 -
						 gsr.d4 * t4;
		tempe = tempe + gsr.bstar * gsr.cc5 * (SIN(mm) -
						 gsr.sinmao);
		templ = templ + gsr.t3cof * t3 + t4 * (gsr.t4cof +
						 gsr.t * gsr.t5cof);
	}

	nm = gsr.no_unkozai;
	em = gsr.ecco;
	inclm = gsr.inclo;

	if (gsr.method)
	{
		tc = gsr.t;
/*
		dspace_simple
			(
			  gsr.irez,
			  gsr.d2201, gsr.d2211, gsr.d3210,
			  gsr.d3222, gsr.d4410, gsr.d4422,
			  gsr.d5220, gsr.d5232, gsr.d5421,
			  gsr.d5433, gsr.dedt, gsr.del1,
			  gsr.del2, gsr.del3, gsr.didt,
			  gsr.dmdt, gsr.dnodt, gsr.domdt,
			  gsr.argpo, gsr.argpdot, gsr.t, tc,
			  gsr.gsto, gsr.xfact, gsr.xlamo,
			  gsr.no_unkozai, CSGP4_REF(gsr.atime),
			  CSGP4_REF(em), CSGP4_REF(argpm), CSGP4_REF(inclm), CSGP4_REF(gsr.xli), CSGP4_REF(mm), CSGP4_REF(gsr.xni),
			  CSGP4_REF(nodem), CSGP4_REF(dndt), CSGP4_REF(nm)
			);
*/
		const SGPF twopi = 2.0 * SGPPI;
		//int iretn;  //, iret;
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
		dndt = 0.0;
		theta = FMOD( (gsr.gsto + tc * rptim), twopi );
		em = em + gsr.dedt * gsr.t;

		inclm = inclm + gsr.didt * gsr.t;
		argpm = argpm + gsr.domdt * gsr.t;
		nodem = nodem + gsr.dnodt * gsr.t;
		mm = mm + gsr.dmdt * gsr.t;

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
		if (gsr.irez != 0)
		{
			// sgp4fix streamline check
			if ((gsr.atime == 0.0) || (gsr.t * gsr.atime <= 0.0) || (FABS(gsr.t) < FABS(gsr.atime)))
			{
				gsr.atime = 0.0;
				gsr.xni = gsr.no_unkozai;
				gsr.xli = gsr.xlamo;
			}
			SGPF xli_cache = gsr.xli;
			// sgp4fix move check outside loop
			if (gsr.t > 0.0)
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
				if (gsr.irez != 2)
				{
					xndt = gsr.del1 * SIN(xli_cache - fasx2) + gsr.del2 * SIN(2.0 * (xli_cache - fasx4)) +
							gsr.del3 * SIN(3.0 * (xli_cache - fasx6));
					xldot = gsr.xni + gsr.xfact;
					xnddt = gsr.del1 * COS(xli_cache - fasx2) +
							2.0 * gsr.del2 * COS(2.0 * (xli_cache - fasx4)) +
							3.0 * gsr.del3 * COS(3.0 * (xli_cache - fasx6));
					xnddt = xnddt * xldot;
				}
				else
				{
					/* --------- near - half-day resonance terms -------- */
					xomi = gsr.argpo + gsr.argpdot * gsr.atime;
					x2omi = xomi + xomi;
					x2li = xli_cache + xli_cache;
					xndt = gsr.d2201 * SIN(x2omi + xli_cache - g22) + gsr.d2211 * SIN(xli_cache - g22) +
						  gsr.d3210 * SIN(xomi + xli_cache - g32) + gsr.d3222 * SIN(-xomi + xli_cache - g32) +
						  gsr.d4410 * SIN(x2omi + x2li - g44) + gsr.d4422 * SIN(x2li - g44) +
						  gsr.d5220 * SIN(xomi + xli_cache- g52) + gsr.d5232 * SIN(-xomi + gsr.xli - g52) +
						  gsr.d5421 * SIN(xomi + x2li - g54) + gsr.d5433 * SIN(-xomi + x2li - g54);
					xldot = gsr.xni + gsr.xfact;
					xnddt = gsr.d2201 * COS(x2omi + xli_cache - g22) + gsr.d2211 * COS(xli_cache - g22) +
						  gsr.d3210 * COS(xomi + xli_cache - g32) + gsr.d3222 * COS(-xomi + xli_cache - g32) +
						  gsr.d5220 * COS(xomi + xli_cache - g52) + gsr.d5232 * COS(-xomi + xli_cache - g52) +
						  2.0 * (gsr.d4410 * COS(x2omi + x2li - g44) +
						  gsr.d4422 * COS(x2li - g44) + gsr.d5421 * COS(xomi + x2li - g54) +
						  gsr.d5433 * COS(-xomi + x2li - g54));
					xnddt = xnddt * xldot;
				}

				/* ----------------------- integrator ------------------- */
				// sgp4fix move end checks to end of routine
				if (FABS(gsr.t - gsr.atime) >= stepp)
				{
					//iret = 0;
					iretn = true;
				}
				else // exit here
				{
					ft = gsr.t - gsr.atime;
					iretn = false;
				}

				if (iretn)
				{
					xli_cache = xli_cache + xldot * delt + xndt * step2;
					gsr.xni = gsr.xni + xndt * delt + xnddt * step2;
					gsr.atime = gsr.atime + delt;
				}
			} while( iretn );

			nm = gsr.xni + xndt * ft + xnddt * ft * ft * 0.5;
			xl = xli_cache + xldot * ft + xndt * ft * ft * 0.5;
			gsr.xli = xli_cache;
			if (gsr.irez != 1)
			{
				mm = xl - 2.0 * nodem + 2.0 * theta;
				dndt = nm - gsr.no_unkozai;
			}
			else
			{
				mm = xl - nodem - argpm + theta;
				dndt = nm - gsr.no_unkozai;
			}
			nm = gsr.no_unkozai + dndt;
		}


	} // if method = true

	if (nm <= 0.0)
	{
		//		 printf("# error nm %f\n", nm);
		gsr.error = 2;
		// sgp4fix add return
		//return false;
	}

	am = POW((gsr.xke / nm), x2o3) * tempa * tempa;
	nm = gsr.xke / POW(am, 1.5);
	em = em - tempe;

	// fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/ )
	{
		//		 printf("# error em %f\n", em);
		gsr.error = 1;
		// sgp4fix to return if there is an error in eccentricity
		//return false;
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if (em < 1.0e-6)
		em = 1.0e-6;
	mm = mm + gsr.no_unkozai * templ;
	xlm = mm + argpm + nodem;
	emsq = em * em;
	temp = 1.0 - emsq;

	nodem = FMOD((nodem), twopi );
	argpm = FMOD((argpm), twopi );
	xlm = FMOD((xlm), twopi );
	mm = FMOD((xlm - argpm - nodem), twopi);

	// sgp4fix recover singly averaged mean elements
	gsr.am = am;
	gsr.em = em;
	gsr.im = inclm;
	gsr.Om = nodem;
	gsr.om = argpm;
	gsr.mm = mm;
	gsr.nm = nm;

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
	if (gsr.method)
	{
/*
		dpper_simple
			(
			  gsr.e3, gsr.ee2, gsr.peo,
			  gsr.pgho, gsr.pho, gsr.pinco,
			  gsr.plo, gsr.se2, gsr.se3,
			  gsr.sgh2, gsr.sgh3, gsr.sgh4,
			  gsr.sh2, gsr.sh3, gsr.si2,
			  gsr.si3, gsr.sl2, gsr.sl3,
			  gsr.sl4, gsr.t, gsr.xgh2,
			  gsr.xgh3, gsr.xgh4, gsr.xh2,
			  gsr.xh3, gsr.xi2, gsr.xi3,
			  gsr.xl2, gsr.xl3, gsr.xl4,
			  gsr.zmol, gsr.zmos, gsr.inclo,
			  false, CSGP4_REF(ep), CSGP4_REF(xincp), CSGP4_REF(nodep), CSGP4_REF(argpp), CSGP4_REF(mp)
			);
*/

		//dpper_simple
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
			zm = gsr.zmos + zns * gsr.t;
			// be sure that the initial call has time set to zero
			//if (init)
			//	zm = gsr.zmos;
			zf = zm + 2.0 * zes * SIN(zm);
			sinzf = SIN(zf);
			f2 = 0.5 * sinzf * sinzf - 0.25;
			f3 = -0.5 * sinzf * COS(zf);
			ses = gsr.se2 * f2 + gsr.se3 * f3;
			sis = gsr.si2 * f2 + gsr.si3 * f3;
			sls = gsr.sl2 * f2 + gsr.sl3 * f3 + gsr.sl4 * sinzf;
			sghs = gsr.sgh2 * f2 + gsr.sgh3 * f3 + gsr.sgh4 * sinzf;
			shs = gsr.sh2 * f2 + gsr.sh3 * f3;
			zm = gsr.zmol + znl * gsr.t;
			//if (init)
			//	zm = gsr.zmol;
			zf = zm + 2.0 * zel * SIN(zm);
			sinzf = SIN(zf);
			f2 = 0.5 * sinzf * sinzf - 0.25;
			f3 = -0.5 * sinzf * COS(zf);
			sel = gsr.ee2 * f2 + gsr.e3 * f3;
			sil = gsr.xi2 * f2 + gsr.xi3 * f3;
			sll = gsr.xl2 * f2 + gsr.xl3 * f3 + gsr.xl4 * sinzf;
			sghl = gsr.xgh2 * f2 + gsr.xgh3 * f3 + gsr.xgh4 * sinzf;
			shll = gsr.xh2 * f2 + gsr.xh3 * f3;
			pe = ses + sel;
			pinc = sis + sil;
			pl = sls + sll;
			pgh = sghs + sghl;
			ph = shs + shll;

			// not init
			{
				pe = pe - gsr.peo;
				pinc = pinc - gsr.pinco;
				pl = pl - gsr.plo;
				pgh = pgh - gsr.pgho;
				ph = ph - gsr.pho;
				(xincp) = (xincp) + pinc;
				(ep) = (ep) + pe;
				sinip = SIN((xincp));
				cosip = COS((xincp));

				/* ----------------- apply periodics directly ------------ */
				//  sgp4fix for lyddane choice
				//  strn3 used original inclination - this is technically feasible
				//  gsfc used perturbed inclination - also technically feasible
				//  probably best to readjust the 0.2 limit value and limit discontinuity
				//  0.2 rad = 11.45916 deg
				//  use next line for original strn3 approach and original inclination
				//  if (inclo >= 0.2)
				//  use next line for gsfc version and perturbed inclination
				if ((xincp) >= 0.2)
				{
					ph = ph / sinip;
					pgh = pgh - cosip * ph;
					(argpp) = (argpp) + pgh;
					(nodep) = (nodep) + ph;
					(mp) = (mp) + pl;
				}
				else
				{
					/* ---- apply periodics with lyddane modification ---- */
					sinop = SIN((nodep));
					cosop = COS((nodep));
					alfdp = sinip * sinop;
					betdp = sinip * cosop;
					dalf = ph * cosop + pinc * cosip * sinop;
					dbet = -ph * sinop + pinc * cosip * cosop;
					alfdp = alfdp + dalf;
					betdp = betdp + dbet;
					(nodep) = FMOD( (nodep), twopi );
					//  sgp4fix for afspc written intrinsic functions
					// nodep used without a trigonometric function ahead
					if (((nodep) < 0.0)) // && (opsmode == 'a'))
						(nodep) = (nodep) + twopi;
					xls = (mp) + (argpp) + cosip * (nodep);
					dls = pl + pgh - pinc * (nodep) * sinip;
					xls = xls + dls;
					xnoh = (nodep);
					(nodep) = ATAN2(alfdp, betdp);
					//  sgp4fix for afspc written intrinsic functions
					// nodep used without a trigonometric function ahead
					if (((nodep) < 0.0)) // && (opsmode == 'a'))
						(nodep) = (nodep) + twopi;
					if (FABS(xnoh - (nodep)) > SGPPI)
					{
						if ((nodep) < xnoh)
							(nodep) = (nodep) + twopi;
						else
							(nodep) = (nodep) - twopi;
					}
					(mp) = (mp) + pl;
					(argpp) = xls - (mp) - cosip * (nodep);
				}
			}   // if init == 'n'
		//end dpper_simple


		if (xincp < 0.0)
		{
			xincp = -xincp;
			nodep = nodep + SGPPI;
			argpp = argpp - SGPPI;
		}
		if ((ep < 0.0) || (ep > 1.0))
		{
			//			printf("# error ep %f\n", ep);
			gsr.error = 3;
			// sgp4fix add return
			//return false;
		}

	} // if method = true

	/* -------------------- long period periodics ------------------ */
	if (gsr.method)
	{
		sinip = SIN(xincp);
		cosip = COS(xincp);
		gsr.aycof = -0.5 * gsr.j3oj2 * sinip;
		// sgp4fix for divide by zero for xincp = 180 deg
		if (FABS(cosip + 1.0) > 1.5e-12)
			gsr.xlcof = -0.25 * gsr.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
		else
			gsr.xlcof = -0.25 * gsr.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
	}
	axnl = ep * COS(argpp);
	temp = 1.0 / (am * (1.0 - ep * ep));
	aynl = ep * SIN(argpp) + temp * gsr.aycof;
	xl = mp + argpp + nodep + temp * gsr.xlcof * axnl;

	/* --------------------- solve kepler's equation --------------- */
	u = FMOD( (xl - nodep), twopi );
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
		{
			tem5 = tem5 > 0.0 ? 0.95 : -0.95;
		}
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
		gsr.error = 4;
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
		temp1 = 0.5 * gsr.j2 * temp;
		temp2 = temp1 * temp;

		/* -------------- update for short period periodics ------------ */
		if (gsr.method  ) // == 'd')
		{
			cosisq = cosip * cosip;
			gsr.con41 = 3.0 * cosisq - 1.0;
			gsr.x1mth2 = 1.0 - cosisq;
			gsr.x7thm1 = 7.0 * cosisq - 1.0;
		}
		mrt = rl * (1.0 - 1.5 * temp2 * betal * gsr.con41) +
				0.5 * temp1 * gsr.x1mth2 * cos2u;
		su = su - 0.25 * temp2 * gsr.x7thm1 * sin2u;
		xnode = nodep + 1.5 * temp2 * cosip * sin2u;
		xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
		mvt = rdotl - nm * temp1 * gsr.x1mth2 * sin2u / gsr.xke;
		rvdot = rvdotl + nm * temp1 * (gsr.x1mth2 * cos2u +
				1.5 * gsr.con41) / gsr.xke;

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
		r[0] = (mrt * ux) * gsr.radiusearthkm;
		r[1] = (mrt * uy) * gsr.radiusearthkm;
		r[2] = (mrt * uz) * gsr.radiusearthkm;
		v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
		v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
		v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
	}  // if pl > 0

	// sgp4fix for decaying satellites
	if (mrt < 1.0)
	{
		//		 printf("# decay condition %11.6f \n",mrt);
		gsr.error = 6;
		//return false;
	}

	return gsr.error;
	//#include "debug7.cpp"
	//return true;

	//END SGP4_SIMPLE


	//gsr.init = 'n';

	//#include "debug6.cpp"
	//sgp4fix return boolean. gsr.error contains any error codes
	//return true;
}  // end sgp4init

#endif

