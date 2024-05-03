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

CSGP4_DECORATOR int sgp4init_simple
	 (
	   /*enum gravconsttype whichconst,*/ SGPF epoch,
	   SGPF xbstar, SGPF xndot, SGPF xnddot, SGPF xecco, SGPF xargpo,
	   SGPF xinclo, SGPF xmo, SGPF xno_kozai,
	   SGPF xnodeo, SGPF initial_time, SGPF * r, SGPF * v,
	   SGPF * a_alta_altp
	 )
{
	//struct elsetrec_simple gsr = { 0 };
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
	SGPF am, em, im, Om, mm, nm;
	// sgp4fix add constant parameters to eliminate mutliple calls during execution
	SGPF tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;


	SGPF a /* Not actually used in algo but fun to look at */;
	SGPF no_kozai;
//	SGPF ndot /* Not actually used in algo*/;
//	SGPF nddot /* Not actually used in algo*/;
	SGPF alta /* Not used in algo, but cool anyway to look at */;
	SGPF altp /* Not used in algo, but cool anyway to look at */;

	/* --------------------- local variables ------------------------ */
	SGPF ao, con42, cosio, sinio, cosio2, eccsq,
		 omeosq, posq, rp, rteosq,
		 cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq,
		 cc2, cc3, coef, coef1, cosio4, day, dndt,
		 emsq, eeta, etasq, gam, argpm, nodem,
		 inclm, perige, pinvsq, psisq, qzms24,
		 rtemsq, s1, s2, s3, s4, s5, s6,
		 s7, sfour, ss1, ss2, ss3, ss4, ss5,
		 ss6, ss7, sz1, sz2, sz3, sz11, sz12,
		 sz13, sz21, sz22, sz23, sz31, sz32, sz33,
		 tc, temp, temp1, temp2, temp3, tsi, xpidot,
		 xhdot1, z1, z2, z3, z11, z12, z13,
		 z21, z22, z23, z31, z32, z33,
		 qzms2t, ss,
		 delmotemp, qzms2ttemp, qzms24temp;

	/* ------------------------ initialization --------------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(Math.PI-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	const SGPF temp4 = 1.5e-12;
	const SGPF twopi = 2.0 * SGPPI;
	const SGPF x2o3 = 2.0 / 3.0;

#if 0
	/* ----------- set all near earth variables to zero ------------ */
	isimp = false; aycof = 0.0;
	con41 = 0.0; cc1 = 0.0; cc4 = 0.0;
	cc5 = 0.0; d2 = 0.0; d3 = 0.0;
	d4 = 0.0; delmo = 0.0; eta = 0.0;
	argpdot = 0.0; omgcof = 0.0; siao = 0.0;
	t = 0.0; t2cof = 0.0; t3cof = 0.0;
	t4cof = 0.0; t5cof = 0.0; x1mth2 = 0.0;
	x7thm1 = 0.0; mdot = 0.0; nodedot = 0.0;
	xlcof = 0.0; xmcof = 0.0; nodecf = 0.0;

	/* ----------- set all deep space variables to zero ------------ */
	irez = 0; d2201 = 0.0; d2211 = 0.0;
	d3210 = 0.0; d3222 = 0.0; d4410 = 0.0;
	d4422 = 0.0; d5220 = 0.0; d5232 = 0.0;
	d5421 = 0.0; d5433 = 0.0; dedt = 0.0;
	del1 = 0.0; del2 = 0.0; del3 = 0.0;
	didt = 0.0; dmdt = 0.0; dnodt = 0.0;
	domdt = 0.0; e3 = 0.0; ee2 = 0.0;
	peo = 0.0; pgho = 0.0; pho = 0.0;
	pinco = 0.0; plo = 0.0; se2 = 0.0;
	se3 = 0.0; sgh2 = 0.0; sgh3 = 0.0;
	sgh4 = 0.0; sh2 = 0.0; sh3 = 0.0;
	si2 = 0.0; si3 = 0.0; sl2 = 0.0;
	sl3 = 0.0; sl4 = 0.0; gsto = 0.0;
	xfact = 0.0; xgh2 = 0.0; xgh3 = 0.0;
	xgh4 = 0.0; xh2 = 0.0; xh3 = 0.0;
	xi2 = 0.0; xi3 = 0.0; xl2 = 0.0;
	xl3 = 0.0; xl4 = 0.0; xlamo = 0.0;
	zmol = 0.0; zmos = 0.0; atime = 0.0;
	xli = 0.0; xni = 0.0;
#endif

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow aernate values
	// this is now the only call for the constants
//	getgravconst_simple(/*whichconst,*/ CSGP4_REF(tumin), CSGP4_REF(mu), CSGP4_REF(radiusearthkm), CSGP4_REF(xke),
//				  CSGP4_REF(j2), CSGP4_REF(j3), CSGP4_REF(j4), CSGP4_REF(j3oj2));
	//Hard-coded WGS72
	(mu) = 398600.8;			// in km3 / s2
	(radiusearthkm) = 6378.135;	 // km
	(xke) = 60.0 / SQRT((radiusearthkm) * (radiusearthkm) * (radiusearthkm) / (mu));
	(tumin) = 1.0 / (xke);
	(j2) = 0.001082616;
	(j3) = -0.00000253881;
	(j4) = -0.00000165597;
	(j3oj2) = (j3) / (j2);


	//-------------------------------------------------------------------------
	error = 0;
	//operationmode = opsmode;
	//strcpy( satnum, satn );

	// sgp4fix - note the following variables are also passed directly via 
	// it is possible to streamline the sgp4init call by deleting the 'x'
	// variables, but the user would need to set the satrec* values first. we
	// include the additional assignments in case twoline2rv is not used.
	bstar = xbstar;
	// sgp4fix allow additional parameters in the struct
//	ndot = xndot;
//	nddot = xnddot;
	ecco = xecco;
	argpo = xargpo;
	inclo = xinclo;
	mo = xmo;
	// sgp4fix rename variables to clarify which mean motion is intended
	no_kozai = xno_kozai;
	nodeo = xnodeo;

	// single averaged mean elements
	am = em = im = Om = mm = nm = 0.0;

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow alternate values no longer needed
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	ss = 78.0 / radiusearthkm + 1.0;
	// sgp4fix use multiply for speed instead of Math.POW
	qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
	qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

	//init = 'y';
	t = 0.0;

	// sgp4fix remove satn as it is not needed in initl
	//initl_simple
	//	(xke, j2, ecco, epoch, inclo, no_kozai, CSGP4_REF(method),
	//	  CSGP4_REF(ainv), CSGP4_REF(ao), CSGP4_REF(con41), CSGP4_REF(con42), CSGP4_REF(cosio), CSGP4_REF(cosio2), CSGP4_REF(eccsq), CSGP4_REF(omeosq),
	//	  CSGP4_REF(posq), CSGP4_REF(rp), CSGP4_REF(rteosq), CSGP4_REF(sinio), CSGP4_REF(gsto), CSGP4_REF(no_unkozai) );
	//initl_simple


		/* --------------------- local variables ------------------------ */
		SGPF ak, d1, del, adel, po;

		// sgp4fix use old way of finding gst
		SGPF ds70;
		SGPF ts70, tfrac, c1, thgr70, fk5r, c1p2p;


		/* ------------- calculate auxillary epoch quantities ---------- */
		(eccsq) = ecco * ecco;
		(omeosq) = 1.0 - (eccsq);
		(rteosq) = SQRT((omeosq));
		(cosio) = COS(inclo);
		(cosio2) = (cosio) * (cosio);

		/* ------------------ un-kozai the mean motion ----------------- */
		ak = POW(xke / no_kozai, x2o3);
		d1 = 0.75 * j2 * (3.0 * (cosio2) - 1.0) / ((rteosq) * (omeosq));
		del = d1 / (ak * ak);
		adel = ak * (1.0 - del * del - del *
				(1.0 / 3.0 + 134.0 * del * del / 81.0));
		del = d1 / (adel * adel);
		(no_unkozai) = no_kozai / (1.0 + del);

		(ao) = POW(xke / ((no_unkozai)), x2o3);
		(sinio) = SIN(inclo);
		po = (ao) * (omeosq);
		(con42) = 1.0 - 5.0 * (cosio2);
		(con41) = -(con42) - (cosio2) - (cosio2);
		//(ainv) = 1.0 / (ao);
		(posq) = po * po;
		(rp) = (ao) * (1.0 - ecco);
		(method) = false;
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
			(gsto) = FMOD(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
			if ((gsto) < 0.0)
				(gsto) = (gsto) + twopi;
		}

	// end initl_simple



	a = POW(no_unkozai * tumin, (-2.0 / 3.0));
	alta = a * (1.0 + ecco) - 1.0;
	altp = a * (1.0 - ecco) - 1.0;
	error = 0;

	if( a_alta_altp )
	{
		a_alta_altp[0] = a;
		a_alta_altp[1] = alta;
		a_alta_altp[2] = altp;
	}

	// sgp4fix remove this check as it is unnecessary
	// the mrt check in sgp4 handles decaying satellite cases even if the starting
	// condition is below the surface of te earth
	//	 if (rp < 1.0)
	//	   {
	//		 printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
	//		 error = 5;
	//	   }

	if ((omeosq >= 0.0) || (no_unkozai >= 0.0))
	{
		isimp = false;
		if (rp < (220.0 / radiusearthkm + 1.0))
			isimp = true;
		sfour = ss;
		qzms24 = qzms2t;
		perige = (rp - 1.0) * radiusearthkm;

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if (perige < 156.0)
		{
			sfour = perige - 78.0;
			if (perige < 98.0)
				sfour = 20.0;
			// sgp4fix use multiply for speed instead of Math.POW
			qzms24temp = (120.0 - sfour) / radiusearthkm;
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
			sfour = sfour / radiusearthkm + 1.0;
		}
		pinvsq = 1.0 / posq;

		tsi = 1.0 / (ao - sfour);
		eta = ao * ecco * tsi;
		etasq = eta * eta;
		eeta = ecco * eta;
		psisq = FABS(1.0 - etasq);
		coef = qzms24 * POW(tsi, 4.0);
		coef1 = coef / POW(psisq, 3.5);
		cc2 = coef1 * no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
					   (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * con41 *
					   (8.0 + 3.0 * etasq * (8.0 + etasq)));
		cc1 = bstar * cc2;
		cc3 = 0.0;
		if (ecco > 1.0e-4)
			cc3 = -2.0 * coef * tsi * j3oj2 * no_unkozai * sinio / ecco;
		x1mth2 = 1.0 - cosio2;
		cc4 = 2.0 * no_unkozai * coef1 * ao * omeosq *
						  (eta * (2.0 + 0.5 * etasq) + ecco *
						  (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
						  (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq *
						  (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 *
						  (2.0 * etasq - eeta * (1.0 + etasq)) * COS(2.0 * argpo)));
		cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
					   (etasq + eeta) + eeta * etasq);
		cosio4 = cosio2 * cosio2;
		temp1 = 1.5 * j2 * pinvsq * no_unkozai;
		temp2 = 0.5 * temp1 * j2 * pinvsq;
		temp3 = -0.46875 * j4 * pinvsq * pinvsq * no_unkozai;
		mdot = no_unkozai + 0.5 * temp1 * rteosq * con41 + 0.0625 *
						   temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
		argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
							(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
							temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
		xhdot1 = -temp1 * cosio;
		nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
							 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
		xpidot = argpdot + nodedot;
		omgcof = bstar * cc3 * COS(argpo);
		xmcof = 0.0;
		if (ecco > 1.0e-4)
			xmcof = -x2o3 * coef * bstar / eeta;
		nodecf = 3.5 * omeosq * xhdot1 * cc1;
		t2cof = 1.5 * cc1;
		// sgp4fix for divide by zero with xinco = 180 deg
		if (FABS(cosio + 1.0) > 1.5e-12)
			xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
		else
			xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
		aycof = -0.5 * j3oj2 * sinio;
		// sgp4fix use multiply for speed instead of Math.POW
		delmotemp = 1.0 + eta * COS(mo);
		delmo = delmotemp * delmotemp * delmotemp;
		sinmao = SIN(mo);
		x7thm1 = 7.0 * cosio2 - 1.0;

		/* --------------- deep space initialization ------------- */
		if ((2 * SGPPI / no_unkozai) >= 225.0)
		{
			method = true;
			isimp = true;
			tc = 0.0;
			inclm = inclo;

			//dscom_simple begin.

				SGPF nodep = nodeo;
				SGPF argpp = argpo;
				SGPF inclp = inclo;
				/* -------------------------- constants ------------------------- */
				const SGPF zes = 0.01675;
				const SGPF zel = 0.05490;
				const SGPF c1ss = 2.9864797e-6;
				const SGPF c1l = 4.7968065e-7;
				const SGPF zsinis = 0.39785416;
				const SGPF zcosis = 0.91744867;
				const SGPF zcosgs = 0.1945905;
				const SGPF zsings = -0.98088458;

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

				(nm) = no_unkozai;
				(em) = ecco;
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
				(peo) = 0.0;
				(pinco) = 0.0;
				(plo) = 0.0;
				(pgho) = 0.0;
				(pho) = 0.0;
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

				(zmol) = FMOD( (4.7199672 + 0.22997150 * (day) - (gam)), twopi );
				(zmos) = FMOD( (6.2565837 + 0.017201977 * (day)), twopi );

				/* ------------------------ do solar terms ---------------------- */
				(se2) = 2.0 * (ss1) * (ss6);
				(se3) = 2.0 * (ss1) * (ss7);
				(si2) = 2.0 * (ss2) * (sz12);
				(si3) = 2.0 * (ss2) * ((sz13) - (sz11));
				(sl2) = -2.0 * (ss3) * (sz2);
				(sl3) = -2.0 * (ss3) * ((sz3) - (sz1));
				(sl4) = -2.0 * (ss3) * (-21.0 - 9.0 * (emsq)) * zes;
				(sgh2) = 2.0 * (ss4) * (sz32);
				(sgh3) = 2.0 * (ss4) * ((sz33) - (sz31));
				(sgh4) = -18.0 * (ss4) * zes;
				(sh2) = -2.0 * (ss2) * (sz22);
				(sh3) = -2.0 * (ss2) * ((sz23) - (sz21));

				/* ------------------------ do lunar terms ---------------------- */
				(ee2) = 2.0 * (s1) * (s6);
				(e3) = 2.0 * (s1) * (s7);
				(xi2) = 2.0 * (s2) * (z12);
				(xi3) = 2.0 * (s2) * ((z13) - (z11));
				(xl2) = -2.0 * (s3) * (z2);
				(xl3) = -2.0 * (s3) * ((z3) - (z1));
				(xl4) = -2.0 * (s3) * (-21.0 - 9.0 * (emsq)) * zel;
				(xgh2) = 2.0 * (s4) * (z32);
				(xgh3) = 2.0 * (s4) * ((z33) - (z31));
				(xgh4) = -18.0 * (s4) * zel;
				(xh2) = -2.0 * (s2) * (z22);
				(xh3) = -2.0 * (s2) * ((z23) - (z21));


			//dscom_simple end.

			// Original code ran dpper_simple( ... true ... ) here. But we don't actually need to.

			argpm = 0.0;
			nodem = 0.0;
			mm = 0.0;

/*
			dsinit_simple
				(
				  xke,
				  cosim, emsq, argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
				  ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, t, tc,
				  gsto, mo, mdot, no_unkozai, nodeo,
				  nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
				  ecco, eccsq, CSGP4_REF(em), CSGP4_REF(argpm), CSGP4_REF(inclm), CSGP4_REF(mm), CSGP4_REF(nm), CSGP4_REF(nodem),
				  CSGP4_REF(irez), CSGP4_REF(atime),
				  CSGP4_REF(d2201), CSGP4_REF(d2211), CSGP4_REF(d3210), CSGP4_REF(d3222),
				  CSGP4_REF(d4410), CSGP4_REF(d4422), CSGP4_REF(d5220), CSGP4_REF(d5232),
				  CSGP4_REF(d5421), CSGP4_REF(d5433), CSGP4_REF(dedt), CSGP4_REF(didt),
				  CSGP4_REF(dmdt), CSGP4_REF(dndt), CSGP4_REF(dnodt), CSGP4_REF(domdt),
				  CSGP4_REF(del1), CSGP4_REF(del2), CSGP4_REF(del3), CSGP4_REF(xfact),
				  CSGP4_REF(xlamo), CSGP4_REF(xli), CSGP4_REF(xni)
				);
*/

			//dsinit_simple

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
				(irez) = 0;
				if (((nm) < 0.0052359877) && ((nm) > 0.0034906585))
					(irez) = 1;
				if (((nm) >= 8.26e-3) && ((nm) <= 9.24e-3) && ((em) >= 0.5))
					(irez) = 2;

				/* ------------------------ do solar terms ------------------- */
				ses = ss1 * zns * ss5;
				sis = ss2 * zns * (sz11 + sz13);
				sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
				sghs = ss4 * zns * (sz31 + sz33 - 6.0);
				shs = -zns * ss2 * (sz21 + sz23);
				// sgp4fix for 180 deg incl
				if (((inclm) < 5.2359877e-2) || ((inclm) > SGPPI - 5.2359877e-2))
					shs = 0.0;
				if (sinim != 0.0)
					shs = shs / sinim;
				sgs = sghs - cosim * shs;

				/* ------------------------- do lunar terms ------------------ */
				(dedt) = ses + s1 * znl * s5;
				(didt) = sis + s2 * znl * (z11 + z13);
				(dmdt) = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
				sghl = s4 * znl * (z31 + z33 - 6.0);
				shll = -znl * s2 * (z21 + z23);
				// sgp4fix for 180 deg incl
				if (((inclm) < 5.2359877e-2) || ((inclm) > SGPPI - 5.2359877e-2))
					shll = 0.0;
				(domdt) = sgs + sghl;
				(dnodt) = shs;
				if (sinim != 0.0)
				{
					(domdt) = (domdt) - cosim / sinim * shll;
					(dnodt) = (dnodt) + shll / sinim;
				}

				/* ----------- calculate deep space resonance effects -------- */
				(dndt) = 0.0;
				theta = FMOD((gsto + tc * rptim), twopi );
				(em) = (em) + (dedt) * t;
				(inclm) = (inclm) + (didt) * t;
				(argpm) = (argpm) + (domdt) * t;
				(nodem) = (nodem) + (dnodt) * t;
				(mm) = (mm) + (dmdt) * t;
				//   sgp4fix for negative inclinations
				//   the following if statement should be commented out
				//if (inclm < 0.0)
				//  {
				//	inclm  = -inclm;
				//	argpm  = argpm - Math.PI;
				//	nodem = nodem + Math.PI;
				//  }

				/* -------------- initialize the resonance terms ------------- */
				if ((irez) != 0)
				{
					aonv = POW((nm) / xke, x2o3);

					/* ---------- geopotential resonance for 12 hour orbits ------ */
					if ((irez) == 2)
					{
						cosisq = cosim * cosim;
						emo = (em);
						(em) = ecco;
						emsqo = emsq;
						emsq = eccsq;
						eoc = (em) * emsq;
						g201 = -0.306 - ((em) - 0.64) * 0.440;

						if ((em) <= 0.65)
						{
							g211 = 3.616 - 13.2470 * (em) + 16.2900 * emsq;
							g310 = -19.302 + 117.3900 * (em) - 228.4190 * emsq + 156.5910 * eoc;
							g322 = -18.9068 + 109.7927 * (em) - 214.6334 * emsq + 146.5816 * eoc;
							g410 = -41.122 + 242.6940 * (em) - 471.0940 * emsq + 313.9530 * eoc;
							g422 = -146.407 + 841.8800 * (em) - 1629.014 * emsq + 1083.4350 * eoc;
							g520 = -532.114 + 3017.977 * (em) - 5740.032 * emsq + 3708.2760 * eoc;
						}
						else
						{
							g211 = -72.099 + 331.819 * (em) - 508.738 * emsq + 266.724 * eoc;
							g310 = -346.844 + 1582.851 * (em) - 2415.925 * emsq + 1246.113 * eoc;
							g322 = -342.585 + 1554.908 * (em) - 2366.899 * emsq + 1215.972 * eoc;
							g410 = -1052.797 + 4758.686 * (em) - 7193.992 * emsq + 3651.957 * eoc;
							g422 = -3581.690 + 16178.110 * (em) - 24462.770 * emsq + 12422.520 * eoc;
							if (em > 0.715)
								g520 = -5149.66 + 29936.92 * (em) - 54087.36 * emsq + 31324.56 * eoc;
							else
								g520 = 1464.74 - 4664.75 * (em) + 3763.64 * emsq;
						}
						if (em < 0.7)
						{
							g533 = -919.22770 + 4988.6100 * (em) - 9064.7700 * emsq + 5542.21 * eoc;
							g521 = -822.71072 + 4568.6173 * (em) - 8491.4146 * emsq + 5337.524 * eoc;
							g532 = -853.66600 + 4690.2500 * (em) - 8624.7700 * emsq + 5341.4 * eoc;
						}
						else
						{
							g533 = -37995.780 + 161616.52 * (em) - 229838.20 * emsq + 109377.94 * eoc;
							g521 = -51752.104 + 218913.95 * (em) - 309468.16 * emsq + 146349.42 * eoc;
							g532 = -40023.880 + 170470.89 * (em) - 242699.48 * emsq + 115605.82 * eoc;
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
						xno2 = (nm) * (nm);
						ainv2 = aonv * aonv;
						temp1 = 3.0 * xno2 * ainv2;
						temp = temp1 * root22;
						(d2201) = temp * f220 * g201;
						(d2211) = temp * f221 * g211;
						temp1 = temp1 * aonv;
						temp = temp1 * root32;
						(d3210) = temp * f321 * g310;
						(d3222) = temp * f322 * g322;
						temp1 = temp1 * aonv;
						temp = 2.0 * temp1 * root44;
						(d4410) = temp * f441 * g410;
						(d4422) = temp * f442 * g422;
						temp1 = temp1 * aonv;
						temp = temp1 * root52;
						(d5220) = temp * f522 * g520;
						(d5232) = temp * f523 * g532;
						temp = 2.0 * temp1 * root54;
						(d5421) = temp * f542 * g521;
						(d5433) = temp * f543 * g533;
						(xlamo) = FMOD( (mo + nodeo + nodeo - theta - theta), twopi );
						(xfact) = mdot + (dmdt) + 2.0 * (nodedot + (dnodt) - rptim) - no_unkozai;
						(em) = emo;
						emsq = emsqo;
					}

					/* ---------------- synchronous resonance terms -------------- */
					if ((irez) == 1)
					{
						g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
						g310 = 1.0 + 2.0 * emsq;
						g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
						f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
						f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
						f330 = 1.0 + cosim;
						f330 = 1.875 * f330 * f330 * f330;
						(del1) = 3.0 * (nm) * (nm) * aonv * aonv;
						(del2) = 2.0 * (del1) * f220 * g200 * q22;
						(del3) = 3.0 * (del1) * f330 * g300 * q33 * aonv;
						(del1) = (del1) * f311 * g310 * q31 * aonv;
						(xlamo) = FMOD(mo + nodeo + argpo - theta, twopi );
						(xfact) = mdot + xpidot - rptim + (dmdt) + (domdt) + (dnodt) - no_unkozai;
					}

					/* ------------ for sgp4, initialize the integrator ---------- */
					(xli) = (xlamo);
					(xni) = no_unkozai;
					(atime) = 0.0;
					(nm) = no_unkozai + (dndt);
				}

			//end dsinit_simple
		}

		/* ----------- set variables if not deep space ----------- */
		if (isimp != true)
		{
			cc1sq = cc1 * cc1;
			d2 = 4.0 * ao * tsi * cc1sq;
			temp = d2 * tsi * cc1 / 3.0;
			d3 = (17.0 * ao + sfour) * temp;
			d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
							 cc1;
			t3cof = d2 + 2.0 * cc1sq;
			t4cof = 0.25 * (3.0 * d3 + cc1 *
							 (12.0 * d2 + 10.0 * cc1sq));
			t5cof = 0.2 * (3.0 * d4 +
							 12.0 * cc1 * d3 +
							 6.0 * d2 * d2 +
							 15.0 * cc1sq * (2.0 * d2 + cc1sq));
		}
	} // if omeosq = 0 ...


	/* finally propagate to zero epoch to initialize all others. */
	// sgp4fix take out check to let satellites process until they are actually below earth surface
	//	   if(error == 0)
	SGPF r_default[3];
	SGPF v_default[3];
	if( !r ) r = r_default;
	if( !v ) v = v_default;

	//sgp4_simple(satrec, initial_time, initial_r, initial_v);
	// SGP4_SIMPLE

	SGPF tsince = initial_time;

	SGPF axnl, aynl, betal, cnod,
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
		//twopi,  //, j2, j3, tumin, j4, xke, j3oj2, radiusearthkm,
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
	//twopi = 2.0 * SGPPI;
	//x2o3 = 2.0 / 3.0;
	// sgp4fix identify constants and allow alternate values
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	vkmpersec = radiusearthkm * xke / 60.0;

	/* --------------------- clear sgp4 error flag ----------------- */
	t = tsince;
	error = 0;

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf = mo + mdot * t;
	argpdf = argpo + argpdot * t;
	nodedf = nodeo + nodedot * t;
	argpm = argpdf;
	mm = xmdf;
	t2 = t * t;
	nodem = nodedf + nodecf * t2;
	tempa = 1.0 - cc1 * t;
	tempe = bstar * cc4 * t;
	templ = t2cof * t2;

	if (isimp != true)
	{
		delomg = omgcof * t;
		// sgp4fix use mutliply for speed instead of Math.POW
		delmtemp = 1.0 + eta * COS(xmdf);
		delm = xmcof *
				 (delmtemp * delmtemp * delmtemp -
				 delmo);
		temp = delomg + delm;
		mm = xmdf + temp;
		argpm = argpdf - temp;
		t3 = t2 * t;
		t4 = t3 * t;
		tempa = tempa - d2 * t2 - d3 * t3 -
						 d4 * t4;
		tempe = tempe + bstar * cc5 * (SIN(mm) -
						 sinmao);
		templ = templ + t3cof * t3 + t4 * (t4cof +
						 t * t5cof);
	}

	nm = no_unkozai;
	em = ecco;
	inclm = inclo;

	if (method)
	{
		tc = t;
/*
		dspace_simple
			(
			  irez,
			  d2201, d2211, d3210,
			  d3222, d4410, d4422,
			  d5220, d5232, d5421,
			  d5433, dedt, del1,
			  del2, del3, didt,
			  dmdt, dnodt, domdt,
			  argpo, argpdot, t, tc,
			  gsto, xfact, xlamo,
			  no_unkozai, CSGP4_REF(atime),
			  CSGP4_REF(em), CSGP4_REF(argpm), CSGP4_REF(inclm), CSGP4_REF(xli), CSGP4_REF(mm), CSGP4_REF(xni),
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
		theta = FMOD( (gsto + tc * rptim), twopi );
		em = em + dedt * t;

		inclm = inclm + didt * t;
		argpm = argpm + domdt * t;
		nodem = nodem + dnodt * t;
		mm = mm + dmdt * t;

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
			if ((atime == 0.0) || (t * atime <= 0.0) || (FABS(t) < FABS(atime)))
			{
				atime = 0.0;
				xni = no_unkozai;
				xli = xlamo;
			}
			SGPF xli_cache = xli;
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
					xldot = xni + xfact;
					xnddt = del1 * COS(xli_cache - fasx2) +
							2.0 * del2 * COS(2.0 * (xli_cache - fasx4)) +
							3.0 * del3 * COS(3.0 * (xli_cache - fasx6));
					xnddt = xnddt * xldot;
				}
				else
				{
					/* --------- near - half-day resonance terms -------- */
					xomi = argpo + argpdot * atime;
					x2omi = xomi + xomi;
					x2li = xli_cache + xli_cache;
					xndt = d2201 * SIN(x2omi + xli_cache - g22) + d2211 * SIN(xli_cache - g22) +
						  d3210 * SIN(xomi + xli_cache - g32) + d3222 * SIN(-xomi + xli_cache - g32) +
						  d4410 * SIN(x2omi + x2li - g44) + d4422 * SIN(x2li - g44) +
						  d5220 * SIN(xomi + xli_cache- g52) + d5232 * SIN(-xomi + xli - g52) +
						  d5421 * SIN(xomi + x2li - g54) + d5433 * SIN(-xomi + x2li - g54);
					xldot = xni + xfact;
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
				if (FABS(t - atime) >= stepp)
				{
					//iret = 0;
					iretn = true;
				}
				else // exit here
				{
					ft = t - atime;
					iretn = false;
				}

				if (iretn)
				{
					xli_cache = xli_cache + xldot * delt + xndt * step2;
					xni = xni + xndt * delt + xnddt * step2;
					atime = atime + delt;
				}
			} while( iretn );

			nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
			xl = xli_cache + xldot * ft + xndt * ft * ft * 0.5;
			xli = xli_cache;
			if (irez != 1)
			{
				mm = xl - 2.0 * nodem + 2.0 * theta;
				dndt = nm - no_unkozai;
			}
			else
			{
				mm = xl - nodem - argpm + theta;
				dndt = nm - no_unkozai;
			}
			nm = no_unkozai + dndt;
		}


	} // if method = true

	if (nm <= 0.0)
	{
		//		 printf("# error nm %f\n", nm);
		error = 2;
		// sgp4fix add return
		//return false;
	}

	am = POW((xke / nm), x2o3) * tempa * tempa;
	nm = xke / POW(am, 1.5);
	em = em - tempe;

	// fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/ )
	{
		//		 printf("# error em %f\n", em);
		error = 1;
		// sgp4fix to return if there is an error in eccentricity
		//return false;
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if (em < 1.0e-6)
		em = 1.0e-6;
	mm = mm + no_unkozai * templ;
	xlm = mm + argpm + nodem;
	emsq = em * em;
	temp = 1.0 - emsq;

	nodem = FMOD((nodem), twopi );
	argpm = FMOD((argpm), twopi );
	xlm = FMOD((xlm), twopi );
	mm = FMOD((xlm - argpm - nodem), twopi);

	// sgp4fix recover singly averaged mean elements
	am = am;
	em = em;
	im = inclm;
	Om = nodem;
//	om = argpm;
	mm = mm;
	nm = nm;

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
	if (method)
	{
/*
		dpper_simple
			(
			  e3, ee2, peo,
			  pgho, pho, pinco,
			  plo, se2, se3,
			  sgh2, sgh3, sgh4,
			  sh2, sh3, si2,
			  si3, sl2, sl3,
			  sl4, t, xgh2,
			  xgh3, xgh4, xh2,
			  xh3, xi2, xi3,
			  xl2, xl3, xl4,
			  zmol, zmos, inclo,
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
			zm = zmos + zns * t;
			// be sure that the initial call has time set to zero
			//if (init)
			//	zm = zmos;
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
			//if (init)
			//	zm = zmol;
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

			// not init
			{
				pe = pe - peo;
				pinc = pinc - pinco;
				pl = pl - plo;
				pgh = pgh - pgho;
				ph = ph - pho;
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
			error = 3;
			// sgp4fix add return
			//return false;
		}

	} // if method = true

	/* -------------------- long period periodics ------------------ */
	if (method)
	{
		sinip = SIN(xincp);
		cosip = COS(xincp);
		aycof = -0.5 * j3oj2 * sinip;
		// sgp4fix for divide by zero for xincp = 180 deg
		if (FABS(cosip + 1.0) > 1.5e-12)
			xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
		else
			xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
	}
	axnl = ep * COS(argpp);
	temp = 1.0 / (am * (1.0 - ep * ep));
	aynl = ep * SIN(argpp) + temp * aycof;
	xl = mp + argpp + nodep + temp * xlcof * axnl;

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
		error = 4;
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
		temp1 = 0.5 * j2 * temp;
		temp2 = temp1 * temp;

		/* -------------- update for short period periodics ------------ */
		if (method  ) // == 'd')
		{
			cosisq = cosip * cosip;
			con41 = 3.0 * cosisq - 1.0;
			x1mth2 = 1.0 - cosisq;
			x7thm1 = 7.0 * cosisq - 1.0;
		}
		mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) +
				0.5 * temp1 * x1mth2 * cos2u;
		su = su - 0.25 * temp2 * x7thm1 * sin2u;
		xnode = nodep + 1.5 * temp2 * cosip * sin2u;
		xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
		mvt = rdotl - nm * temp1 * x1mth2 * sin2u / xke;
		rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u +
				1.5 * con41) / xke;

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
		r[0] = (mrt * ux) * radiusearthkm;
		r[1] = (mrt * uy) * radiusearthkm;
		r[2] = (mrt * uz) * radiusearthkm;
		v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
		v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
		v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
	}  // if pl > 0

	// sgp4fix for decaying satellites
	if (mrt < 1.0)
	{
		//		 printf("# decay condition %11.6f \n",mrt);
		error = 6;
		//return false;
	}

	return error;
	//#include "debug7.cpp"
	//return true;

	//END SGP4_SIMPLE


	//init = 'n';

	//#include "debug6.cpp"
	//sgp4fix return boolean. error contains any error codes
	//return true;
}  // end sgp4init

#endif

