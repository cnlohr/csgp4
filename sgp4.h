#ifndef _SGP4_H
#define _SGP4_H

#include <math.h>

// Single-file C header version of david vallado's SGP4Lib.cs, version: "SGP4 Version 2020-03-12";

enum gravconsttype { wgs72old, wgs72, wgs84 }; // wgs72 is the standard and should be used with JSPOC TLEs

#define SGPPI 3.1415926535897932384626433

struct elsetrec
{
    // change to scripting to handle upcoming numbering changes, either alpha 5 or 9 digit
    char satnum[32];
    int epochyr, epochtynumrev;
    int error;
    char operationmode;
    char init, method;

    /* Near Earth */
    int isimp;
    double aycof, con41, cc1, cc4, cc5, d2, d3, d4,
           delmo, eta, argpdot, omgcof, sinmao, t, t2cof, t3cof,
           t4cof, t5cof, x1mth2, x7thm1, mdot, nodedot, xlcof, xmcof,
           nodecf;

    /* Deep Space */
    int irez;
    double d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232,
           d5421, d5433, dedt, del1, del2, del3, didt, dmdt,
           dnodt, domdt, e3, ee2, peo, pgho, pho, pinco,
           plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3,
           si2, si3, sl2, sl3, sl4, gsto, xfact, xgh2,
           xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3,
           xl4, xlamo, zmol, zmos, atime, xli, xni;

    double a, altp, alta, epochdays, jdsatepoch, jdsatepochF, nddot, ndot,
           bstar, rcse, inclo, nodeo, ecco, argpo, mo,
           no_kozai;
    // sgp4fix add new variables from tle
    char classification;
    char intldesg[32];  //10
    int ephtype;
    long elnum, revnum;
    // sgp4fix add unkozai'd variable
    double no_unkozai;
    // sgp4fix add singly averaged variables
    double am, em, im, Om, om, mm, nm;
    // sgp4fix add constant parameters to eliminate mutliple calls during execution
    double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;

      //       Additional elements to capture relevant TLE and object information:        |
    long dia_mm; // RSO dia in mm
    double period_sec; // Period in seconds
    char active[32]; // "Active S/C" flag (0=n, 1=y) 
    char not_orbital[32]; // "Orbiting S/C" flag (0=n, 1=y)  
    double rcs_m2; // "RCS (m^2)" storage 
};


/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satn        - satellite number - not needed, placed in satrec
*    xke         - reciprocal of tumin                       
*    j2          - j2 zonal harmonic
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 Math.Cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst- no longer used
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

static void initl
     (
    // sgp4fix satn not needed. include in satrec in case needed later  
    // int satn,      
    // sgp4fix just pass in xke and j2
    // gravconsttype whichconst, 
       double xke, double j2,
       double ecco, double epoch, double inclo, double no_kozai,

	// Output params:
       char * method,
       double * ainv, double * ao, double * con41, double * con42, double * cosio,
       double * cosio2, double * eccsq, double * omeosq, double * posq,
       double * rp, double * rteosq, double * sinio, double * gsto,
       char opsmode, double * no_unkozai
     )
{
    /* --------------------- local variables ------------------------ */
    double ak, d1, del, adel, po, x2o3;

    // sgp4fix use old way of finding gst
    double ds70;
    double ts70, tfrac, c1, thgr70, fk5r, c1p2p;
    const double twopi = 2.0 * SGPPI;

    /* ----------------------- earth constants ---------------------- */
    // sgp4fix identify constants and allow alternate values
    // only xke and j2 are used here so pass them in directly
    // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
    x2o3 = 2.0 / 3.0;

    /* ------------- calculate auxillary epoch quantities ---------- */
    *eccsq = ecco * ecco;
    *omeosq = 1.0 - *eccsq;
    *rteosq = sqrt(*omeosq);
    *cosio = cos(inclo);
    *cosio2 = *cosio * *cosio;

    /* ------------------ un-kozai the mean motion ----------------- */
    ak = pow(xke / no_kozai, x2o3);
    d1 = 0.75 * j2 * (3.0 * *cosio2 - 1.0) / (*rteosq * *omeosq);
    del = d1 / (ak * ak);
    adel = ak * (1.0 - del * del - del *
            (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    *no_unkozai = no_kozai / (1.0 + del);

    *ao = pow(xke / (*no_unkozai), x2o3);
    *sinio = sin(inclo);
    po = *ao * *omeosq;
    *con42 = 1.0 - 5.0 * *cosio2;
    *con41 = -*con42 - *cosio2 - *cosio2;
    *ainv = 1.0 / *ao;
    *posq = po * po;
    *rp = *ao * (1.0 - ecco);
    *method = 'n';

    // sgp4fix modern approach to finding sidereal time
    if (opsmode == 'a')
    {
        // sgp4fix use old way of finding gst
        // count integer number of days from 0 jan 1970
        ts70 = epoch - 7305.0;
        ds70 = floor(ts70 + 1.0e-8);
        tfrac = ts70 - ds70;
        // find greenwich location at epoch
        c1 = 1.72027916940703639e-2;
        thgr70 = 1.7321343856509374;
        fk5r = 5.07551419432269442e-15;
        c1p2p = c1 + twopi;
        *gsto = fmod(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
        if (*gsto < 0.0)
            *gsto = *gsto + twopi;
    }
    else
        *gsto = gstime(epoch + 2433281.5);


    //#include "debug5.cpp"
}  // end initl



        /*-----------------------------------------------------------------------------
        *
        *                             procedure sgp4init
        *
        *  this procedure initializes variables for sgp4.
        *
        *  author        : david vallado                  719-573-2600   28 jun 2005
        *
        *  inputs        :
        *    opsmode     - mode of operation afspc or improved 'a', 'i'
        *    whichconst  - which set of constants to use  72, 84
        *    satn        - satellite number
        *    bstar       - sgp4 type drag coefficient              kg/m2er
        *    ecco        - eccentricity
        *    epoch       - epoch time in days from jan 0, 1950. 0 hr
        *    argpo       - argument of perigee (output if ds)
        *    inclo       - inclination
        *    mo          - mean anomaly (output if ds)
        *    no          - mean motion
        *    nodeo       - right ascension of ascending node
        *
        *  outputs       :
        *    satrec      - common values for subsequent calls
        *    return code - non-zero on error.
        *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
        *                   2 - mean motion less than 0.0
        *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
        *                   4 - semi-latus rectum < 0.0
        *                   5 - epoch elements are sub-orbital
        *                   6 - satellite has decayed
        *
        *  locals        :
        *    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
        *    cc1sq  , cc2    , cc3
        *    coef   , coef1
        *    cosio4      -
        *    day         -
        *    dndt        -
        *    em          - eccentricity
        *    emsq        - eccentricity squared
        *    eeta        -
        *    etasq       -
        *    gam         -
        *    argpm       - argument of perigee
        *    nodem       -
        *    inclm       - inclination
        *    mm          - mean anomaly
        *    nm          - mean motion
        *    perige      - perigee
        *    pinvsq      -
        *    psisq       -
        *    qzms24      -
        *    rtemsq      -
        *    s1, s2, s3, s4, s5, s6, s7          -
        *    sfour       -
        *    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
        *    sz1, sz2, sz3
        *    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
        *    tc          -
        *    temp        -
        *    temp1, temp2, temp3       -
        *    tsi         -
        *    xpidot      -
        *    xhdot1      -
        *    z1, z2, z3          -
        *    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
        *
        *  coupling      :
        *    getgravconst-
        *    initl       -
        *    dscom       -
        *    dpper       -
        *    dsinit      -
        *    sgp4        -
        *
        *  references    :
        *    hoots, roehrich, norad spacetrack report #3 1980
        *    hoots, norad spacetrack report #6 1986
        *    hoots, schumacher and glover 2004
        *    vallado, crawford, hujsak, kelso  2006
          ----------------------------------------------------------------------------*/

        public void sgp4init
             (
               gravconsttype whichconst, char opsmode, string satn, double epoch,
               double xbstar, double xndot, double xnddot, double xecco, double xargpo,
               double xinclo, double xmo, double xno_kozai,
               double xnodeo, elsetrec * satrec
             )
        {
            /* --------------------- local variables ------------------------ */
            double ao, ainv, con42, cosio, sinio, cosio2, eccsq,
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
            double r[3];
            double v[3];

            /* ------------------------ initialization --------------------- */
            // sgp4fix divisor for divide by zero check on inclination
            // the old check used 1.0 + Math.Cos(Math.PI-1.0e-9), but then compared it to
            // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
            const double temp4 = 1.5e-12;

            /* ----------- set all near earth variables to zero ------------ */
            satrec.isimp = 0; satrec.method = 'n'; satrec.aycof = 0.0;
            satrec.con41 = 0.0; satrec.cc1 = 0.0; satrec.cc4 = 0.0;
            satrec.cc5 = 0.0; satrec.d2 = 0.0; satrec.d3 = 0.0;
            satrec.d4 = 0.0; satrec.delmo = 0.0; satrec.eta = 0.0;
            satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao = 0.0;
            satrec.t = 0.0; satrec.t2cof = 0.0; satrec.t3cof = 0.0;
            satrec.t4cof = 0.0; satrec.t5cof = 0.0; satrec.x1mth2 = 0.0;
            satrec.x7thm1 = 0.0; satrec.mdot = 0.0; satrec.nodedot = 0.0;
            satrec.xlcof = 0.0; satrec.xmcof = 0.0; satrec.nodecf = 0.0;

            /* ----------- set all deep space variables to zero ------------ */
            satrec.irez = 0; satrec.d2201 = 0.0; satrec.d2211 = 0.0;
            satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
            satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
            satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt = 0.0;
            satrec.del1 = 0.0; satrec.del2 = 0.0; satrec.del3 = 0.0;
            satrec.didt = 0.0; satrec.dmdt = 0.0; satrec.dnodt = 0.0;
            satrec.domdt = 0.0; satrec.e3 = 0.0; satrec.ee2 = 0.0;
            satrec.peo = 0.0; satrec.pgho = 0.0; satrec.pho = 0.0;
            satrec.pinco = 0.0; satrec.plo = 0.0; satrec.se2 = 0.0;
            satrec.se3 = 0.0; satrec.sgh2 = 0.0; satrec.sgh3 = 0.0;
            satrec.sgh4 = 0.0; satrec.sh2 = 0.0; satrec.sh3 = 0.0;
            satrec.si2 = 0.0; satrec.si3 = 0.0; satrec.sl2 = 0.0;
            satrec.sl3 = 0.0; satrec.sl4 = 0.0; satrec.gsto = 0.0;
            satrec.xfact = 0.0; satrec.xgh2 = 0.0; satrec.xgh3 = 0.0;
            satrec.xgh4 = 0.0; satrec.xh2 = 0.0; satrec.xh3 = 0.0;
            satrec.xi2 = 0.0; satrec.xi3 = 0.0; satrec.xl2 = 0.0;
            satrec.xl3 = 0.0; satrec.xl4 = 0.0; satrec.xlamo = 0.0;
            satrec.zmol = 0.0; satrec.zmos = 0.0; satrec.atime = 0.0;
            satrec.xli = 0.0; satrec.xni = 0.0;

            /* ------------------------ earth constants ----------------------- */
            // sgp4fix identify constants and allow alternate values
            // this is now the only call for the constants
            getgravconst(whichconst, out satrec.tumin, out satrec.mu, out satrec.radiusearthkm, out satrec.xke,
                          out satrec.j2, out satrec.j3, out satrec.j4, out satrec.j3oj2);
            //-------------------------------------------------------------------------
            satrec.error = 0;
            satrec.operationmode = opsmode;
            satrec.satnum = satn;

            // sgp4fix - note the following variables are also passed directly via satrec.
            // it is possible to streamline the sgp4init call by deleting the 'x'
            // variables, but the user would need to set the satrec.* values first. we
            // include the additional assignments in case twoline2rv is not used.
            satrec.bstar = xbstar;
            // sgp4fix allow additional parameters in the struct
            satrec.ndot = xndot;
            satrec.nddot = xnddot;
            satrec.ecco = xecco;
            satrec.argpo = xargpo;
            satrec.inclo = xinclo;
            satrec.mo = xmo;
            // sgp4fix rename variables to clarify which mean motion is intended
            satrec.no_kozai = xno_kozai;
            satrec.nodeo = xnodeo;

            // single averaged mean elements
            satrec.am = satrec.em = satrec.im = satrec.Om = satrec.mm = satrec.nm = 0.0;

            /* ------------------------ earth constants ----------------------- */
            // sgp4fix identify constants and allow alternate values no longer needed
            // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
            ss = 78.0 / satrec.radiusearthkm + 1.0;
            // sgp4fix use multiply for speed instead of Math.Pow
            qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
            qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
            x2o3 = 2.0 / 3.0;

            satrec.init = 'y';
            satrec.t = 0.0;

            // sgp4fix remove satn as it is not needed in initl
            initl
                (satrec.xke, satrec.j2, satrec.ecco, epoch, satrec.inclo, ref satrec.no_kozai, out satrec.method,
                  &ainv, &ao, &satrec.con41, &con42, &cosio, &cosio2, &eccsq, &omeosq,
                  &posq, &rp, &rteosq, &sinio, &satrec.gsto, satrec.operationmode, &satrec.no_unkozai);
            satrec.a = Math.Pow(satrec.no_unkozai * satrec.tumin, (-2.0 / 3.0));
            satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
            satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
            satrec.error = 0;

            // sgp4fix remove this check as it is unnecessary
            // the mrt check in sgp4 handles decaying satellite cases even if the starting
            // condition is below the surface of te earth
            //     if (rp < 1.0)
            //       {
            //         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
            //         satrec.error = 5;
            //       }

            if ((omeosq >= 0.0) || (satrec.no_unkozai >= 0.0))
            {
                satrec.isimp = 0;
                if (rp < (220.0 / satrec.radiusearthkm + 1.0))
                    satrec.isimp = 1;
                sfour = ss;
                qzms24 = qzms2t;
                perige = (rp - 1.0) * satrec.radiusearthkm;

                /* - for perigees below 156 km, s and qoms2t are altered - */
                if (perige < 156.0)
                {
                    sfour = perige - 78.0;
                    if (perige < 98.0)
                        sfour = 20.0;
                    // sgp4fix use multiply for speed instead of Math.Pow
                    qzms24temp = (120.0 - sfour) / satrec.radiusearthkm;
                    qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                    sfour = sfour / satrec.radiusearthkm + 1.0;
                }
                pinvsq = 1.0 / posq;

                tsi = 1.0 / (ao - sfour);
                satrec.eta = ao * satrec.ecco * tsi;
                etasq = satrec.eta * satrec.eta;
                eeta = satrec.ecco * satrec.eta;
                psisq = Math.Abs(1.0 - etasq);
                coef = qzms24 * Math.Pow(tsi, 4.0);
                coef1 = coef / Math.Pow(psisq, 3.5);
                cc2 = coef1 * satrec.no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta *
                               (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 *
                               (8.0 + 3.0 * etasq * (8.0 + etasq)));
                satrec.cc1 = satrec.bstar * cc2;
                cc3 = 0.0;
                if (satrec.ecco > 1.0e-4)
                    cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no_unkozai * sinio / satrec.ecco;
                satrec.x1mth2 = 1.0 - cosio2;
                satrec.cc4 = 2.0 * satrec.no_unkozai * coef1 * ao * omeosq *
                                  (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
                                  (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (ao * psisq) *
                                  (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                                  (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                                  (2.0 * etasq - eeta * (1.0 + etasq)) * Math.Cos(2.0 * satrec.argpo)));
                satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                               (etasq + eeta) + eeta * etasq);
                cosio4 = cosio2 * cosio2;
                temp1 = 1.5 * satrec.j2 * pinvsq * satrec.no_unkozai;
                temp2 = 0.5 * temp1 * satrec.j2 * pinvsq;
                temp3 = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no_unkozai;
                satrec.mdot = satrec.no_unkozai + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 *
                                   temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
                satrec.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
                                    (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                                    temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
                xhdot1 = -temp1 * cosio;
                satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                                     2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
                xpidot = satrec.argpdot + satrec.nodedot;
                satrec.omgcof = satrec.bstar * cc3 * Math.Cos(satrec.argpo);
                satrec.xmcof = 0.0;
                if (satrec.ecco > 1.0e-4)
                    satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
                satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
                satrec.t2cof = 1.5 * satrec.cc1;
                // sgp4fix for divide by zero with xinco = 180 deg
                if (Math.Abs(cosio + 1.0) > 1.5e-12)
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
                else
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
                satrec.aycof = -0.5 * satrec.j3oj2 * sinio;
                // sgp4fix use multiply for speed instead of Math.Pow
                delmotemp = 1.0 + satrec.eta * Math.Cos(satrec.mo);
                satrec.delmo = delmotemp * delmotemp * delmotemp;
                satrec.sinmao = Math.Sin(satrec.mo);
                satrec.x7thm1 = 7.0 * cosio2 - 1.0;

                /* --------------- deep space initialization ------------- */
                if ((2 * SGPPI / satrec.no_unkozai) >= 225.0)
                {
                    satrec.method = 'd';
                    satrec.isimp = 1;
                    tc = 0.0;
                    inclm = satrec.inclo;

                    dscom
                        (
                          epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo,
                          satrec.no_unkozai, out snodm, out cnodm, out sinim, out cosim, out sinomm, out cosomm,
                          out day, out satrec.e3, out satrec.ee2, out em, out emsq, out gam,
                          out satrec.peo, out satrec.pgho, out satrec.pho, out satrec.pinco,
                          out satrec.plo, out rtemsq, out satrec.se2, out satrec.se3,
                          out satrec.sgh2, out satrec.sgh3, out satrec.sgh4,
                          out satrec.sh2, out  satrec.sh3, out satrec.si2, out satrec.si3,
                          out satrec.sl2, out satrec.sl3, out satrec.sl4, out s1, out s2, out s3, out s4, out s5,
                          out s6, out s7, out ss1, out ss2, out ss3, out ss4, out ss5, out ss6, out ss7, out sz1, out sz2, out sz3,
                          out sz11, out sz12, out sz13, out sz21, out sz22, out sz23, out sz31, out sz32, out sz33,
                          out satrec.xgh2, out satrec.xgh3, out satrec.xgh4, out satrec.xh2,
                          out satrec.xh3, out satrec.xi2, out satrec.xi3, out satrec.xl2,
                          out satrec.xl3, out satrec.xl4, out nm, out z1, out z2, out z3, out z11,
                          out z12, out z13, out z21, out z22, out z23, out z31, out z32, out z33,
                          out satrec.zmol, out satrec.zmos
                        );
                    dpper
                        (
                          satrec.e3, satrec.ee2, satrec.peo, satrec.pgho,
                          satrec.pho, satrec.pinco, satrec.plo, satrec.se2,
                          satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
                          satrec.sh2, satrec.sh3, satrec.si2, satrec.si3,
                          satrec.sl2, satrec.sl3, satrec.sl4, satrec.t,
                          satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
                          satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
                          satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos, inclm, satrec.init,
                          ref satrec.ecco, ref satrec.inclo, ref satrec.nodeo, ref satrec.argpo, ref satrec.mo,
                          satrec.operationmode
                        );

                    argpm = 0.0;
                    nodem = 0.0;
                    mm = 0.0;

                    dsinit
                        (
                          satrec.xke,
                          cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                          ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
                          satrec.gsto, satrec.mo, satrec.mdot, satrec.no_unkozai, satrec.nodeo,
                          satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                          satrec.ecco, eccsq, ref em, ref argpm, ref inclm, ref mm, ref nm, ref nodem,
                          out satrec.irez, ref satrec.atime,
                          ref satrec.d2201, ref satrec.d2211, ref satrec.d3210, ref satrec.d3222,
                          ref satrec.d4410, ref satrec.d4422, ref satrec.d5220, ref satrec.d5232,
                          ref satrec.d5421, ref satrec.d5433, out satrec.dedt, out satrec.didt,
                          out satrec.dmdt, out dndt, out satrec.dnodt, out satrec.domdt,
                          ref satrec.del1, ref satrec.del2, ref satrec.del3, ref satrec.xfact,
                          ref satrec.xlamo, ref satrec.xli, ref satrec.xni
                        );
                }

                /* ----------- set variables if not deep space ----------- */
                if (satrec.isimp != 1)
                {
                    cc1sq = satrec.cc1 * satrec.cc1;
                    satrec.d2 = 4.0 * ao * tsi * cc1sq;
                    temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
                    satrec.d3 = (17.0 * ao + sfour) * temp;
                    satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
                                     satrec.cc1;
                    satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
                    satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
                                     (12.0 * satrec.d2 + 10.0 * cc1sq));
                    satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
                                     12.0 * satrec.cc1 * satrec.d3 +
                                     6.0 * satrec.d2 * satrec.d2 +
                                     15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
                }
            } // if omeosq = 0 ...

            /* finally propagate to zero epoch to initialize all others. */
            // sgp4fix take out check to let satellites process until they are actually below earth surface
            //       if(satrec.error == 0)
            sgp4(ref satrec, 0.0, r, v);

            satrec.init = 'n';

            //#include "debug6.cpp"
            //sgp4fix return boolean. satrec.error contains any error codes
            //return true;
        }  // end sgp4init




        /*-----------------------------------------------------------------------------
        *
        *                             procedure sgp4
        *
        *  this procedure is the sgp4 prediction model from space command. this is an
        *    updated and combined version of sgp4 and sdp4, which were originally
        *    published separately in spacetrack report #3. this version follows the
        *    methodology from the aiaa paper (2006) describing the history and
        *    development of the code.
        *
        *  author        : david vallado                  719-573-2600   28 jun 2005
        *
        *  inputs        :
        *    satrec	 - initialised structure from sgp4init() call.
        *    tsince	 - time since epoch (minutes)
        *
        *  outputs       :
        *    r           - position vector                     km
        *    v           - velocity                            km/sec
        *  return code - non-zero on error.
        *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
        *                   2 - mean motion less than 0.0
        *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
        *                   4 - semi-latus rectum < 0.0
        *                   5 - epoch elements are sub-orbital
        *                   6 - satellite has decayed
        *
        *  locals        :
        *    am          -
        *    axnl, aynl        -
        *    betal       -
        *    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
        *    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
        *    cosisq  , cossu   , sinsu   , cosu    , sinu
        *    delm        -
        *    delomg      -
        *    dndt        -
        *    eccm        -
        *    emsq        -
        *    ecose       -
        *    el2         -
        *    eo1         -
        *    eccp        -
        *    esine       -
        *    argpm       -
        *    argpp       -
        *    omgadf      -
        *    pl          -
        *    r           -
        *    rtemsq      -
        *    rdotl       -
        *    rl          -
        *    rvdot       -
        *    rvdotl      -
        *    su          -
        *    t2  , t3   , t4    , tc
        *    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
        *    u   , ux   , uy    , uz     , vx     , vy     , vz
        *    inclm       - inclination
        *    mm          - mean anomaly
        *    nm          - mean motion
        *    nodem       - right asc of ascending node
        *    xinc        -
        *    xincp       -
        *    xl          -
        *    xlm         -
        *    mp          -
        *    xmdf        -
        *    xmx         -
        *    xmy         -
        *    nodedf      -
        *    xnode       -
        *    nodep       -
        *    np          -
        *
        *  coupling      :
        *    getgravconst- no longer used. Variables are conatined within satrec
        *    dpper
        *    dpspace
        *
        *  references    :
        *    hoots, roehrich, norad spacetrack report #3 1980
        *    hoots, norad spacetrack report #6 1986
        *    hoots, schumacher and glover 2004
        *    vallado, crawford, hujsak, kelso  2006
          ----------------------------------------------------------------------------*/

        public void sgp4
             (
               ref elsetrec satrec, double tsince,
               double[] r, double[] v
             )
        {
            double am, axnl, aynl, betal, cosim, cnod,
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
            // the old check used 1.0 + Math.Cos(Math.PI-1.0e-9), but then compared it to
            // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
            const double temp4 = 1.5e-12;
            twopi = 2.0 * SGPPI;
            x2o3 = 2.0 / 3.0;
            // sgp4fix identify constants and allow alternate values
            // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
            vkmpersec = satrec.radiusearthkm * satrec.xke / 60.0;

            /* --------------------- clear sgp4 error flag ----------------- */
            satrec.t = tsince;
            satrec.error = 0;

            /* ------- update for secular gravity and atmospheric drag ----- */
            xmdf = satrec.mo + satrec.mdot * satrec.t;
            argpdf = satrec.argpo + satrec.argpdot * satrec.t;
            nodedf = satrec.nodeo + satrec.nodedot * satrec.t;
            argpm = argpdf;
            mm = xmdf;
            t2 = satrec.t * satrec.t;
            nodem = nodedf + satrec.nodecf * t2;
            tempa = 1.0 - satrec.cc1 * satrec.t;
            tempe = satrec.bstar * satrec.cc4 * satrec.t;
            templ = satrec.t2cof * t2;

            if (satrec.isimp != 1)
            {
                delomg = satrec.omgcof * satrec.t;
                // sgp4fix use mutliply for speed instead of Math.Pow
                delmtemp = 1.0 + satrec.eta * Math.Cos(xmdf);
                delm = satrec.xmcof *
                         (delmtemp * delmtemp * delmtemp -
                         satrec.delmo);
                temp = delomg + delm;
                mm = xmdf + temp;
                argpm = argpdf - temp;
                t3 = t2 * satrec.t;
                t4 = t3 * satrec.t;
                tempa = tempa - satrec.d2 * t2 - satrec.d3 * t3 -
                                 satrec.d4 * t4;
                tempe = tempe + satrec.bstar * satrec.cc5 * (Math.Sin(mm) -
                                 satrec.sinmao);
                templ = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
                                 satrec.t * satrec.t5cof);
            }

            nm = satrec.no_unkozai;
            em = satrec.ecco;
            inclm = satrec.inclo;
            if (satrec.method == 'd')
            {
                tc = satrec.t;
                dspace
                    (
                      satrec.irez,
                      satrec.d2201, satrec.d2211, satrec.d3210,
                      satrec.d3222, satrec.d4410, satrec.d4422,
                      satrec.d5220, satrec.d5232, satrec.d5421,
                      satrec.d5433, satrec.dedt, satrec.del1,
                      satrec.del2, satrec.del3, satrec.didt,
                      satrec.dmdt, satrec.dnodt, satrec.domdt,
                      satrec.argpo, satrec.argpdot, satrec.t, tc,
                      satrec.gsto, satrec.xfact, satrec.xlamo,
                      satrec.no_unkozai, ref satrec.atime,
                      ref em, ref argpm, ref inclm, ref satrec.xli, ref mm, ref satrec.xni,
                      ref nodem, out dndt, ref nm
                    );
            } // if method = d

            if (nm <= 0.0)
            {
                //         printf("# error nm %f\n", nm);
                satrec.error = 2;
                // sgp4fix add return
                //return false;
            }
            am = Math.Pow((satrec.xke / nm), x2o3) * tempa * tempa;
            nm = satrec.xke / Math.Pow(am, 1.5);
            em = em - tempe;

            // fix tolerance for error recognition
            // sgp4fix am is fixed from the previous nm check
            if ((em >= 1.0) || (em < -0.001)/* || (am < 0.95)*/ )
            {
                //         printf("# error em %f\n", em);
                satrec.error = 1;
                // sgp4fix to return if there is an error in eccentricity
                //return false;
            }
            // sgp4fix fix tolerance to avoid a divide by zero
            if (em < 1.0e-6)
                em = 1.0e-6;
            mm = mm + satrec.no_unkozai * templ;
            xlm = mm + argpm + nodem;
            emsq = em * em;
            temp = 1.0 - emsq;

            nodem = (nodem) % twopi;
            argpm = (argpm) % twopi;
            xlm = (xlm) % twopi;
            mm = (xlm - argpm - nodem) % twopi;

            // sgp4fix recover singly averaged mean elements
            satrec.am = am;
            satrec.em = em;
            satrec.im = inclm;
            satrec.Om = nodem;
            satrec.om = argpm;
            satrec.mm = mm;
            satrec.nm = nm;

            /* ----------------- compute extra mean quantities ------------- */
            sinim = Math.Sin(inclm);
            cosim = Math.Cos(inclm);

            /* -------------------- add lunar-solar periodics -------------- */
            ep = em;
            xincp = inclm;
            argpp = argpm;
            nodep = nodem;
            mp = mm;
            sinip = sinim;
            cosip = cosim;
            if (satrec.method == 'd')
            {
                dpper
                    (
                      satrec.e3, satrec.ee2, satrec.peo,
                      satrec.pgho, satrec.pho, satrec.pinco,
                      satrec.plo, satrec.se2, satrec.se3,
                      satrec.sgh2, satrec.sgh3, satrec.sgh4,
                      satrec.sh2, satrec.sh3, satrec.si2,
                      satrec.si3, satrec.sl2, satrec.sl3,
                      satrec.sl4, satrec.t, satrec.xgh2,
                      satrec.xgh3, satrec.xgh4, satrec.xh2,
                      satrec.xh3, satrec.xi2, satrec.xi3,
                      satrec.xl2, satrec.xl3, satrec.xl4,
                      satrec.zmol, satrec.zmos, satrec.inclo,
                      'n', ref ep, ref xincp, ref nodep, ref argpp, ref mp, satrec.operationmode
                    );
                if (xincp < 0.0)
                {
                    xincp = -xincp;
                    nodep = nodep + SGPPI;
                    argpp = argpp - SGPPI;
                }
                if ((ep < 0.0) || (ep > 1.0))
                {
                    //            printf("# error ep %f\n", ep);
                    satrec.error = 3;
                    // sgp4fix add return
                    //return false;
                }
            } // if method = d

            /* -------------------- long period periodics ------------------ */
            if (satrec.method == 'd')
            {
                sinip = Math.Sin(xincp);
                cosip = Math.Cos(xincp);
                satrec.aycof = -0.5 * satrec.j3oj2 * sinip;
                // sgp4fix for divide by zero for xincp = 180 deg
                if (Math.Abs(cosip + 1.0) > 1.5e-12)
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
                else
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
            }
            axnl = ep * Math.Cos(argpp);
            temp = 1.0 / (am * (1.0 - ep * ep));
            aynl = ep * Math.Sin(argpp) + temp * satrec.aycof;
            xl = mp + argpp + nodep + temp * satrec.xlcof * axnl;

            /* --------------------- solve kepler's equation --------------- */
            u = (xl - nodep) % twopi;
            eo1 = u;
            tem5 = 9999.9;
            ktr = 1;
            // sgp4fix for c# intiialize
            coseo1 = 0.0;
            sineo1 = 0.0;
            //   sgp4fix for kepler iteration
            //   the following iteration needs better limits on corrections
            while ((Math.Abs(tem5) >= 1.0e-12) && (ktr <= 10))
            {
                sineo1 = Math.Sin(eo1);
                coseo1 = Math.Cos(eo1);
                tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
                tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
                if (Math.Abs(tem5) >= 0.95)
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
                //         printf("# error pl %f\n", pl);
                satrec.error = 4;
                // sgp4fix add return
                //return false;
            }
            else
            {
                rl = am * (1.0 - ecose);
                rdotl = Math.Sqrt(am) * esine / rl;
                rvdotl = Math.Sqrt(pl) / rl;
                betal = Math.Sqrt(1.0 - el2);
                temp = esine / (1.0 + betal);
                sinu = am / rl * (sineo1 - aynl - axnl * temp);
                cosu = am / rl * (coseo1 - axnl + aynl * temp);
                su = Math.Atan2(sinu, cosu);
                sin2u = (cosu + cosu) * sinu;
                cos2u = 1.0 - 2.0 * sinu * sinu;
                temp = 1.0 / pl;
                temp1 = 0.5 * satrec.j2 * temp;
                temp2 = temp1 * temp;

                /* -------------- update for short period periodics ------------ */
                if (satrec.method == 'd')
                {
                    cosisq = cosip * cosip;
                    satrec.con41 = 3.0 * cosisq - 1.0;
                    satrec.x1mth2 = 1.0 - cosisq;
                    satrec.x7thm1 = 7.0 * cosisq - 1.0;
                }
                mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) +
                        0.5 * temp1 * satrec.x1mth2 * cos2u;
                su = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
                xnode = nodep + 1.5 * temp2 * cosip * sin2u;
                xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
                mvt = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke;
                rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
                        1.5 * satrec.con41) / satrec.xke;

                /* --------------------- orientation vectors ------------------- */
                sinsu = Math.Sin(su);
                cossu = Math.Cos(su);
                snod = Math.Sin(xnode);
                cnod = Math.Cos(xnode);
                sini = Math.Sin(xinc);
                cosi = Math.Cos(xinc);
                xmx = -snod * cosi;
                xmy = cnod * cosi;
                ux = xmx * sinsu + cnod * cossu;
                uy = xmy * sinsu + snod * cossu;
                uz = sini * sinsu;
                vx = xmx * cossu - cnod * sinsu;
                vy = xmy * cossu - snod * sinsu;
                vz = sini * cossu;

                /* --------- position and velocity (in km and km/sec) ---------- */
                r[0] = (mrt * ux) * satrec.radiusearthkm;
                r[1] = (mrt * uy) * satrec.radiusearthkm;
                r[2] = (mrt * uz) * satrec.radiusearthkm;
                v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
                v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
                v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
            }  // if pl > 0

            // sgp4fix for decaying satellites
            if (mrt < 1.0)
            {
                //         printf("# decay condition %11.6f \n",mrt);
                satrec.error = 6;
                //return false;
            }

            //#include "debug7.cpp"
            //return true;
        }  // end sgp4




#endif

