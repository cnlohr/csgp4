# csgp4

<IMG SRC=https://github.com/cnlohr/csgp4/actions/workflows/build.yml/badge.svg>

C, header-only port of David Vallado's SGP4Lib, for use in embedded, and unusual situations.

Approximate values per-core on a AMD 5800X:

| Operation | Run Time |
| --- | --- |
| Init | 0.3993 us/iteration |
| Run |  0.2270 us/iteration |
| Init and Run (ISS) | 0.4059 us/iteration |
| Init and Run (Beyond Geostationary) | 1.0818 us/iteration |

This is kind of rough, so USE AT YOUR OWN RISK. All I've done to validate is that the orbital tracks betwen this and the Python SGP4 lib come up similiarly.

## How to use

Check both [checkProg.c](https://github.com/cnlohr/csgp4/blob/master/checkProg.c) and [checkProgSimple.c](https://github.com/cnlohr/csgp4/blob/master/checkProgSimple.c) for general usage. More detailed information follows:

### How to use csgp4.h

Just `#include csgp4.h` at the top of your program, from there... 

The idea is that you can either read directly from a TLE file, which is acquirable from celestrak, you can see general usage here https://github.com/cnlohr/csgp4/blob/master/checkProg.c#L8-L29 - but in general, you open the file, and pass it to `ParseFileOrString`.  With the same function, you can also also pass in a string containing the TLEs you want to load instead.

Once the TLEs are loaded, you can use `sgp4` to compute the position of your desired satellite at the desired time.

Output is in the TEME (Earth Intertial Frame), and input times are the relative time difference between the epoch and the time in seconds.

### How to use csgp4_simple.h

csgp4_simple.h is a lot like csgp4.h, however, it assumes you have already parsed the TLE records, and it also just computes the object's position at a time and does not initialize.  It is geared for situations where you can already pull an object's TLE.  For instance, on a microcontroller pulling data from a server that already runs the full csgp4 processing, or on the GPU.

## Floating point performance

We also test single precision floating point performance, and in general, over the course of a day for most satellites, it is off by a few meters, or up to 5 or 6 km over the course of a month, as compared to double precision.

Care should be taken surroinding the date format sent in.  In many places in code, we split times into "days" and "fractional days", i.e. `jd` and `jdfrac` (Julian days being days since 4713 bc).  Only once you have subtracted your epoch and this time should you squish the values together into one float.  This is because when floats become very large, their precision for handling small numbers gets worse and worse.  For instance, at 8 million, you can discern 8,000,000 and 8,000,000.5 but, at 17,000,000, 17,000,000.5 just resolves to 17,000,000.

## TODO
 * Figure out why Deep Space is off some.
 * Can we avoid extra call to sgp4 in sgp4init?
 * Is meanMotion1, meanMotion2 used at all?  Can we just remove them?  (and ndot)
 * Make it so `sgp4init_simple` can take epoch as epocf and epoch fraction.


## Resources

 * SGP4 Transliterated from https://celestrak.org/software/vallado-sw.php, specifically https://celestrak.org/software/vallado/cs.zip
 * https://en.wikipedia.org/wiki/Two-line_element_set
 * https://spaceaware.io/
 * https://cdn.digitalarsenal.io/celestrak/NORAD/elements/catalog.txt
 * https://celestrak.org/NORAD/elements/
 * Space Stations: https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle
 *  https://en.wikipedia.org/wiki/Julian_day

