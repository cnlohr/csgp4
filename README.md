# csgp4

C, header-only port of David Vallado's SGP4Lib, for use in embedded, and unusual situations.

Approximate values per-core on a AMD 5800X:

| Operation | Run Time |
| --- | --- |
| Init | 0.4692 us/iteration |
| Run |  0.2261 us/iteration |

This is kind of rough, so USE AT YOUR OWN RISK. All I've done to validate is that the orbital tracks betwen this and the Python SGP4 lib come up similiarly.

## TODO
 * Figure out why Deep Space is off some.
 * Can we avoid extra call to sgp4 in sgp4init?
 * Is meanMotion1, meanMotion2 used at all?  Can we just remove them?  (and ndot)

## Resources

 * SGP4 Transliterated from https://celestrak.org/software/vallado-sw.php, specifically https://celestrak.org/software/vallado/cs.zip
 * https://en.wikipedia.org/wiki/Two-line_element_set
 * https://spaceaware.io/
 * https://cdn.digitalarsenal.io/celestrak/NORAD/elements/catalog.txt
 * https://celestrak.org/NORAD/elements/
 * Space Stations: https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle


