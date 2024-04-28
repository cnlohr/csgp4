# sattrack

C, header-only port of David Vallado's SGP4Lib, for use in embedded, and unusual situations.

Approximate values per-core on a AMD 5800X:
Init: 0.469150 us/iteration
Run:  0.226054 us/iteration

This is kind of rough, so USE AT YOUR OWN RISK.

All I've done to validate is that the orbital tracks betwen this and the Python SGP4 lib come up similiarly.

## Resources
 * SGP4 Transliterated from https://celestrak.org/software/vallado-sw.php, specifically https://celestrak.org/software/vallado/cs.zip

 * https://en.wikipedia.org/wiki/Two-line_element_set

 * https://spaceaware.io/
 * https://cdn.digitalarsenal.io/celestrak/NORAD/elements/catalog.txt


 * https://celestrak.org/NORAD/elements/
 * Space Stations: https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle


