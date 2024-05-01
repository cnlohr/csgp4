#include <stdio.h>

#define CSGP4_USE_FLOAT 0
#define CSGP4_INIT 0

#include "csgp4.h"

int main( int argc, char ** argv )
{
	double ro[3];
	double vo[3];
	struct elsetrec iss;
	double startmfe = (1714240800 - /* iss.epoch */ 1000000)/60.0;
	sgp4 (&iss, startmfe, ro,  vo);
	return 0;
}
