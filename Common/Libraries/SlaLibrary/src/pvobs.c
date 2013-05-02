#include "slalib.h"
#include "slamac.h"
void slaPvobs ( double p, double h, double stl, double pv[6] )
/*
**  - - - - - - - - -
**   s l a P v o b s
**  - - - - - - - - -
**
**  Position and velocity of an observing station.
**
**  (double precision)
**
**  Given:
**     p     double     latitude (geodetic, radians)
**     h     double     height above reference spheroid (geodetic, metres)
**     stl   double     local apparent sidereal time (radians)
**
**  Returned:
**     pv    double[6]  position/velocity 6-vector (au, au/s, true
**                                         equator and equinox of date)
**
**  Notes:
**
**   1) IAU 1976 constants are used.
**
**   2) Polar motion is neglected.
**
**   3) If local Earth rotation angle is supplied instead of local
**      apparent sidereal time, the result is with respect to the
**      celestial intermediate origin instead of the equinox.
**
**  Called:  slaGeoc
**
**  Last revision:   18 March 2013
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

#define SR 7.292115855306589e-5  /* Mean sidereal rate (at J2000)
                                    in radians per (UT1) second */

{
   double r, z, s, c, v;

/* Geodetic to geocentric conversion */
   slaGeoc ( p, h, &r, &z );

/* Functions of ST */
   s = sin ( stl );
   c = cos ( stl );

/* Speed */
   v = SR * r;

/* Position */
   pv[0] = r * c;
   pv[1] = r * s;
   pv[2] = z;

/* Velocity */
   pv[3] = - v * s;
   pv[4] = v * c;
   pv[5] = 0.0;
}
