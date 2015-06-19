#include "slalib.h"
#include "slamac.h"
double slaDt ( double epoch )
/*
**  - - - - - -
**   s l a D t
**  - - - - - -
**
**  Estimate the offset between dynamical time and Universal Time
**  for a given historical epoch.
**
**  (double precision)
**
**  Given:
**     epoch    double    (Julian) epoch (e.g. 1850.0)
**
**  The result is a rough estimate of ET-UT (after 1984, TT-UT1) at
**  the given epoch, in seconds.
**
**  Notes:
**
**  1  Depending on the epoch, one of three parabolic approximations
**     is used:
**
**      before 979    Stephenson & Morrison's 390 BC to AD 948 model
**      979 to 1708   Stephenson & Morrison's 948 to 1600 model
**      after 1708    McCarthy & Babcock's post-1650 model
**
**     The breakpoints are chosen to ensure continuity:  they occur
**     at places where the adjacent models give the same answer as
**     each other.
**
**  2  The accuracy is modest, with errors of up to 20 sec during
**     the interval since 1650, rising to perhaps 30 min by 1000 BC.
**     Comparatively accurate values from AD 1600 are tabulated in
**     the Astronomical Almanac (see section K8 of the 1995 AA).
**
**  3  The use of double-precision for both argument and result is
**     purely for compatibility with other SLALIB time functions.
**
**  4  The models used are based on a lunar tidal acceleration value
**     of -26.00 arcsec per century.
**
**  Reference:  Explanatory Supplement to the Astronomical Almanac,
**              ed P.K.Seidelmann, University Science Books (1992),
**              section 2.553, p83.  This contains references to
**              the Stephenson & Morrison and McCarthy & Babcock
**              papers.
**
**  Last revision:   22 October 2006
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double t, w, s;


/* Centuries since 1800 */
   t = ( epoch - 1800.0 ) / 100.0;

/* Select model */
   if ( epoch >= 1708.185161980887 ) {

   /* Post-1708: use McCarthy & Babcock */
      w = t - 0.19;
      s = 5.156 + 13.3066 * w * w;
   } else {
      if ( epoch >= 979.0258204760233 ) {

   /* 979-1708: use Stephenson & Morrison's 948-1600 model */
         s = 25.5 * t * t;
      } else {

      /* Pre-979: use Stephenson & Morrison's 390 BC to AD 948 model */
         s = 1360.0 + ( 320.0 + 44.3 * t ) * t;
      }
   }

/* Result */
   return s;
}
