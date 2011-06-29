#include "slalib.h"
#include "slamac.h"
void slaMap ( double rm, double dm, double pr, double pd,
              double px, double rv, double eq, double date,
              double *ra, double *da )
/*
**  - - - - - - -
**   s l a M a p
**  - - - - - - -
**
**  Transform star RA,Dec from mean place to geocentric apparent.
**
**  The reference frames and timescales used are post IAU 1976.
**
**  References:
**     1984 Astronomical Almanac, pp B39-B41.
**     (also Lederle & Schwan, Astron. Astrophys. 134, 1-6, 1984)
**
**  Given:
**     rm,dm    double     mean RA,Dec (rad)
**     pr,pd    double     proper motions:  RA,Dec changes per Julian year
**     px       double     parallax (arcsec)
**     rv       double     radial velocity (km/sec, +ve if receding)
**     eq       double     epoch and equinox of star data (Julian)
**     date     double     TDB for apparent place (JD-2400000.5)
**
**  Returned:
**     *ra,*da  double     apparent RA,Dec (rad)
**
**  Called:
**     slaMappa       star-independent parameters
**     slaMapqk       quick mean to apparent
**
**  Notes:
**
**  1)  eq is the Julian epoch specifying both the reference frame and
**      the epoch of the position - usually 2000.  For positions where
**      the epoch and equinox are different, use the routine slaPm to
**      apply proper motion corrections before using this routine.
**
**  2)  The distinction between the required TDB and TT is always
**      negligible.  Moreover, for all but the most critical
**      applications UTC is adequate.
**
**  3)  The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt.
**
**  4)  This routine may be wasteful for some applications because it
**      recomputes the Earth position/velocity and the precession-
**      nutation matrix each time, and because it allows for parallax
**      and proper motion.  Where multiple transformations are to be
**      carried out for one epoch, a faster method is to call the
**      slaMappa routine once and then either the slaMapqk routine
**      (which includes parallax and proper motion) or slaMapqkz (which
**      assumes zero parallax and proper motion).
**
**  5)  The accuracy is limited by imperfections in the IAU 1976/1980
**      models for precession and nutation.  Corrections are tabulated
**      in IERS Bulletin B and at the present epoch are of order 50 mas.
**      An improved precession-nutation model can be introduced by
**      using slaMappa and slaMapqk (see the previous note) and
**      replacing the precession-nutation matrix into the parameter
**      array directly.
**
**  6)  The accuracy is further limited by the routine slaEvp, called
**      by slaMappa, which computes the Earth position and velocity
**      using the methods of Stumpff.  The maximum error is about
**      0.3 mas.
**
**  Last revision:   8 May 2000
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double amprms[21];

/* Star-independent parameters */
   slaMappa ( eq, date, amprms );

/* Mean to apparent */
   slaMapqk ( rm, dm, pr, pd, px, rv, amprms, ra, da );
}
