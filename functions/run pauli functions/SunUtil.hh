
#ifndef SRC_UTIL_SUNUTIL_HH_
#define SRC_UTIL_SUNUTIL_HH_

#include <TVector3.h>

namespace RAT {
  /*
   * Returns the direction from the detector to the Sun for a given date.
   * Original implementation by Phil Harvey.
   *
   * This method implements the calculation of the direction from the
   * detector to the Sun. The zenith angle is not enough to generate the
   * direction from the Sun and therefore this function implementation by
   * Phil Harvey is used. Phil's original notes:
   *
   * Calculates the direction to the sun in detector coordinates for the given time in UTC with zero time of Jan 0,1975.
   * (ref: http://hotel04.ausys.se/pausch/comp/ppcomp.html)
   *  Verification: This method agrees with US Navy data to within 0.6 degrees between 1990 and 2020,
   *  and agrees with the SNOMAN calculation to within 0.03 degrees for the year 2000.
   *  (US Navy data from http://aa.usno.navy.mil/AA/data/docs/AltAz.html)
   *
   * This method has then been updated to use the zero time definition of SNO+ (00:00:00 Jan 1 2010).
   *
   * http://hotel04.ausys.se/pausch/comp/ppcomp.html, http://aa.usno.navy.mil/AA/data/docs/AltAz.html
   *
   * UTdays   Number of days since SNOMAN  time zero (31 December 1974).
   * UTsecs   Number of seconds since the start of the present day.
   * UTnsecs  Number of nanoseconds since start of present second.
   * labtwist Rotation angle to align the X axis with the detector construction North.
   * Returns: Direction from the detector to the Sun.
   * Attention: The direction is from the detector TO the Sun. For use in the generator it has to be inverted.
   */

  TVector3 SunDirection(const int UTdays, const int UTsecs, const int UTnsecs, const double labtwist=-49.58);
}
#endif /* SRC_UTIL_SUNUTIL_HH_ */
