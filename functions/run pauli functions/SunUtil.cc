/*
 * Miscellaneous utilities relating to solar neutrinos and to the Sun
 * in general. The utilities use ROOT declarations (instead of G4) to
 * facilitate its use from ROOT, where some of these utilities might be
 * useful.
 *
 * Author: N. Barros <nfbarros@hep.upenn.edu>
 *
 *  Created on: Feb 10, 2015
 *      Author: nbarros
 */

#include <RAT/SunUtil.hh>
#include <RAT/Log.hh>
#include <cmath>

#include <CLHEP/Units/PhysicalConstants.h>
using CLHEP::pi;
using CLHEP::twopi;

#include <TRotation.h>


namespace RAT {

TVector3 SunDirection(const int UTdays, const int UTsecs, const int UTnsecs, const double labtwist) {
    // - N. Barros 13/01/2012
    //
    // Applied a correction for the new RAT time zero definition, which is 00:00 01 Jan 2010
    // Internally the calculations still use a time zero definition of Jan 0 2000, which
    // is the time zero for the astronomical constants.

    // Calculates the direction to the sun in detector coordinates
    // for the given time in UTC with zero time of Jan 0,1975.
    // (ref: http://hotel04.ausys.se/pausch/comp/ppcomp.html)
    //
    // Verification: This method agrees with US Navy data to within
    // 0.6 degrees between 1990 and 2020, and agrees with the SNOMAN
    // calculation to within 0.03 degrees for the year 2000.
    // (US Navy data from http://aa.usno.navy.mil/AA/data/docs/AltAz.html)
    //
    // - P. Harvey 06/16/00

    // constant to convert to radians
    const double to_rad = pi / 180.0;

    //2442412.500000 -> Jan 0 1975
    //2451543.500000 -> Jan 0 2000
    //2455197.500000 -> Jan 1 2010
    // offset in days from SNOMAN time zero (Jan 0,1975) to time zero
    // for the astronomical constants (Jan 0,2000)
    //  const double days_offset = 9131;

    // offset between SNO+ time zero and time zero for astronomical constants.
    // This calculation leads to a negative number. The reason is to keep the calculations
    // below untouched. The offset between 2010 and 2000 has a different sign than the offset between 1975 and 2000
    // as the calculation was originally using.
    const double days_offset =  2451543.5 /*31 Dec 1999*/- 2455197.5 /*Jan 1 2010*/;

    // Don't touch anything below this point.

    // the number of days since Jan 0, 2000 (time zero for astronomical constants)
    double days = UTdays + UTsecs/86400.0 + UTnsecs/(1.0e9*86400.0) - days_offset;
    // the obliquity of the ecliptic
    double ecl = (23.4393 - 3.563E-7 * days) * to_rad;
    // the argument of perihelion
    double w = (282.9404 + 4.70935E-5 * days) * to_rad;
    // the eccentricity of earth's orbit
    double e = 0.016709 - 1.151E-9 * days;
    // the mean anomaly
    double ma = (356.0470 + 0.9856002585 * days) * to_rad;
    // the eccentric anomaly
    double ea = ma + e * sin(ma) * (1.0 + e * cos(ma));

    // the direction to the sun relative to the major axis of the elliptical orbit
    double xv = cos(ea) - e;
    double yv = sqrt(1.0 - e*e) * sin(ea);

    // the true anomaly (angle between position and perihelion)
    double v = atan2(yv, xv);

    // the Sun's true longitude in ecliptic rectangular geocentric coordinates
    double sun_longitude = v + w;

    // calculate the direction to sun including earth tilt
    // (x is east, y is north, z is away from the earth)
    double xs = cos(sun_longitude);
    double ys = sin(sun_longitude);
    TVector3 vec1(ys*cos(ecl), ys*sin(ecl), xs);

    // calculate rotation of the earth about it's axis for this time of day
    // (k0 is empirically generated to agree with US Navy data - PH 06/16/00)
    double k0 = 0.27499;               // rotational position at time zero
    double k1 = 1.0 + 1.0 / 365.2425;  // rotational period of the earth in days
    double spin = k0 + k1 * days;      // spin of the earth (number of revolutions)
    // compute the spin angle (reducing to the range -2*PI to 2*PI)
    double spin_angle = (spin - (long)spin) * twopi;

    double longitude = (81+12./60+5./3600) * to_rad;
    double latitude  = (46+28./60+31.0/3600) * to_rad;

    TRotation rotator;
    rotator.RotateY(longitude - spin_angle);
    rotator.RotateX(latitude);
    rotator.RotateZ(labtwist * to_rad);

    TVector3 sun_dir = rotator*vec1;

#ifdef RATDEBUG
    detail << "[SunUtil]::SunDirection :  Setting direction to Sun as [ " << sun_dir.X() << " , " << sun_dir.Y() << " , " << sun_dir.Z() << " ]."<< newline;
#endif

    return sun_dir;

}

}
