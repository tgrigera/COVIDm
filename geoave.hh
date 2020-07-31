/*
 * geoave.hh -- Class to average in geometrically growing windows
 *
 * Copied from the glsim project
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef GEOAVE_HH
#define GEOAVE_HH

#include <ostream>
#include <vector>
#include <stdexcept>
#include <errno.h>
#include <math.h>

/** \class Geoave
    \ingroup Analysis
    \brief Average time series in geometrically growing windows

    This class takes pairs of numbers \f$x,y\f$ (through push()) and
    averages together all points whose x-coordinate (assume
    it is a time) falls within a window.  The first window starts
    at \f$t_0\f$ and is of length \f$b\f$, successive windows grow
    geometrically by factor \f$w\f$.  In other words, window
    boundaries are located at

    \f[ t_n = t_0 + b \frac{w^n - 1}{w-1} \f]

    When returnting the averages (get_aves()) or when printing to a
    stream (operator<<), the  time assigned to the
    average value is the center of the window, i.e. the
    abscissas \f$s_n\f$ are

    \f[  s_n = t_n + \frac{1}{2} b w^n \f]

    Points may be given in any order, and averages are available at
    any time.  Points can be added after querying for the accumulated
    averages.

    The case \f$w=1\f$ is supported and handled as a special case,
    yielding windows of fixed width equal to \f$b\f$. \f$w<1\f$ is not
    recommended.
*/
class Geoave {
public:
  /** Constructor.
      \param  t0  Time origin (\f$t_0\f$)
      \param  wfactor Window factor (\f$w\f$)
      \param  base  Size of first window (\f$b\f$)
  */
  Geoave(double t0=0,double wfactor=1.5,double base=1);

  void push(double time,double e); ///< Add pair (first argument is interpreted as time)
  void get_aves(std::vector<double>& time,std::vector<double>& ave,
		std::vector<double>& var) const; ///< Returns the abcissas \f$s_0\f$, averages and variance (vectors given are first cleared)
  double ave(int i) const {return rave[i];}
  double var(int i) const {return rvarn[i]/(count[i]-1); }
  int Nsamp(int i) const {return count[i];}

private:
  double              base,t0,wfactor;
  double              logwf,read_fb;

  std::vector<unsigned long>   count;
  std::vector<double> rave,rvarn;

  /// Write times, averages and average error (sqrt(variance/N)) to stream (three columns)
  friend std::ostream& operator<<(std::ostream&,const Geoave&); 
} ;

#endif /* GEOAVE_HH */

