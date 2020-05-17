/*
 * avevar.hh -- Class for average and variance with West recurrence
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

#ifndef AVEVAR_HH
#define AVEVAR_HH

#include <limits>
#include <vector>

/** \class AveVar
    \ingroup Analysis
    \brief Compute average and variance using West recurrence
*/
template <bool maxmin=true>
class AveVar {
public:
  AveVar() {clear();};
  AveVar& push(double);
  AveVar& push(std::vector<double> v) {for (auto d : v) push(d);}
  AveVar& clear();

  double  ave() const {return ave_;}
  double  var() const {return var_/(N_-1);}
  long    N() const {return N_;}
  double  max() const {return max_;}
  double  min() const {return min_;}
  
private:
  double ave_,var_;
  long   N_;
  double max_,min_;

} ;

template <bool maxmin>
inline AveVar<maxmin>& AveVar<maxmin>::clear()
{
  ave_=0;
  var_=0;
  N_=0;
  max_=std::numeric_limits<double>::min();
  min_=std::numeric_limits<double>::max();
  return *this;
}
  

template <bool maxmin>
inline AveVar<maxmin>& AveVar<maxmin>::push(double x)
{
  double Q,R;

  N_++;
  Q=x-ave_;
  R=Q/N_;
  ave_+=R;
  var_+=Q*R*(N_-1);

  if (maxmin) {
    if (x<min_) min_=x;
    if (x>max_) max_=x;
  }

  return *this;
}

#endif /* AVEVAR_HH */
