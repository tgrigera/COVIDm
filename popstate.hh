/*
 * popstate.hh -- structures to record, print and average the
 *                epidemiological state of the whole population.
 *
 * This file is part of COVIDm.
 *
 * COVIDm is copyright (C) 2020 by the authors (see file AUTHORS)
 * 
 * COVIDm is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * COVIDm is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE.
 *
 */

#ifndef POPSTATE_HH
#define POPSTATE_HH

#include <ostream>
#include <vector>

#include "geoave.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Population_state

class Population_state {
public:
  Population_state(int NS=1,int NE=1,int NI=1,int NR=1);
  virtual ~Population_state() {}

  virtual const char *header();
  virtual void print(std::ostream&,bool print_time=true);

  double time;
  std::vector<double> S,E,I,R;

protected:
  void addSIRhdr();
  std::string hdr;
} ;
	       
inline Population_state::Population_state(int NS,int NE,int NI,int NR) :
  S(NS,0.), E(NE,0.), I(NI,0.), R(NR,0.)
{}

inline std::ostream& operator<<(std::ostream& o,Population_state& state)
{
  state.print(o);
  return o;
}

///////////////////////////////////////////////////////////////////////////////
//
// SIRstate

class SIRstate : public Population_state {
public:
  SIRstate() :
    Population_state(1,0,1,1) {}
  virtual void push(double time,double S,double I,double R);
} ;

///////////////////////////////////////////////////////////////////////////////
//
// SIRstate_av

class SIRstate_av : public SIRstate {
public:
  SIRstate_av(double deltat=1.) :
    SIRstate(),
    Sav(-0.5*deltat,1.,deltat),
    Iav(-0.5*deltat,1.,deltat),
    Rav(-0.5*deltat,1.,deltat)
  {}

  const char* header();
  void print(std::ostream&,bool print_time=true);
  void push(double time,double S,double I,double R);

private:
  Geoave Sav,Iav,Rav;

} ;

#endif /* POPSTATE_HH */
