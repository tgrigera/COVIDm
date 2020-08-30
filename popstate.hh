/*
 * popstate.hh -- structures to record, print and average the
 *                (global) epidemiological state of the whole population.
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

//
// Print source context (for debugging)
#define I_AM_HERE \
  std::cerr << "(DD) Reached " << __FILE__ << ":" << __LINE__ << " (in function " \
  <<	 __func__ << ")\n";

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
  void create_colnums(int);
  void addSIRhdr();
  std::string hdr;
  std::string colnums;
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

struct SIRistate {
  double S,I,R;
} ;

class SIRstate : public Population_state {
public:
  SIRstate() :
    Population_state(1,0,1,1) {}
  virtual void push(double time,SIRistate &istate);
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
  void push(double time,SIRistate &istate);

private:
  Geoave Sav,Iav,Rav;

} ;

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRstate

struct SEEIIRistate {
  int    N;                  // Total population
  int    S,E1,E2,I1,I2,R;    // Total in states SEIR
  // cumulative infections by type
  int    inf_imported,inf_close,inf_community;
  int    Eacc;
  double beta_out;
  double tinf;
} ;


class SEEIIRstate : public Population_state {
public:
  SEEIIRstate() :
    Population_state(1,2,2,1),
    time0(0), Eacc0(0), I0(0)
  {time=-10;}
  const char* header();
  virtual void push(double time,SEEIIRistate &s);

private:
  double time0;
  int    Eacc0,I0;
} ;


///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRstate_av

class SEEIIRstate_av : public SEEIIRstate {
public:
  SEEIIRstate_av(double deltat=1.) :
    SEEIIRstate(),
    Sav(-0.5*deltat,1.,deltat),
    E1av(-0.5*deltat,1.,deltat),
    E2av(-0.5*deltat,1.,deltat),
    I1av(-0.5*deltat,1.,deltat),
    I2av(-0.5*deltat,1.,deltat),
    Rav(-0.5*deltat,1.,deltat),
    Impav(-0.5*deltat,1.,deltat),
    Closeav(-0.5*deltat,1.,deltat),
    Commav(-0.5*deltat,1.,deltat),
    Nav(-0.5*deltat,1.,deltat),
    betaav(-0.5*deltat,1.,deltat),
    RRav(-0.5*deltat,1.,deltat),
    time0(0), Eacc0(0), I0(0)
  {}

  const char* header();
  void print(std::ostream&,bool print_time=true);
  void push(double time,SEEIIRistate &s);

private:
  Geoave Sav,E1av,E2av,I1av,I2av,Rav,Impav,Closeav,Commav,Nav,betaav,RRav;
private:
  double time0;
  int    Eacc0,I0;
} ;


#endif /* POPSTATE_HH */
