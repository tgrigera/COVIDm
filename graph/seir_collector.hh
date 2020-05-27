/*
 * seir_collector.hh --  collector classes for SEIR models
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

#ifndef SEIR_COLLECTOR_HH
#define SEIR_COLLECTOR_HH

#include "esampler.hh"
#include "sirmodel.hh"

///////////////////////////////////////////////////////////////////////////////
//
// SEIRcollector_base

class SEIRcollector_base : public Collector {
public:
  SEIRcollector_base(int NS=1,int NE=1,int NI=1,int NR=1);

  virtual const char *header();
  virtual void print(std::ostream&,bool print_time=true);
  virtual ~SEIRcollector_base() {}

  double time;
  std::vector<double> S,E,I,R;

protected:
  void create_colnums(int);
  void addSIRhdr();
  std::string hdr;
  std::string colnums;
} ;
	       
inline SEIRcollector_base::SEIRcollector_base(int NS,int NE,int NI,int NR) :
  S(NS,0.), E(NE,0.), I(NI,0.), R(NR,0.)
{}

inline std::ostream& operator<<(std::ostream& o,SEIRcollector_base& state)
{
  state.print(o);
  return o;
}

///////////////////////////////////////////////////////////////////////////////
//
// SIRcollector

template <typename IGraph>
class SIRcollector : public SEIRcollector_base {
public:
  SIRcollector(SIR_model<IGraph> &model) :
    SEIRcollector_base(1,0,1,1),
    model(model)
  {}
  void collect(double time);

protected:
  SIR_model<IGraph> &model;
} ;

template <typename IGraph>
void SIRcollector<IGraph>::collect(double time_)
{
  time=time_;
  S[0]=model.NS;
  I[0]=model.NI;
  R[0]=model.NR;

  std::cout << *this << '\n';
}

// class SIRcollector_av : public SIRcollector {
// public:
//   SIRstate_av(double deltat=1.) :
//     SIRstate(),
//     Sav(-0.5*deltat,1.,deltat),
//     Iav(-0.5*deltat,1.,deltat),
//     Rav(-0.5*deltat,1.,deltat)
//   {}

//   const char* header();
//   void print(std::ostream&,bool print_time=true);
//   void push(double time,SIRistate &istate);

// private:
//   Geoave Sav,Iav,Rav;

// } ;



#endif /* SEIR_COLLECTOR_HH */
