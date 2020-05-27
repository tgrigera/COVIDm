/*
 * esampler.hh -- classes to obtain samples and statistic form graph
 *                epidemiological models
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

#ifndef ESAMPLER_HH
#define ESAMPLER_HH

#include <iostream>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
//
// Sampler
//
// This works as Gillespie_sampler, but is designed with a pure
// virtual push_data metod.  It iis a simple template class which acts
// as a filter between the data generation and actual data output or
// averaging.  It's function is to send to output the state of the
// system at regular intervals, adequately sampling data generated
// with varying time step by the Gillespie algorithm.
//
// The template parameter is a boolean, which controls whether the
// class actually performs the sampling.  If true (default), it
// semples as described above, if false it simply passes everything
// through to the pusher.

class Collector {
public:
  virtual void collect(double time)=0;
} ;

class Sampler {
public:
  Sampler(Collector *collector);
  virtual void sample(double time)=0;

protected:
  Collector* collector;
} ;
  
inline Sampler::Sampler(Collector *collector) :
  collector(collector)
{}

class Passthrough_sampler : public Sampler {
public:
  Passthrough_sampler(double tmax,Collector *collector);
  void sample(double time);

private:
  double tmax;

} ;

inline Passthrough_sampler::Passthrough_sampler(double tmax,Collector *collector) :
  Sampler(collector), tmax(tmax)
{}

inline void Passthrough_sampler::sample(double time)
{
  if (time<=tmax)
    collector->collect(time);
}

///////////////////////////////////////////////////////////////////////////////
//
// Gillespie_sampler

class Gillespie_sampler : public Sampler {
public:
  Gillespie_sampler(double t0,double tmax,double deltat,Collector *collector);
  void sample(double time);

private:
  double     t0,deltat,tmax;
  double     tlast,tnext;
} ;
  
inline Gillespie_sampler::Gillespie_sampler(double t0,double tmax,double deltat,Collector *collector) :
  Sampler(collector),
  t0(t0),
  deltat(deltat),
  tmax(tmax),
  tlast(t0-deltat),
  tnext(t0)
{}
  

inline void Gillespie_sampler::sample(double time)
{
  if (time<=tnext || tlast>=tmax) return;
  double t;
  for (t=tlast+=deltat; t<time && t<=tmax; t+=deltat) 
    collector->collect(t);
  tnext=t;
  tlast=t-deltat;
}


#endif /* ESAMPLER_HH */
