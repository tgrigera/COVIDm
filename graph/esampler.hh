/*
 * esampler.hh -- classes to obtain samples and statistic form kinetic Monte Carlo
 *                (Gillespie) simulaitons
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

/*
 * This file defines classes to take appropriate samples from the time
 * series generated by a kinetic Monte Carlo (also called Gillespie)
 * simulation. The idea is that they act as a filter between the data
 * generation and actual data output or averaging, sending to output
 * data at regular intervals, adequately sampling data generated with
 * varying time step by the Gillespie algorithm.
 *
 * These classes only decide when to take the samples; actual data
 * collection is done by a class derived from Collector.  In this way
 * we abstract the picking of the right time to sample from the actual
 * data collection: the classes in this file have no knowledge of the
 * kind of data they are supposed to collect, that information is
 * managed by the Collector descendent.
 *
 * To use these classes one thus derives from Collector a class (with
 * access to the actual data), which includes a collect method, which
 * will be called and passed the time at which the sample must be
 * taken.  Then one creates one of the Sampler derived classes
 * (described below), passing the constructor the Collector object.
 * Then, during the simulation, one calls Sampler::sample(time) after
 * the new time is computed (but *before* changing the configuration).
 * See the description of the Sampler classes for details of when the
 * sample is taken.
 *
 * Example code:
 *
 * SimulationData SimData;          // structure containing the data to be sampled
 * MyCollector collector(&SimData); // Derived from collector, will read
 *                                  // the data from the data structure
 * Gillespie_sampler sampler(0,tmax,1,&collector); // One of the Sampler descendats
 * double time=0;
 * while (time<=tmax) {
 *   (compute transition probabilities)
 *   (advance time)
 *   sampler.sample(time);
 *   (perform transition, updating SimData)
 * } 
 *
 */


#ifndef ESAMPLER_HH
#define ESAMPLER_HH

#include <iostream>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
//
// class Collector
//
// Abstract base class for objects that collect the actual data and
// perform output, compute statistical averages, etc

class Collector {
public:
  virtual void collect(double time)=0;
} ;


///////////////////////////////////////////////////////////////////////////////
//
// Sampler
//
// Abstract base class for the different types of samplers (currently
// Passthrough_sampler and Gillespie_sampler).  The run() method of
// the simulaiton should use an objecto of this abstract type and call
// Sampler::sample(time) as shown in the example above.

class Sampler {
public:
  Sampler(Collector *collector);
  virtual void sample(double time)=0;
  virtual ~Sampler() {}

protected:
  Collector* collector;
} ;
  
inline Sampler::Sampler(Collector *collector) :
  collector(collector)
{}

///////////////////////////////////////////////////////////////////////////////
//
// class Passthrough_sampler
//
// This class does not do any sampling, it calls the collector every
// time it gets called.  The only manipulation is to pass to the
// collector the time of the previous call, so that times will be
// recorded correctly when the sampler is called as in the example
// code above.


class Passthrough_sampler : public Sampler {
public:
  Passthrough_sampler(double t0,double tmax,Collector *collector);
  void sample(double time);

private:
  double tprev;
  double tmax;
} ;

inline Passthrough_sampler::Passthrough_sampler(double t0,double tmax,Collector *collector) :
  Sampler(collector),
  tprev(t0),
  tmax(tmax)
{}

inline void Passthrough_sampler::sample(double time)
{
  if (time<=tmax)
    collector->collect(tprev);
  tprev=time;
}

///////////////////////////////////////////////////////////////////////////////
//
// class Gillespie_sampler
//
// This class acts a filter, sending to the collector times sampled at
// regular intervals instead of the random sequence of time generated
// by the Gillespie algorithm.  This means generating less calls to
// the collector if the time intervals of the Gillespie dynamics are
// shorter than the deltat parameter (see the constructor), or more
// calls if the Gillespie time interval is larger than deltat (in this
// last case, the same sample will be appear in the output more than
// once).

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
