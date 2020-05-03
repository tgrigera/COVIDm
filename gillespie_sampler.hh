/*
 * gillespie_sampler.hh
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

#ifndef GILLESPIE_SMPLER_HH
#define GILLESPIE_SAMPLER_HH

///////////////////////////////////////////////////////////////////////////////
//
// Gillespie_sampler
//
// This is a simple template class which acts as a filter between the
// data generation and actual data output or averaging.  It's function
// is to send to output the state of the system at regular intervals,
// adequately sampling data generated with varying time step by the
// Gillespie algorithm.
//
// Template parameters are Pusher (a class with a method push(double
// time,Data &data) that does the output or average, and Data, the
// class/struct that the Pusher class expects to receive the data.
//
// Usage: create a Gillespie_sampler object giving the pusher, start
// time, end time and desired incremet.  Then there are two methods to
// call: push_time and push_data.
//
// You must call push_time with the time for the next step *before*
// the data has been updated, then call push_data() giving the data
// for the new time.  This is because push_data just stores a pointer
// to the new data, but cannot know whether it must be output before
// knowing what the next Gillespie time is.
//
// The following is an example of how to use Gillespie_sampler; see also sir_m.cc
//
// void run_gillespie(SIRstate *state)
// {
//   Gillespie_sampler<SIRstate,SIRistate> gsamp(*state,0,options.steps,1);
//
//   SIRistate gstate;
//   (load data for starting state)
//   gsamp.push_data(gstate);
//
//   double time=0;
//   while (time<max_time) {
//     (compute transition probabilities)
//     (advance time)
//     gsamp.push_time(time);
//     (perform transition)
//     (load new data)
//     gsamp.push_data(gstate);
//   }
// }


template <typename Pusher,typename Data>
class Gillespie_sampler {
public:
  Gillespie_sampler(Pusher& pusher,double t0,double tmax,double deltat);
  void push_time(double time);
  void push_data(Data& d) {data=&d;}
  
private:
  Pusher &pusher;
  double t0,deltat,tmax;
  double tlast,tnext;
  Data*  data;
} ;

template <typename Pusher,typename Data>
Gillespie_sampler<Pusher,Data>::Gillespie_sampler(Pusher& pusher,double t0,double tmax,double deltat) :
  pusher(pusher),
  t0(t0),
  deltat(deltat),
  tmax(tmax),
  tlast(t0-deltat),
  tnext(t0),
  data(0)
{
}

template <typename Pusher,typename Data>
void Gillespie_sampler<Pusher,Data>::push_time(double time)
{
  if (time<=tnext || tlast>=tmax) return;
  double t;
  for (t=tlast+=deltat; t<time && t<=tmax; t+=deltat) 
    pusher.push(t,*data);
  tnext=t;
  tlast=t-deltat;
}

#endif /* GILLESPIE_SAMPLER_HH */
