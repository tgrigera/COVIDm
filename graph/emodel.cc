/*
 * emodel.cc --  driver for epidemiological models on graphs
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

#include <limits>

#include "emodel.hh"
#include "../qdrandom.hh"
#include "../bsearch.hh"
#include "esampler.hh"


///////////////////////////////////////////////////////////////////////////////
//
// Epidemiological_model

void Epidemiological_model::update_cumulative_rates()
{
  cumulative_rates.clear();
  cumulative_rates.reserve(transitions.size()+1);
  cumulative_rates.push_back(0.);
  double cr=0;
  for (transition& t: transitions) {
    cr += t.rate;
    cumulative_rates.push_back(cr);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// simulation driver: uses a given Epidemiological_model to implement
// Gillespie dynamics.  Output through a Gillespie_sampler object

void run(Epidemiological_model *model,Sampler* sampler,event_queue_t& events,double tmax)
{
  Exponential_distribution rexp;
  Uniform_real ran(0,1.);
  double deltat,time=0;

  event_queue_t levents=events;
  levents.push(new Event(std::numeric_limits<double>::max()));
  model->set_all_susceptible();
  model->compute_all_rates();
  while (time<=tmax) {

    // get transition probability and advance time
    model->update_cumulative_rates();
    double mutot=model->cumulative_rates.back();
    deltat=rexp(1./mutot);
    time+=deltat;

    if (time>=levents.front()->time) {              // external event: imported infections, etc

      time=levents.front()->time;
      sampler->sample(time);
      if (levents.size()==1) {
	delete levents.front();
	levents.pop();
	break;
      }
      levents.front()->apply(model);
      levents.pop();

    } else {

      sampler->sample(time);
      // choose the transition and apply it
      double r=ran()*mutot;
      int e=bsearch(r,model->cumulative_rates);
      model->apply_transition(e);
	
    }
  }
}
