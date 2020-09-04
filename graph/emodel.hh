/*
 * emodel.hh -- classes for epidemiological models on graphs
 *
 * The classes here implement different epidemiological models.  They
 * are coded generically and can be used over any graph provided by
 * egraph.hh.
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

#ifndef EMODEL_HH
#define EMODEL_HH

#include "esampler.hh"
#include "egraph.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Epidemiological_model

class Forced_transition;

class Epidemiological_model {
public:
  virtual ~Epidemiological_model() {}
  virtual void apply_transition(int)=0;
  virtual void compute_all_rates()=0;
  virtual void set_all_susceptible()=0;
  virtual void add_imported(Forced_transition*)=0;
  void         update_cumulative_rates();

  std::vector<double> cumulative_rates;

protected:
  struct transition {
    int           nodeid;
    double        rate;
    short         type;
  } ;
  std::vector<transition> transitions;

} ;

#include "eevents.hh"

template <typename EGraph>
class Epidemiological_model_graph_base : public Epidemiological_model {
public:
  Epidemiological_model_graph_base(EGraph &egraph);
  void         compute_all_rates();
  virtual void compute_rates(typename EGraph::igraph_t::Node)=0;

protected:
  EGraph&  egraph;

} ;

template <typename EGraph>
inline Epidemiological_model_graph_base<EGraph>::Epidemiological_model_graph_base(EGraph& egraph) :
  egraph(egraph)
{}

template <typename EGraph>
void Epidemiological_model_graph_base<EGraph>::compute_all_rates()
{
  for (typename EGraph::igraph_t::NodeIt node(egraph.igraph); node!=lemon::INVALID; ++node)
    compute_rates(node);
  update_cumulative_rates();
}

void run(Epidemiological_model *model,Sampler* sampler,event_queue_t& events,double tmax);


#endif /* EMODEL_HH */
