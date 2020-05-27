/*
 * eevents.hh -- structures for external events
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

#ifndef EEVENTS_HH
#define EEVENTS_HH

#include <queue>
#include <lemon/concepts/digraph.h>

class Event {
public:
  Event(double time) : time(time)  {}
  virtual void apply(Epidemiological_model*) {}

  double  time;
} ;

class Imported_infection : public Event {
public:
  Imported_infection(double time,int new_cases) :
    Event(time), new_cases(new_cases) {}
  void apply(Epidemiological_model*);

  int  new_cases;
} ;

inline void Imported_infection::apply(Epidemiological_model* em)
{
  em->add_imported_infections(this);
}

typedef std::queue<Event*> event_queue_t;


#endif /* EEVENTS_HH */
