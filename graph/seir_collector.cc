/*
 * seir_collector.cc -- collector classes for SEIR models
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

#include <numeric>

#include "seir_collector.hh"

///////////////////////////////////////////////////////////////////////////////
//
// SEIRcollector_base

static void print_detail(std::ostream&,const std::vector<double>&);

const char *SEIRcollector_base::header()
{
  hdr.clear();
  hdr="#      time ";
  addSIRhdr();
  return hdr.c_str();
}  

void SEIRcollector_base::addSIRhdr()
{
  char buf[50];
  int  i;

  if (E.size()>0)
    sprintf(buf,"          S           E           I           R ");
  else
    sprintf(buf,"          S           I           R ");
  hdr+=buf;

  if (S.size()>1)
    for (i=1; i<=S.size(); ++i) {
      sprintf(buf,"         S%d ",i);
      hdr+=buf;
    }
  if (E.size()>1)
    for (i=1; i<=E.size(); ++i) {
      sprintf(buf,"         E%d ",i);
      hdr+=buf;
    }
  if (I.size()>1)
    for (i=1; i<=I.size(); ++i) {
      sprintf(buf,"         I%d ",i);
      hdr+=buf;
    }
  if (R.size()>1)
    for (i=1; i<=R.size(); ++i) {
      sprintf(buf,"         R%d ",i);
      hdr+=buf;
    }
}

void SEIRcollector_base::create_colnums(int nc)
{
  char buf[30];
  colnums="#     ( 1)|";
  for (int i=2; i<=nc; ++i) {
    sprintf(buf," |     (%2d)|",i);
    colnums+=buf;
  }
  colnums+="\n";
}

void SEIRcollector_base::print(std::ostream& o,bool print_time)
{
  static char buf[200];
  int    i;

  double St=std::accumulate(S.begin(),S.end(),0.);
  double Et=std::accumulate(E.begin(),E.end(),0.);
  double It=std::accumulate(I.begin(),I.end(),0.);
  double Rt=std::accumulate(R.begin(),R.end(),0.);

  if (print_time) {
    sprintf(buf,"%11.6g ",time);
    o << buf;
  }
  if (E.size()>0)
    sprintf(buf,"%11.6g %11.6g %11.6g %11.6g ",St,Et,It,Rt);
  else
    sprintf(buf,"%11.6g %11.6g %11.6g ",St,It,Rt);
  o << buf;

  print_detail(o,S);
  print_detail(o,E);
  print_detail(o,I);
  print_detail(o,R);
}

static void print_detail(std::ostream& o,const std::vector<double>& v)
{
  if (v.size()<=1) return;

  char buf[50];

  for (int i=0; i<v.size(); ++i) {
    sprintf(buf,"%11.6g ",v[i]);
    o << buf;
  }
}

