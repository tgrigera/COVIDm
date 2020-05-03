/*
 * popstate.cc -- structures to record, print and average the
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

#include <iostream>
#include <numeric>

#include "popstate.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Population_state

static void print_detail(std::ostream&,const std::vector<double>&);

const char *Population_state::header()
{
  hdr.clear();
  hdr="#      time ";
  addSIRhdr();
  return hdr.c_str();
}  

void Population_state::addSIRhdr()
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

void Population_state::create_colnums(int nc)
{
  char buf[30];
  colnums="#     ( 1)|";
  for (int i=2; i<=nc; ++i) {
    sprintf(buf," |     (%2d)|",i);
    colnums+=buf;
  }
  colnums+="\n";
}

void Population_state::print(std::ostream& o,bool print_time)
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

///////////////////////////////////////////////////////////////////////////////
//
// SIRstate

void SIRstate::push(double time_,SIRistate &s)
{
  time=time_;
  S[0]=s.S;
  I[0]=s.I;
  R[0]=s.R;

  std::cout << *this << '\n';
}

///////////////////////////////////////////////////////////////////////////////
//
// SIRstate_av

const char *SIRstate_av::header()
{
  hdr.clear();
  hdr="#           |----------- Average -------------| |------------ Variance -----------|\n";
  hdr+="#      time ";
  addSIRhdr();
  hdr+=" ";
  addSIRhdr();
  return hdr.c_str();
}  


void SIRstate_av::push(double time,SIRistate &s)
{
  Sav.push(time,s.S);
  Iav.push(time,s.I);
  Rav.push(time,s.R);
}

void SIRstate_av::print(std::ostream& o,bool print_time)
{
  std::vector<double> tim,Sa,Sv,Ia,Iv,Ra,Rv;
  Sav.get_aves(tim,Sa,Sv);
  Iav.get_aves(tim,Ia,Iv);
  Rav.get_aves(tim,Ra,Rv);

  for (int i=0; i<Sv.size(); ++i) {
    time=tim[i];
    S[0]=Sa[i];
    I[0]=Ia[i];
    R[0]=Ra[i];
    Population_state::print(o,print_time);
    S[0]=Sv[i];
    I[0]=Iv[i];
    R[0]=Rv[i];
    Population_state::print(o,false);
    o << '\n';
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRstate

const char *SEEIIRstate::header()
{
  hdr.clear();
  create_colnums(14);
  hdr=colnums;
  hdr+="#      time ";
  addSIRhdr();
  hdr+="   Imported  CloseCntct   Community           N    beta_out";
  return hdr.c_str();
}  

void SEEIIRstate::push(double time_,SEEIIRistate &s)
{
  time=time_;
  S[0]=s.S;
  E[0]=s.E1;
  E[1]=s.E2;
  I[0]=s.I1;
  I[1]=s.I2;
  R[0]=s.R;
  Population_state::print(std::cout,true);
  char buf[200];
  sprintf(buf,"%11d %11d %11d %11d %11.6g",s.inf_imported,s.inf_close,s.inf_community,
	  s.N,s.beta_out);
  std::cout << buf << '\n';
}

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRstate_av

const char *SEEIIRstate_av::header()
{
  hdr.clear();
  create_colnums(27);
  hdr=colnums;
  hdr+="#           |------------------------------------------------------------------------ Average ------------------------------------------------------------------------| |----------------------------------------------------------------------- Variance ------------------------------------------------------------------------|\n";
  hdr+="#      time ";
  addSIRhdr();
  hdr+="   Imported  CloseCntct   Community           N    beta_out";
  addSIRhdr();
  hdr+="    Imported  CloseCntct   Community           N    beta_out";
  return hdr.c_str();
}  

void SEEIIRstate_av::push(double time,SEEIIRistate &s)
{
  Nav.push(time,s.N);
  Sav.push(time,s.S);
  E1av.push(time,s.E1);
  E2av.push(time,s.E2);
  I1av.push(time,s.I1);
  I2av.push(time,s.I2);
  Rav.push(time,s.R);
  Impav.push(time,s.inf_imported);
  Closeav.push(time,s.inf_close);
  Commav.push(time,s.inf_community);
  betaav.push(time,s.beta_out);
}

void SEEIIRstate_av::print(std::ostream& o,bool print_time)
{
  std::vector<double> tim,Sa,Sv,E1a,E1v,E2a,E2v,I1a,I1v,I2a,I2v,Ra,Rv,
    impa,impv,closea,closev,comma,commv,betaa,betav,Na,Nv;
  Nav.get_aves(tim,Na,Nv);
  Sav.get_aves(tim,Sa,Sv);
  E1av.get_aves(tim,E1a,E1v);
  E2av.get_aves(tim,E2a,E2v);
  I1av.get_aves(tim,I1a,I1v);
  I2av.get_aves(tim,I2a,I2v);
  Rav.get_aves(tim,Ra,Rv);
  Impav.get_aves(tim,impa,impv);
  Closeav.get_aves(tim,closea,closev);
  Commav.get_aves(tim,comma,commv);
  betaav.get_aves(tim,betaa,betav);

  char buf[200];
  for (int i=0; i<Sv.size(); ++i) {
    time=tim[i];
    S[0]=Sa[i];
    E[0]=E1a[i];
    E[1]=E2a[i];
    I[0]=I1a[i];
    I[1]=I2a[i];
    R[0]=Ra[i];
    Population_state::print(o,print_time);
    sprintf(buf,"%11.6g %11.6g %11.6g %11.6g %11.6g",impa[i],closea[i],comma[i],Na[i],betaa[i]);
    std::cout << buf;

    S[0]=Sv[i];
    E[0]=E1v[i];
    E[1]=E2v[i];
    I[0]=I1v[i];
    I[1]=I2v[i];
    R[0]=Rv[i];
    Population_state::print(o,false);
    sprintf(buf," %11.6g %11.6g %11.6g %11.6g %11.6g",impv[i],closev[i],commv[i],Nv[i],betav[i]);
    std::cout << buf << '\n';
  }
}
