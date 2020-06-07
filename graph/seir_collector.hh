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

#include "../geoave.hh"
#include "esampler.hh"
#include "sirmodel.hh"
#include "seirmodel.hh"

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

template <typename EGraph>
class SIRcollector : public SEIRcollector_base {
public:
  typedef SIR_model<EGraph> model_t;
  
  SIRcollector(model_t &model) :
    SEIRcollector_base(1,0,1,1),
    model(model)
  {}
  void collect(double time);

protected:
  SIR_model<EGraph> &model;
} ;

template <typename EGraph>
void SIRcollector<EGraph>::collect(double time_)
{
  time=time_;
  typename model_t::aggregate_node *anode=model.anodemap[model.hroot];
  S[0]=anode->NS;
  I[0]=anode->NI;
  R[0]=anode->NR;

  std::cout << *this << '\n';
}

///////////////////////////////////////////////////////////////////////////////
//
// SIRcollector_av

template <typename EGraph>
class SIRcollector_av : public SIRcollector<EGraph> {
public:
  typedef SIR_model<EGraph> model_t;

  SIRcollector_av(SIR_model<EGraph> &model,double deltat=1.) :
    SIRcollector<EGraph>(model),
    Sav(-0.5*deltat,1.,deltat),
    Iav(-0.5*deltat,1.,deltat),
    Rav(-0.5*deltat,1.,deltat)
  {}

  const char* header();
  void print(std::ostream&,bool print_time=true);
  void collect(double time);

private:
  Geoave Sav,Iav,Rav;
  using SEIRcollector_base::hdr;
  using SEIRcollector_base::time;
  using SEIRcollector_base::S;
  using SEIRcollector_base::I;
  using SEIRcollector_base::R;
  using SEIRcollector_base::addSIRhdr;
  using SIRcollector<EGraph>::model;
} ;

template <typename EGraph>
const char *SIRcollector_av<EGraph>::header()
{
  hdr.clear();
  hdr="#           |----------- Average -------------| |------------ Variance -----------|\n";
  hdr+="#      time ";
  addSIRhdr();
  hdr+=" ";
  addSIRhdr();
  return hdr.c_str();
}  

template <typename EGraph>
void SIRcollector_av<EGraph>::collect(double time)
{
  typename model_t::aggregate_node *anode=model.anodemap[model.hroot];
  Sav.push(time,anode->NS);
  Iav.push(time,anode->NI);
  Rav.push(time,anode->NR);
}

template <typename EGraph>
void SIRcollector_av<EGraph>::print(std::ostream& o,bool print_time)
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
    SEIRcollector_base::print(o,print_time);
    S[0]=Sv[i];
    I[0]=Iv[i];
    R[0]=Rv[i];
    SEIRcollector_base::print(o,false);
    o << '\n';
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRcollector

template <typename EGraph>
class SEEIIRcollector : public SEIRcollector_base {
public:
  typedef SEEIIR_model<EGraph> model_t;
  
  SEEIIRcollector(model_t &model) :
    SEIRcollector_base(1,2,2,1),
    model(model)
  {}
  void collect(double time);

protected:
  SEEIIR_model<EGraph> &model;
} ;

template <typename EGraph>
void SEEIIRcollector<EGraph>::collect(double time_)
{
  time=time_;
  typename model_t::aggregate_data *anode=model.anodemap[model.hroot];
  S[0]=anode->NS;
  E[0]=anode->NE1;
  E[1]=anode->NE2;
  I[0]=anode->NI1;
  I[1]=anode->NI2;
  R[0]=anode->NR;

  std::cout << *this << '\n';
}

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIRcollector_av

template <typename EGraph>
class SEEIIRcollector_av : public SEEIIRcollector<EGraph> {
public:
  typedef SEEIIR_model<EGraph> model_t;

  SEEIIRcollector_av(SEEIIR_model<EGraph> &model,double deltat=1.) :
    SEEIIRcollector<EGraph>(model),
    Sav(-0.5*deltat,1.,deltat),
    E1av(-0.5*deltat,1.,deltat),
    E2av(-0.5*deltat,1.,deltat),
    I1av(-0.5*deltat,1.,deltat),
    I2av(-0.5*deltat,1.,deltat),
    Rav(-0.5*deltat,1.,deltat),
    Impav(-0.5*deltat,1.,deltat),
    Closeav(-0.5*deltat,1.,deltat),
    Commav(-0.5*deltat,1.,deltat),
    Nav(-0.5*deltat,1.,deltat)
  {}

  const char* header();
  void print(std::ostream&,bool print_time=true);
  void collect(double time);

private:
  Geoave Sav,E1av,E2av,I1av,I2av,Rav,Impav,Closeav,Commav,Nav;

  using SEIRcollector_base::hdr;
  using SEIRcollector_base::time;
  using SEIRcollector_base::S;
  using SEIRcollector_base::E;
  using SEIRcollector_base::I;
  using SEIRcollector_base::R;
  using SEIRcollector_base::addSIRhdr;
  using SEIRcollector_base::create_colnums;
  using SEIRcollector_base::colnums;
  using SEEIIRcollector<EGraph>::model;
} ;

template <typename EGraph>
const char *SEEIIRcollector_av<EGraph>::header()
{
  hdr.clear();
  create_colnums(23);
  hdr=colnums;
  hdr+="#           |----------------------------------------------------------- Average -------------------------------------------------------------| |------------------------------------------------------------ Variance -----------------------------------------------------------|\n";
  hdr+="#      time ";
  addSIRhdr();
  hdr+="   Imported  CloseCntct   Community";
  addSIRhdr();
  hdr+="    Imported  CloseCntct   Community";
  return hdr.c_str();
}  

template <typename EGraph>
void SEEIIRcollector_av<EGraph>::collect(double time_)
{
  typename model_t::aggregate_data *anode=model.anodemap[model.hroot];
  time=time_;
  Nav.push(time_,anode->Ntot);
  Sav.push(time_,anode->NS);
  E1av.push(time_,anode->NE1);
  E2av.push(time_,anode->NE2);
  I1av.push(time_,anode->NI1);
  I2av.push(time_,anode->NI2);
  Rav.push(time_,anode->NR);
  Impav.push(time,anode->Nimported);
  Closeav.push(time,anode->Nclose);
  Commav.push(time,anode->Ncommunity);

  S[0]=anode->NS;
  E[0]=anode->NE1;
  E[1]=anode->NE2;
  I[0]=anode->NI1;
  I[1]=anode->NI2;
  R[0]=anode->NR;
}

template <typename EGraph>
void SEEIIRcollector_av<EGraph>::print(std::ostream& o,bool print_time)
{
  std::vector<double> tim,Sa,Sv,E1a,E1v,E2a,E2v,I1a,I1v,I2a,I2v,Ra,Rv,
    impa,impv,closea,closev,comma,commv,Na,Nv;
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

  char buf[200];
  for (int i=0; i<Sv.size(); ++i) {
    time=tim[i];
    S[0]=Sa[i];
    E[0]=E1a[i];
    E[1]=E2a[i];
    I[0]=I1a[i];
    I[1]=I2a[i];
    R[0]=Ra[i];
    SEIRcollector_base::print(o,print_time);
    sprintf(buf,"%11.6g %11.6g %11.6g %11.6g",impa[i],closea[i],comma[i],Na[i]);
    std::cout << buf;

    S[0]=Sv[i];
    E[0]=E1v[i];
    E[1]=E2v[i];
    I[0]=I1v[i];
    I[1]=I2v[i];
    R[0]=Rv[i];
    SEIRcollector_base::print(o,false);
    sprintf(buf," %11.6g %11.6g %11.6g %11.6g",impv[i],closev[i],commv[i],Nv[i]);
    std::cout << buf << '\n';
  }
}


#endif /* SEIR_COLLECTOR_HH */
