/*
 * seeiir_i3.cc
 *
 * An implementation (class SEIRPopulation and main())
 *
 * Stochastic SEIR with two E and to I states.  Population separated
 * in families.  Simulated in continuous time (Gillespie algorithm).
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

#include <assert.h>

///////////////////////////////////////////////////////////////////////////////
//
// Simulation

/*
 * Global and family state
 *
 */

struct Family {
  int  M;
  int  S,E1,E2,I1,I2,R;
  int  first_susceptible_in_list;
  int  infected_families_in_list;

  Family(int M=0) : M(M), S(M), E1(0), E2(0), I1(0), I2(0), R(0),
		    first_susceptible_in_list(-1), infected_families_in_list(-1)
  {}
} ;

std::ostream& operator<<(std::ostream& o,const Family &f)
{
  std::cout << "M = " << f.M  << " S,E1,E2,I1,I2,R = "
	    << f.S << ' ' << f.E1 << ' ' << f.E2 << ' '
	    << f.I1 << ' ' << f.I2 << ' ' << f.R;
  return o;
}

/*
 * class SEIRPopulation holds per family information, computes rates and
 * performs individual state switchs (function event)
 *
 */
class SEIRPopulation {
public:
  SEIRPopulation(int NFamilies,double beta_in,double beta_out,double sigma,
		 double gamma,int Mmax, double PM[]);
  ~SEIRPopulation();
  
  void rebuild_families();
  void set_all_S();            // reset all individuals to S
  void compute_rates();        // compute all rates from scratch
  void event(int f,double r);  // perform an event in family f, update family's rates
  void add_imported(int I);    // add infected (E1) at random so that the number
                               // of imported cases becomes I
  void set_beta_out(double b); // call to change beta_out during simulation (invalidates rates)
  
  void local_infection(int fn);
  void global_infection();
  void E1E2();
  void E2I1();
  void I1I2();
  void I2R();

  SEEIIRistate        gstate;
  std::vector<double> cumrate;
  double              total_rate;

  std::vector<int>    families_infected;

private:
  double beta_in,beta_out,sigma,gamma;
  int    NFamilies,Mmax;
  Uniform_integer       *ran;;
  Discrete_distribution *Mdist;

  std::vector<Family*>     families;
  std::vector<int>         listS,listE1,listE2,listI1,listI2;

  void erase_susceptible(int family_number);
} ;
		 
SEIRPopulation::SEIRPopulation(int NFamilies,double beta_in,double beta_out,double sigma,
			       double gamma,int Mmax, double *P) :
  beta_in(beta_in),
  beta_out(beta_out),
  sigma(sigma),
  gamma(gamma),
  NFamilies(NFamilies),
  Mmax(Mmax),
  ran(0),
  Mdist(0)
{
  Mdist = new Discrete_distribution(Mmax+1,P);
  ran = new Uniform_integer;
  rebuild_families();
}

SEIRPopulation::~SEIRPopulation()
{
  for (auto f: families) delete f;
  delete ran;
  delete Mdist;
}

void SEIRPopulation::rebuild_families()
{
  gstate.N=gstate.S=gstate.E1=gstate.E2=gstate.I1=gstate.I2=gstate.R=0;
  gstate.inf_close=gstate.inf_community=gstate.inf_imported=0;
  gstate.beta_out=beta_out;

  Family *fam;
  families.clear();
  families_infected.clear();
  listS.clear();
  listE1.clear();
  listE2.clear();
  listI1.clear();
  listI2.clear();

  for (int f=0; f<NFamilies; ++f) {
    fam = new Family( (*Mdist)() );
    gstate.N+=fam->M;
    gstate.S+=fam->S;
    fam->infected_families_in_list=-1;
    families.push_back(fam);
    fam->first_susceptible_in_list=listS.size();
    for (int iS=0; iS<fam->S; iS++)
      listS.push_back(families.size()-1);
  }
}

void SEIRPopulation::set_all_S()
{
  gstate.S=gstate.N;
  gstate.E1=gstate.E2=gstate.I1=gstate.I2=0;gstate.R=0;
  gstate.inf_close=gstate.inf_community=gstate.inf_imported=0;

  families_infected.clear();;
  listS.clear();
  listE1.clear();
  listE2.clear();
  listI1.clear();
  listI2.clear();
  for (int fn=0; fn<families.size(); ++fn) {
    Family* f=families[fn];
    f->S=f->M;
    f->E1=f->E2=f->I1=f->I2=0;
    f->infected_families_in_list=-1;
    f->first_susceptible_in_list=listS.size();
    for (int iS=0; iS<f->S; iS++)
      listS.push_back(fn);
  }
}

// call to change beta_out during simulation (invalidates rates)
inline void SEIRPopulation::set_beta_out(double beta)
{
  beta_out=beta;
  gstate.beta_out=beta;
}

/*
 * compute_rates()
 *
 * This computes the cumulative rates of all events and updates the
 * first_susceptible_in_list field of the families vector
 *
 */
void SEIRPopulation::compute_rates()
{
  int N1=gstate.N-1;
  int f;
  double cr=0;
  
  cumrate.clear();
  cumrate.reserve(gstate.I1+gstate.I2+10);
  cumrate.push_back(0.);

  // First: local infections (S->E1)
  for (int fn: families_infected ) {
    cr += families[fn]->S * beta_in*( families[fn]->I1 + families[fn]->I2);
    cumrate.push_back(cr);
  }
  // Global infections (S->E1)
  cr += gstate.S*beta_out*(gstate.I1+gstate.I2)/(gstate.N-1) ;
  cumrate.push_back(cr);
  // E1->E2
  cr +=  gstate.E1*2*sigma;
  cumrate.push_back(cr);
  // E2->I1
  cr +=  gstate.E2*2*sigma;
  cumrate.push_back(cr);
  // I1->I2
  cr += gstate.I1*2*gamma;
  cumrate.push_back(cr);
  // I2->R
  cr += gstate.I2*2*gamma;
  cumrate.push_back(cr);
  
  total_rate=cumrate.back();
}

/*
 * erase a susceptible from family f and update lists
 *
 */
void SEIRPopulation::erase_susceptible(int family_number)
{
  Family* f=families[family_number];
  f->S--;
  gstate.S--;
  auto Sp=listS.begin() + f->first_susceptible_in_list;
  listS.erase(Sp);

  // update pointers to susceptibles in the following families
  // the pointer will be invalid if the family has no susceptibles,
  // but in that case it should never be looked up any way
  for (int fn=family_number+1; fn<families.size(); ++fn)
    families[fn]->first_susceptible_in_list--;
}

/*
 * methods to apply particular events
 *
 */
void SEIRPopulation::local_infection(int fn) {
  families[fn]->E1++;
  gstate.E1++;
  gstate.inf_close++;

  listE1.push_back(fn);
  erase_susceptible(fn);
}

void SEIRPopulation::global_infection() {
  int Si=(*ran)(gstate.S);         // choose a susceptible with equal probability
  int fn=listS[Si];                // find its family
  families[fn]->E1++;
  gstate.E1++;

  listE1.push_back(fn);
  erase_susceptible(fn);
}

void SEIRPopulation::E1E2() {
  int Ei=(*ran)(gstate.E1);         // choose E1 with equal probability
  auto Ep=listE1.begin()+Ei;        // find its family
  int  fn=*Ep;
  families[fn]->E1--;               // update family and gstate
  families[fn]->E2++;
  gstate.E1--;
  gstate.E2++;
  listE1.erase(Ep);                  // update lists
  listE2.push_back(fn);
}

void SEIRPopulation::E2I1() {
  int Ei=(*ran)(gstate.E2);          // choose E2 with equal probability
  auto Ep=listE2.begin()+Ei;         // find its family
  int  fn=*Ep;
  if (families[fn]->I1+families[fn]->I2+families[fn]->R>0) gstate.inf_close++;
  else gstate.inf_community++;

  families[fn]->E2--;                // update family and gstate
  families[fn]->I1++;
  gstate.E2--;
  gstate.I1++;
  listE2.erase(Ep);                  // update lists
  listI1.push_back(fn);

  if (families[fn]->I1 + families[fn]->I2 == 1) {          // first infection in family, record list
    families[fn]->infected_families_in_list=families_infected.size();
    families_infected.push_back(fn);
  }
}

void SEIRPopulation::I1I2() {
  int Ei=(*ran)(gstate.I1);          // choose E2 with equal probability
  auto Ep=listI1.begin()+Ei;         // find its family
  int  fn=*Ep;
  families[fn]->I1--;                           // update family and gstate
  families[fn]->I2++;
  gstate.I1--;
  gstate.I2++;
  listI1.erase(Ep);                  // update lists
  listI2.push_back(fn);
}

void SEIRPopulation::I2R() {
  int Ei=(*ran)(gstate.I2);          // choose E2 with equal probability
  auto Ep=listI2.begin()+Ei;         // find its family
  int  fn=*Ep;
  families[fn]->I2--;                           // update family and gstate
  families[fn]->R++;
  gstate.I2--;
  gstate.R++;
  listI2.erase(Ep);                  // update lists

  if (families[fn]->I1 + families[fn]->I2 == 0) {          // no more infections, remove family from list
    auto fip=families_infected.begin() + families[fn]->infected_families_in_list;
    families[fn]->infected_families_in_list=-1;
    for (auto fip1=fip+1; fip1!=families_infected.end(); ++fip1)
      families[*fip1]->infected_families_in_list--;
    families_infected.erase(fip);
  }
}

void SEIRPopulation::add_imported(int I)
{
  I-=gstate.inf_imported;  // This is the number of new cases

  if (I>gstate.S)
    {std::cerr << "Cannot add imported, too few suscetibles\n"; exit(1);}
  if (I<0)
    {std::cerr << "Error in imported infections file: external infections must be monotnically increasing\n"; exit(1);}

  // Randomly choose and infect I individuals
  for (int infn=0; infn<I; ++infn) {
    int Sn=(*ran)(gstate.S);
    // find in family and infect in state I1
    int fn=listS[Sn];        // find its family
    families[fn]->I1++;
    gstate.I1++;
    listI1.push_back(fn);
    erase_susceptible(fn);
    if (families[fn]->I1 + families[fn]->I2 == 1) {          // first infection in family, record list
      families[fn]->infected_families_in_list=families_infected.size();
      families_infected.push_back(fn);
    }
  }
  gstate.inf_imported+=I;

}

///////////////////////////////////////////////////////////////////////////////
//
// main simulation driver (Gillespie)

void run(SEIRPopulation &pop,SEEIIRstate *state)
{
  Exponential_distribution rexp;
  Uniform_real ran(0,1.);
  double deltat,time=0;
  double last=-10;

  event_queue_t events=event_queue;

  while (time<options.steps) {

    // compute transition probabilities
    pop.compute_rates();
    double mutot=pop.total_rate;

    // advance time
    deltat=rexp(1./mutot);
    time+=deltat;

    if (time>=events.front().time) {              // imported infections or beta change

      if (events.size()==1) break;
      time=events.front().time;

      switch (events.front().kind) {
      case event::infection:
	pop.add_imported(events.front().I);
	break;
      case event::beta_change:
	pop.set_beta_out(events.front().beta);
	break;
      }
      events.pop();

    } else {

      // choose the transition and apply it
      double r=ran()*mutot;
      int e=bsearch(r,pop.cumrate);      // choose event
      if (e<pop.families_infected.size())       // local infection
	pop.local_infection(pop.families_infected[e]);
      else {                             // global proceses
	e-=pop.families_infected.size();
	switch(e) {
	case 0:	pop.global_infection(); break;
	case 1: pop.E1E2(); break;
	case 2: pop.E2I1(); break;
	case 3: pop.I1I2(); break;
	case 4: pop.I2R(); break;
	}
      }
	
    }

    if (time>=last+1.) {  // print or accumulate average
      last=time;
      state->push(time,pop.gstate);
    }

  }
}
