/*
 * seeiir_i1.cc
 *
 * An implementation (class SEIRPopulation and main())
 *
 * This is the driver for SEEIIR simulation (main and parameter
 * reading).  Link with a SEEIIR imlementation (seeiir_i1.cc or
 * seeiir_i2.cc) for a funcitoning simulation.
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

  Family(int M=0) : M(M), S(M), E1(0), E2(0), I1(0), I2(0), R(0)
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
  ~SEIRPopulation() { delete ran; delete Mdist;}
  
  void rebuild_families();
  void set_all_S();            // reset all individuals to S
  void compute_rates();        // compute all rates from scratch
  void event(int f,double r);  // perform an event in family f, update family's rates
  void add_imported(int I);    // add infected (E1) at random so that the number
                               // of imported cases becomes I
  void set_beta_out(double b); // call to change beta_out during simulation (invalidates rates)
  
  double              total_rate;
  SEEIIRistate        gstate;
  std::vector<double> cumrate;

private:
  double beta_in,beta_out,sigma,gamma;
  int    NFamilies,Mmax;
  Uniform_integer       *ran;
  Discrete_distribution *Mdist;

  std::vector<Family> families;

  std::vector<int>    suscount;
  
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
  ran = new Uniform_integer;
  Mdist = new Discrete_distribution(Mmax+1,P);
  rebuild_families();
}

void SEIRPopulation::rebuild_families()
{
  gstate.N=gstate.S=gstate.E1=gstate.E2=gstate.I1=gstate.I2=0;gstate.R=0;
  gstate.inf_close=gstate.inf_community=gstate.inf_imported=0;

  Family fam;
  families.clear();

  for (int f=0; f<NFamilies; ++f) {
    fam.M=(*Mdist)();
    fam.S=fam.M;
    gstate.N+=fam.M;
    gstate.S+=fam.S;
    families.push_back(fam);
  }
  gstate.beta_out=beta_out;
  cumrate.resize(families.size()+1,0.);
  suscount.resize(families.size()+1,0.);
}

void SEIRPopulation::set_all_S()
{
  gstate.S=gstate.N;
  gstate.E1=gstate.E2=gstate.I1=gstate.I2=0;gstate.R=0;
  gstate.inf_close=gstate.inf_community=gstate.inf_imported=0;

  for (int f=0; f<families.size(); ++f) {
    families[f].S=families[f].M;
    families[f].E1=families[f].E2=families[f].I1=families[f].I2=0;
  }
  cumrate.resize(families.size()+1,0.);
  suscount.resize(families.size()+1,0.);
}

// call to change beta_out during simulation (invalidates rates)
inline void SEIRPopulation::set_beta_out(double beta)
{
  beta_out=beta;
  gstate.beta_out=beta;
}

/*
 * compute_rates() computes the cumulative rates of all events
 *
 */
void SEIRPopulation::compute_rates()
{
  int N1=gstate.N-1;
  int f;
  cumrate[0]=0;
  for (int f=0; f<families.size(); ++f) {
    double famrate =families[f].S *
      ( beta_in*(families[f].I1 + families[f].I2) + beta_out*(gstate.I1+gstate.I2)/N1 ) +
      families[f].E1 * 2*sigma +
      families[f].E2 * 2*sigma +
      families[f].I1 * 2*gamma +
      families[f].I2 * 2*gamma;
    cumrate[f+1]=cumrate[f]+famrate;
  }
  total_rate=cumrate[families.size()];
}

void SEIRPopulation::event(int f,double r)   // Do infection or recovery on family f
                                             // and update the family rate
{
  double rE1E2=families[f].E1 * 2*sigma;
  double rE2I1=rE1E2 + families[f].E2 * 2*sigma;
  double rI1I2=rE2I1 + families[f].I1 * 2*gamma;
  double rI2R=rI1I2 + families[f].I2 * 2*gamma;
  double localrate=rI2R + families[f].S * beta_in * (families[f].I1 + families[f].I2);

  if (r<rE1E2) {
    families[f].E1--;
    gstate.E1--;
    families[f].E2++;
    gstate.E2++;
  } else if (r<rE2I1) {
    families[f].E2--;
    gstate.E2--;
    families[f].I1++;
    gstate.I1++;
  } else if (r<rI1I2) {
    families[f].I1--;
    gstate.I1--;
    families[f].I2++;
    gstate.I2++;
  } else if (r<rI2R) {
    families[f].I2--;
    gstate.I2--;
    families[f].R++;
    gstate.R++;
  } else {
    families[f].S--;
    gstate.S--;
    families[f].E1++;
    gstate.E1++;
    if (r<localrate) gstate.inf_close++;
    else gstate.inf_community++;
  }    
}

void SEIRPopulation::add_imported(int I)
{
  I-=gstate.inf_imported;  // This is the number of new cases
  if (I>gstate.S)
    {std::cerr << "Cannot add imported, too few suscetibles\n"; exit(1);}
  if (I<0)
    {std::cerr << "Error in imported infections file: external infections must be monotnically increasing\n"; exit(1);}

  // to choose in which family to infect, build cumulative
  // count of susceptibles and draw random number
  suscount[0]=0;
  for (int i=0; i<families.size(); ++i)
    suscount[i+1]=suscount[i]+families[i].S;
  if (suscount.back()!=gstate.S)
    {std::cerr << "ugh:  gstate.S " << gstate.S << " but  " << suscount.back() <<"\n";exit(1);}

  // Randomly choose I disticnt individuals
  std::vector<int> inf;
  for (int i=0; i<I; ++i) {
    int Sn;
    do Sn=(*ran)(gstate.S)+1;
    while ( std::find(inf.begin(),inf.end(),Sn)!=inf.end() ) ;
    inf.push_back(Sn);
  }

  // find them in the families and infect
  for (auto Sn: inf) {
    int f=bsearch(Sn,suscount);
    families[f].S--;
    families[f].E1++;
  }
  gstate.S-=I;
  gstate.E1+=I;
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
  double last=-1;

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
      int f=bsearch(r,pop.cumrate);      // choose family
      pop.event(f,r-pop.cumrate[f]);     // choose and apply event

    }

    if (time>=last+1.) {
      last=time;
    // print or accumulate average
      state->push(time,pop.gstate);
    }

  }
}
