/*
 * seeiir_i2.cc
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

struct readS {
  static int& element(Family* f) {return f->S;}
} ;
struct readE1 {
  static int& element(Family* f) {return f->E1;}
} ;
struct readE2 {
  static int& element(Family* f) {return f->E2;}
} ;
struct readI1 {
  static int& element(Family* f) {return f->I1;}
} ;
struct readI2 {
  static int& element(Family* f) {return f->I2;}
} ;

template <typename readerT> class Cumulative_count {
public:
  void reset_count(std::vector<Family*> &vec);
  void update_count();
  std::vector<int> count;

private:
  std::vector<Family*> *vec;
} ;

template <typename readerT>
void Cumulative_count<readerT>::reset_count(std::vector<Family*> &v)
{
  vec=&v;
  count.resize(vec->size()+1);
  update_count();
}

template <typename readerT>
void Cumulative_count<readerT>::update_count()
{
  count[0]=0;
  for (int i=0; i<vec->size(); ++i)
    count[i+1]=count[i]+readerT::element((*vec)[i]);
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
  
  void local_infection(int f);
  void global_infection();
  void E1E2();
  void E2I1();
  void I1I2();
  void I2R();

  double beta_in,beta_out,sigma,gamma;
  int    NFamilies,Mmax;
  Uniform_integer       *ran;;
  Discrete_distribution *Mdist;

  SEEIIRistate             gstate;
  std::vector<Family*>     families;
  int                      families_infected;
  Cumulative_count<readS>  cumS;
  Cumulative_count<readE1> cumE1;
  Cumulative_count<readE2> cumE2;
  Cumulative_count<readI1> cumI1;
  Cumulative_count<readI2> cumI2;

  std::vector<double> cumrate;
  double              total_rate;
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

  Family *fam;
  families.clear();

  for (int f=0; f<NFamilies; ++f) {
    fam = new Family( (*Mdist)() );
    gstate.N+=fam->M;
    gstate.S+=fam->S;
    families.push_back(fam);
  }
  families_infected=0;

  cumS.reset_count(families);
  cumE1.reset_count(families);
  cumE2.reset_count(families);
  cumI1.reset_count(families);
  cumI2.reset_count(families);
}

void SEIRPopulation::set_all_S()
{
  gstate.S=gstate.N;
  gstate.E1=gstate.E2=gstate.I1=gstate.I2=0;gstate.R=0;
  gstate.inf_close=gstate.inf_community=gstate.inf_imported=0;

  for (int f=0; f<families.size(); ++f) {
    families[f]->S=families[f]->M;
    families[f]->E1=families[f]->E2=families[f]->I1=families[f]->I2=0;
  }
  families_infected=0;

  cumS.update_count();
  cumE1.update_count();
  cumE2.update_count();
  cumI1.update_count();
  cumI2.update_count();
}

/*
 * compute_rates() computes the cumulative rates of all events
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
  for (f=0; f<families_infected; ++f) {
    cr += families[f]->S * beta_in*(families[f]->I1 + families[f]->I2);
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
 * methods to apply particular events
 *
 */
void SEIRPopulation::local_infection(int f) {
  families[f]->S--;
  families[f]->E1++;
  gstate.S--;
  gstate.E1++;

  cumS.update_count();
  cumE1.update_count();
  gstate.inf_close++;
}

void SEIRPopulation::global_infection() {
  int Si=(*ran)(gstate.S)+1;       // choose a susceptible with equal probability
  int f=bsearch(Si,cumS.count);    // find its family
  families[f]->S--;
  families[f]->E1++;
  gstate.S--;
  gstate.E1++;

  cumS.update_count();
  cumE1.update_count();
  gstate.inf_community++;
}

void SEIRPopulation::E1E2() {
  int Ei=(*ran)(gstate.E1)+1;          // choose E1 with equal probability
  int f=bsearch(Ei,cumE1.count);       // find its family
  families[f]->E1--;
  families[f]->E2++;
  gstate.E1--;
  gstate.E2++;

  cumE1.update_count();
  cumE2.update_count();
}

void SEIRPopulation::E2I1() {
  int Ei=(*ran)(gstate.E2)+1;      // choose E2 with equal probability
  int f=bsearch(Ei,cumE2.count);   // find its family
  families[f]->E2--;
  families[f]->I1++;
  gstate.E2--;
  gstate.I1++;

  if (families[f]->I1 + families[f]->I2 == 1) {    // first infection in family, move it up in the list
    std::swap(families[f],families[families_infected]);
    families_infected++;
    cumS.update_count();
    cumE1.update_count();
    cumE2.update_count();
    cumI1.update_count();
    cumI2.update_count();
  } else {
    cumE2.update_count();
    cumI1.update_count();
  }
}

void SEIRPopulation::I1I2() {
  int Ei=(*ran)(gstate.I1)+1;        // choose I1 with equal probability
  int f=bsearch(Ei,cumI1.count);     // find its family
  families[f]->I1--;
  families[f]->I2++;
  gstate.I1--;
  gstate.I2++;

  cumI1.update_count();
  cumI2.update_count();
}

void SEIRPopulation::I2R() {
  int Ei=(*ran)(gstate.I2)+1;       // choose I2 with equal probability
  int f=bsearch(Ei,cumI2.count);    // find its family
  families[f]->I2--;
  families[f]->R++;
  gstate.I2--;
  gstate.R++;

  if (families[f]->I1 + families[f]->I2 == 0) {  // no more infections, move family down in the list
    families_infected--;
    std::swap(families[f],families[families_infected]);
    cumS.update_count();
    cumE1.update_count();
    cumE2.update_count();
    cumI1.update_count();
    cumI2.update_count();
  } else
    cumI2.update_count();
}

void SEIRPopulation::add_imported(int I)
{
  I-=gstate.inf_imported;  // This is the number of new cases
  if (I>gstate.S)
    {std::cerr << "Cannot add imported, too few suscetibles\n"; exit(1);}

  // Randomly choose I disticnt individuals
  std::vector<int> inf;
  for (int i=0; i<I; ++i) {
    int Sn;
    do Sn=(*ran)(gstate.S)+1;
    while ( std::find(inf.begin(),inf.end(),Sn)!=inf.end() ) ;
    inf.push_back(Sn);
  }

  // find them in the families and infect
  for (auto Si: inf) {
    int f=bsearch(Si,cumS.count);    // find its family
    families[f]->S--;
    families[f]->E1++;
  }

  gstate.S-=I;
  gstate.E1+=I;
  gstate.inf_imported+=I;

  cumS.update_count();
  cumE1.update_count();
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

  opt::ei imported_end={std::numeric_limits<double>::max(),0};
  opt::imported_infections_t imported=options.imported_infections;
  imported.push(imported_end);
  double t_next_imported=imported.front().time;

  while (time<options.steps) {

    // compute transition probabilities
    pop.compute_rates();
    double mutot=pop.total_rate;

    // advance time
    deltat=rexp(1./mutot);
    time+=deltat;

    if (time>=t_next_imported) {              // it's time to add new imported infections

      if (imported.size()==1) break;
      time=t_next_imported;
      pop.add_imported(imported.front().I);
      imported.pop();
      t_next_imported = imported.front().time;

    } else {

      // choose the transition and apply it
      double r=ran()*mutot;
      int e=bsearch(r,pop.cumrate);      // choose event
      if (e<pop.families_infected)       // local infection
	pop.local_infection(e);
      else {                             // global proceses
	e-=pop.families_infected;
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
