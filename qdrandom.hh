/*
 * qdrandom.hh -- quick and dirty interface to GSL's random numbers
 *
 * This is taken from glsim's interface, with less features (no
 * exceptions, no multiple scopes, no serialization through Boost).
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * Coronel Pringles, 24-jul-2017
 *
 */

#ifndef QDRANDOM_HH
#define QDRANDOM_HH

#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Random_number_generator {
public:
  Random_number_generator(const unsigned long seed=0);
  ~Random_number_generator();
  void set_seed(const unsigned long seed);
  int save(FILE *f);
  int load(FILE *f);
  void save(std::ostream&);
  void load(std::istream&);

  unsigned long raw();
  unsigned long min() const;
  unsigned long max() const;
  unsigned long range() const;

private:
  static gsl_rng *generator;
  static Random_number_generator *glsim_generator;
  
  template<typename T> friend class Random_distribution_base;

}  ;

inline void Random_number_generator::set_seed(const unsigned long seed)
{
  gsl_rng_set(generator,seed);
}

inline unsigned long Random_number_generator::raw()
{
  return gsl_rng_get(generator);
}

inline unsigned long Random_number_generator::min() const
{
  return gsl_rng_min(generator);
}

inline unsigned long Random_number_generator::max() const
{
  return gsl_rng_max(generator);
}

inline unsigned long Random_number_generator::range() const
{
  return gsl_rng_max(generator)-gsl_rng_min(generator);
}

/******************************************************************************/

template <typename ranT>
class Random_distribution_base {
public:
  Random_distribution_base();
  virtual ranT operator()()=0;

  unsigned long raw();
  unsigned long min() const;
  unsigned long max() const;
  unsigned long range() const;

protected:
  gsl_rng                 *generator;
  Random_number_generator *glsim_generator;

} ;

typedef Random_distribution_base<double> rdbase_double;
typedef Random_distribution_base<unsigned long> rdbase_ulong;


template <typename ranT> inline Random_distribution_base<ranT>::
Random_distribution_base()
{
  if (Random_number_generator::generator==0) {
    std::cerr << "Random number generator not initialized\n";
    exit(22);
  }
  generator=Random_number_generator::generator;
  glsim_generator=Random_number_generator::glsim_generator;
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::raw()
{
  return glsim_generator->raw();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::min() const {
  return glsim_generator->min();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::max() const {
  return glsim_generator->max();
}

template <typename ranT>
unsigned long Random_distribution_base<ranT>::range() const {
  return glsim_generator->range();
}

/*****************************************************************************/

class Uniform_integer : public rdbase_ulong {
public:
  Uniform_integer(unsigned long default_m=10);
  unsigned long operator()();
  unsigned long operator()(unsigned long m);

private:
  unsigned long default_m;
} ;

inline
Uniform_integer::Uniform_integer(unsigned long m) :
  default_m(m)
{}

inline unsigned long Uniform_integer::operator()()
{
  return gsl_rng_uniform_int(generator,default_m);
}

inline unsigned long Uniform_integer::operator()(unsigned long m)
{
  return gsl_rng_uniform_int(generator,m);
}

/*****************************************************************************
 *
 * Uniform distribution
 *
 */

class Uniform_real : public rdbase_double {
public:
  Uniform_real(double a=0,double b=1);
  double operator()();

private:
  double a,range;
} ;

inline Uniform_real::Uniform_real(double a_,double b) :
  a(a_)
{
  range=b-a;
}

inline double Uniform_real::operator()()
{
  return a+range*gsl_rng_uniform(generator);
}

/*****************************************************************************
 *
 * Exponential distribution of mean mu
 *
 */

class Exponential_distribution : public rdbase_double {
public:
  Exponential_distribution(double mu_=1) :
    mu(mu_) {}
  double operator()();
  double operator()(double mu);

private:
  double mu;
} ;

inline double Exponential_distribution::operator()()
{
  return gsl_ran_exponential(generator,mu);
}

inline double Exponential_distribution::operator()(double mu_)
{
  return gsl_ran_exponential(generator,mu_);
}

/*****************************************************************************/

// Vectors on the unit sphere

class Spherical3d_distribution : public rdbase_double {
public:
  Spherical3d_distribution() {}

  double operator()();
  void operator()(double *r);
} ;


inline double Spherical3d_distribution::operator()()
{
  std::cerr << "operator() requires 3 arguments in Spherical3d_distribution\n";
  return -1;
}

inline void Spherical3d_distribution::operator()(double *r)
{
  gsl_ran_dir_3d(generator,r,r+1,r+2);
}


#endif /* QDRANDOM_HH */
