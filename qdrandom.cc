/*
 * qdrandom.cc -- quick and dirty interface to GSL's random numbers
 *
 * This is taken from glsim's interface, with less features (no
 * exceptions, no multiple scopes, no serialization through Boost).
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * Coronel Pringles, 24-jul-2017
 *
 */

#include "qdrandom.hh"

Random_number_generator::Random_number_generator(const unsigned long seed)
{
  generator=gsl_rng_alloc(::gsl_rng_mt19937);
  glsim_generator=this;
  set_seed(seed);
}

Random_number_generator::~Random_number_generator()
{
  gsl_rng_free(generator);
  generator=0;
  glsim_generator=0;
}

void Random_number_generator::save(std::ostream& os)
{
  void *state=gsl_rng_state(generator);
  os.write((char*) state,gsl_rng_size(generator));
}

void Random_number_generator::load(std::istream& is)
{
  void *state=gsl_rng_state(generator);
  is.read((char*) state,gsl_rng_size(generator));
}

/*
 * Static data
 */

gsl_rng* Random_number_generator::generator=0;
Random_number_generator* Random_number_generator::glsim_generator=0;
