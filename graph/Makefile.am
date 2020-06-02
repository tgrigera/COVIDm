#
# Makefile.am -- automake script for COVIDm (graph subdirectory)
#
# COVIDm is copyright (c) 2020 by the authors (see AUTHORS)
#

bin_PROGRAMS = sir_fc sir_sq

sir_sq_SOURCES = sir_sq.cc emodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc

sir_fc_SOURCES = sir_sq.cc emodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc
sir_fc_CPPFLAGS = -DDEBUG_FCGRAPH

EXTRA_DIST = emodel.hh esampler.hh egraph.hh eevents.hh seir_collector.hh

