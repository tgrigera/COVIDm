#
# Makefile.am -- automake script for COVIDm (graph subdirectory)
#
# COVIDm is copyright (c) 2020 by the authors (see AUTHORS)
#

bin_PROGRAMS = sir_fc sir_sq seeiir_fc seeiir_fc_altR seeiir_sq

sir_sq_SOURCES = sir_sq.cc emodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc

sir_fc_SOURCES = sir_sq.cc emodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc
sir_fc_CPPFLAGS = -DDEBUG_FCGRAPH

seeiir_fc_SOURCES = seeiir_fc.cc emodel.cc seirmodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc

seeiir_fc_altR_SOURCES = seeiir_fc.cc emodel.cc seirmodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc
seeiir_fc_altR_CPPFLAGS = -DALT_R0

seeiir_sq_SOURCES = seeiir_sq.cc emodel.cc seirmodel.cc seir_collector.cc egraph.cc ../qdrandom.cc ../geoave.cc

EXTRA_DIST = emodel.hh esampler.hh egraph.hh eevents.hh seir_collector.hh

