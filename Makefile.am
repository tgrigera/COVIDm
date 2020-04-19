#
# Makefile.am -- automake script for COVIDm
#
# COVIDm is copyright (c) 2020 by the authors (see AUTHORS)
#


bin_PROGRAMS = sir sir_m sir_f seeiir_i1 seeiir_i2 seeiir_i3

sir_SOURCES = sir.cc qdrandom.cc

sir_f_SOURCES = sir_f.cc bsearch.cc popstate.cc geoave.cc qdrandom.cc

sir_m_SOURCES = sir_m.cc popstate.cc geoave.cc qdrandom.cc

seeiir_i1_SOURCES = seeiir_main.cc bsearch.cc popstate.cc geoave.cc qdrandom.cc
seeiir_i1_CPPFLAGS = -DSEEIIR_IMPLEMENTATION_1

seeiir_i2_SOURCES = seeiir_main.cc bsearch.cc popstate.cc geoave.cc qdrandom.cc
seeiir_i2_CPPFLAGS = -DSEEIIR_IMPLEMENTATION_2

seeiir_i3_SOURCES = seeiir_main.cc bsearch.cc popstate.cc geoave.cc qdrandom.cc
seeiir_i3_CPPFLAGS = -DSEEIIR_IMPLEMENTATION_3

noinst_HEADERS = bsearch.hh qdrandom.hh read_arg.hh popstate.hh geoave.hh
