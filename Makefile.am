#
# Makefila.am -- automake script for COVIDm
#
# COVIDm is copyright (c) 2020 by the authors (see AUTHORS)
#
# Process this file with autoconf to produce a configure script.

bin_PROGRAMS = sir sir_m sir_f

sir_SOURCES = sir.cc qdrandom.cc

sir_f_SOURCES = sir_f.cc bsearch.cc popstate.cc geoave.cc qdrandom.cc

sir_m_SOURCES = sir_m.cc popstate.cc geoave.cc qdrandom.cc

noinst_HEADERS = bsearch.hh qdrandom.hh read_arg.hh popstate.hh geoave.hh
