/*
 * read_arg.hh - simple template functions to aid in argument parsing
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

#include <stdexcept>

template <typename T>
void read_arg(char *argv[],T& opt)
{ throw std::runtime_error("unimplemented argument conversion"); }

template <> void read_arg<double>(char *argv[],double& opt)
{ options.last_arg_read++; opt = atof(argv[options.last_arg_read]); }

template <> void read_arg<int>(char *argv[],int& opt)
{ options.last_arg_read++; opt = atoi(argv[options.last_arg_read]); }

template <> void read_arg<long>(char *argv[],long& opt)
{ options.last_arg_read++; opt = atol(argv[options.last_arg_read]); }

template <> void read_arg<char*>(char *argv[],char*& opt)
{ options.last_arg_read++; opt = argv[options.last_arg_read]; }

