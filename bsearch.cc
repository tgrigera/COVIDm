/*
 * bsearch.cc - binary search
 *
 * Binary search function for faster lookup of cumulative rate table
 * (for Gillespie algorithm).  Overloaded to accept C++ vectors or C
 * arrays but otherwise identical.
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

#include <vector>

/*
 * int bsearch(double r,int N,double crates[])
 *
 * input: 
 *  - crates[0..N-1] : cumulative rates (sorted ascending)
 *  - r              : number between 0 and crates[N-1]
 *
 * return:
 *  - index n such that  crates[n] < r <= crates[n+1]
 *
 */

int bsearch(double r,int N,double crates[])
{
  int i=0;
  int j=N-1;

  while (j!=i+1) {
    int k= i + (j-i)/2;
    if ( crates[k] < r ) i=k;
    else j=k;
  }
  return i;
}

/*
 * int bsearch(double rstd::vector<double> &crates)
 *
 * input: 
 *  - crates : cumulative rates (sorted ascending)
 *  - r      : number between 0 and crates.back()
 *
 * return:
 *  - index n such that  crates[n] < r <= crates[n+1]
 *
 */

int bsearch(double r,std::vector<double> &crates)
{
  int i=0;
  int j=crates.size()-1;

  while (j!=i+1) {
    int k= i + (j-i)/2;
    if ( crates[k] < r ) i=k;
    else j=k;
  }
  return i;
}

