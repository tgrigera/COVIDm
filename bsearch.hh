/*
 * bsearch.hh - header for binary search
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

#ifndef BSEARCH_HH
#define BSEARCH_HH

int bsearch(double r,int N,double crates[]);

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

template <typename T>
size_t bsearch(T r,std::vector<T> &crates)
{
  size_t i=0;
  size_t j=crates.size();

  while (i<j) {
    int k= (i+j)/2;
    if ( crates[k] < r ) i=k+1;
    else j=k;
  }
  return i-1;
}

#endif /* BSEARCH_HH */
