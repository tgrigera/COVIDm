/*
 * gillespieav.hh
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

#ifndef GILLESPIEAV_HH
#define GILLESPIEAV_HH

template <typename Data>
class Gillespie_av {
public:
  Gillespie_av(double t0,double deltat);
  void push(double t,Data& d);
  
private:
  double t0,deltat,tmax;
  double tlast;
  Data   datalast;

  void ipush(double t,Data d) {std::cerr << "**** pushing time " << t << '\n';}

} ;

template <typename Data>
Gillespie_av<Data>::Gillespie_av(double t0,double deltat,double tmax) :
  t0(t0),
  deltat(deltat),
  tmax(tmax),
  tlast(t0-10*deltat)
{
}

template <typename Data>
void Gillespie_av<Data>::push(double time,Data &data)
{
  std::cerr << "received time " << time << "\n";
  if (time==t0) {
    datalast=data;
    tlast=time;
    ipush(tlast,datalast);
  } else if (time<tlast+deltat) {
    datalast=data;
  } else {
    double t;
    for (t=tlast+deltat; t<time && t<=tmax; t+=deltat)
      ipush(t,datalast);
    datalast=data;
    if (time==t) {
      ipush(t,data); tlast=t;}
    else tlast=t-deltat;
  }
}

#endif /* GILLESPIEAV_HH */
