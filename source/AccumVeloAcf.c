/*-----------------------------------------------------------------------------

    MPMD-v2.0 : MULTI POTENTIAL MOLECULAR DYNAMICS-version 2.0 
    A parallel classical molecular dynamics code
    Copyright (C) 2018  Harish Charan, charan.harish@gmail.com

    This program is free software but a proper permission must be taken: you can 
    redistribute it and/or modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
    more details.

    You should have received a copy of the GNU General Public License along 
    with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    See the README file in the top-level MPMD-v2.0 directory.

-----------------------------------------------------------------------------*/



#include<stdio.h>
#include"global.h"
double Integrate(double *, int);
void PrintVeloAcf () ;
void ZeroVeloAcf  () ;
void AccumVeloAcf(){
  int j, nb;
  real fac;
  for(nb = 1; nb <= nBuffAcf; nb ++){
  if (indexVeloAcf[nb] == nValAcf)	{
  for (j = 1; j <= nValAcf; j ++)
  VeloAcfAv[j] += VeloAcf[nb][j];
  indexVeloAcf[nb] = 0 ; countVeloAcfAv ++ ;
  if (countVeloAcfAv == limitAcfAv){
  fac = 1. / (NDIM * nAtom * limitAcfAv);
  diffuseVeloAcfInt = fac * stepAcf * deltaT * Integrate (VeloAcfAv, nValAcf);
  //for(j = 2; j <= nValAcf; j ++)
  //VeloAcfAv[j] = VeloAcfAv[j] / VeloAcfAv[1] ;
  //VeloAcfAv[1] = 1. ;
  PrintVeloAcf () ;
  ZeroVeloAcf  () ;

  }	}  }
}
