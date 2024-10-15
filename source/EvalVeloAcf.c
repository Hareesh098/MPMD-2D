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
  void AccumVeloAcf () ;
  void EvalVeloAcf(){
  int k, n, nb, ni;
  for(nb =1; nb <= nBuffAcf; nb ++ ){
  indexVeloAcf[nb] ++ ;
  if (indexVeloAcf[nb] <= 0 ) continue ;
  if (indexVeloAcf[nb] == 1 ) {
  k = 1 ;
  for(n =1; n <= nAtom; n++){
  VeloAcfOrg[nb][k] = vx[n] - 0.5 * deltaT * ax[n];
  VeloAcfOrg[nb][k+1] = vy[n] - 0.5 * deltaT * ay[n];
  k += 2 ;
  } }
  ni = indexVeloAcf[nb] ; VeloAcf[nb][ni] = 0. ;
  k = 1 ;
  for (n = 1; n <= nAtom; n++){
  VeloAcf[nb][ni] += ( VeloAcfOrg[nb][k] * (vx[n] - 0.5 * deltaT *ax[n]) + VeloAcfOrg[nb][k+1] * (vy[n] - 0.5 * deltaT *ay[n]))  ; 
  k += 2;
  } }
  AccumVeloAcf () ;		
  }
