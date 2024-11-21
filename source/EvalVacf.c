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
void AccumVacf();
void EvalVacf(){
  int n, nb, ni;
  double viscVec = 0.;
  double v[3];
  for(n = 1 ; n <= nAtom ; n ++){
    v[1] = vx[n] - 0.5*ax[n]*deltaT;
    v[2] = vy[n] - 0.5*ay[n]*deltaT;
    viscVec += v[1]*v[2];
  }
  viscVec += rfAtom;
  for(nb = 1 ; nb <= nBuffAcf ; nb ++){
    indexAcf[nb] ++;
    if(indexAcf[nb] <= 0)continue;
    if(indexAcf[nb] == 1){
      viscAcfOrg[nb] = viscVec;
    }
    ni = indexAcf[nb];
    viscAcf[nb][ni] = viscAcfOrg[nb]*viscVec;
  }
  AccumVacf();
}
