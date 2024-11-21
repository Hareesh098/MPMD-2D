/*-----------------------------------------------------------------------------

    MPMD-v2.0 : MULTI POTENTIAL MOLECULAR DYNAMICS 
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
    
    See the README file in the top-level MPMD directory.

-----------------------------------------------------------------------------*/



#include<stdio.h>
#include"global.h"
double Integrate(double *, int);
void PrintVacf();
void ZeroVacf();
void AccumVacf(){
  double fac;
  int j, nb;
  for(nb = 1 ; nb <= nBuffAcf ; nb ++){
    if(indexAcf[nb] == nValAcf){
      for(j = 1 ; j <= nValAcf; j ++){
	viscAcfAv[j] +=  viscAcf[nb][j];
      }
      indexAcf[nb] = 0;
      countAcfAv ++;
      if(countAcfAv == limitAcfAv){
	fac = 1./(kinEnergy*region[1]*region[2]*limitAcfAv);
	viscAcfInt = fac*stepAcf*deltaT*Integrate(viscAcfAv, nValAcf);
	PrintVacf();
	ZeroVacf();
      }
    }
  }
}

