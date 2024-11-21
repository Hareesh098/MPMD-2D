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
#include<stdlib.h>
#include"global.h"
void Close(){
  free(rx);
  free(ry);
  free(vx);
  free(vy);
  free(ax);
  free(ay);
  free(fax);
  free(fay);
  free(cellList);
  free(atomType);

  int n;
//Spacetime correlation
  free(complxFunc2);
  free(indexFluct2);
  free(spacetimeFluctAv2);
  for (n = 0; n <= nBuffCorr; n ++ ) {
   free(spacetimeFluct2[n]);
   free(complxFuncOrg2[n]);
  }
  free(spacetimeFluct2);
  free(complxFuncOrg2);

//Viscocity 
  free(indexAcf);
  free(viscAcfOrg);
  free(viscAcfAv);
  for(n = 0 ; n <= nBuffAcf ; n ++)
    free(viscAcf[n]);
  free(viscAcf);

//Diffusion
  free(xDiffuseAv);
  free(yDiffuseAv);
  free(rDiffuseAv);
 for(n = 0; n <= nBuffAcf; n ++){
   free(rDiffuseTrue[n]);
   free(rDiffuseOrg[n]);
   free(xDiffuse[n]);
   free(yDiffuse[n]);
   free(rDiffuse[n]);
 }  
 free(rDiffuseTrue);
 free(rDiffuseOrg);
 free(xDiffuse);
 free(yDiffuse);
 free(rDiffuse);
 free(indexDiffuse);


//Velocity auto correlation function
 free(VeloAcfAv);
 free(indexVeloAcf);
 for(n = 0; n <= nBuffAcf; n ++){
	free(VeloAcfOrg[n]);
  free(VeloAcf[n]);
  }
  free(VeloAcfOrg);
  free(VeloAcf);


//Spacetime correlation
  free(complxFunc);
  free(indexFluct);
  free(spacetimeFluctAv);
  for (n = 0; n <= nBuffCorr; n ++ ) {
   free(spacetimeFluct[n]);
   free(complxFuncOrg[n]);
  }
  free(spacetimeFluct);
  free(complxFuncOrg);	

//Velocity distribution
  free(histVel);

//Configurational temperature
  free(atomType);
  free(A);
  free(Ax);
  free(Ay);
  free(TValSum);
}



