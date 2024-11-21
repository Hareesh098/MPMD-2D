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
#include<math.h>
#include"global.h"
void AccumSpacetimeFluct2();
void EvalSpacetimeFluct2(){
  real cosV=0., cosV0=0., cosV1=0., cosV2=0., sinV=0., sinV1=0., sinV2=0., w;
  
  int j, m, n, nb, ni, nv, nc;
  real kMin = 2. * M_PI / region[1];
  
  
  
  for (j = 1; j <= 6*nFunCorr2; j++)
    complxFunc2[j] = 0.;

  for (n = 1; n <= nAtom; n++){
    j = 1;
    for (m = 1; m <= nFunCorr2; m++){
      if(m == 1){
	cosV = cos(kMin*rx[n]);
	sinV = sin(kMin*rx[n]);
	cosV0 = cosV;
      }else if(m == 2){
	cosV1 = cosV;
	sinV1 = sinV;
	cosV = 2.*cosV0*cosV1-1;
	sinV = 2.*cosV0*sinV1;
      }else{
	cosV2 = cosV1;
	sinV2 = sinV1;
	cosV1 = cosV;
	sinV1 = sinV;
	cosV = 2.*cosV0*cosV1-cosV2;
	sinV = 2.*cosV0*sinV1-sinV2;
      }
      for(nc = 1; nc <= 3; nc++){
        if(nc == 1) w = vx[n] - 0.5 * deltaT * ax[n] ;
         
        else if(nc == 2) w = vy[n] - 0.5 * deltaT * ay[n] ;
         
        else w = 1.; 
          
       complxFunc2[j] += w*cosV;
       complxFunc2[j+1] += w*sinV;
      j += 2;
     }
    }
  }

  for (nb = 1; nb <= nBuffCorr; nb++){
    indexFluct2[nb] += 1;
    if (indexFluct2[nb] <= 0) continue;
    ni = 3*nFunCorr2 * (indexFluct2[nb] - 1);
    if (indexFluct2[nb] == 1){
      for (j = 1; j <= 6*nFunCorr2; j++)
       complxFuncOrg2[nb][j] = complxFunc2[j];
    }

    for (j = 1; j <= 3*nFunCorr2; j++)
      spacetimeFluct2[nb][ni + j] = 0.;
    
    j = 1;
    for (m = 1; m <= nFunCorr2; m++){
      w = Sqr(kMin * m);
     for(nc = 1; nc <= 3; nc ++){
       nv = 3*m + ni;
       if(nc == 1) nv -= 2; 
 
       else if(nc == 2) nv -= 1;
         
       else{
         nv = nv; w = 1.;
        }

       spacetimeFluct2[nb][nv] += w * ( complxFunc2[j] * complxFuncOrg2[nb][j] + complxFunc2[j + 1] * complxFuncOrg2[nb][j + 1] ) ;
      j += 2;
      }
    }
  }
AccumSpacetimeFluct2();
}
