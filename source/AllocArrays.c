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
void AllocArrays(){
  int n;

  // Spacetime correlation for k = k_x
    spacetimeFluct = (real **) malloc ((nBuffCorr + 1)*sizeof(real *));
    complxFuncOrg = (real **) malloc ((nBuffCorr + 1)*sizeof(real *));
    for (n = 0; n <= nBuffCorr; n ++){
     spacetimeFluct[n] = (real *) malloc ((3*nFunCorr * nValCorr + 1)*sizeof(real ));
     complxFuncOrg[n] = (real *) malloc ((6 * nFunCorr + 1)*sizeof(real )); 
     }
     complxFunc = (real *) malloc ((6 * nFunCorr + 1)*sizeof(real ));
     spacetimeFluctAv = (real *) malloc ((3*nFunCorr * nValCorr + 1)*sizeof(real ));
     indexFluct = (int *) malloc ((nFunCorr * nValCorr + 1)*sizeof(int ));

  // Spacetime correlation for k = k_y
    spacetimeFluct2 = (real **) malloc ((nBuffCorr + 1)*sizeof(real *));
    complxFuncOrg2 = (real **) malloc ((nBuffCorr + 1)*sizeof(real *));
    for (n = 0; n <= nBuffCorr; n ++){
     spacetimeFluct2[n] = (real *) malloc ((3*nFunCorr2 * nValCorr + 1)*sizeof(real ));
     complxFuncOrg2[n] = (real *) malloc ((6 * nFunCorr2 + 1)*sizeof(real )); 
     }
     complxFunc2 = (real *) malloc ((6 * nFunCorr2 + 1)*sizeof(real ));
     spacetimeFluctAv2 = (real *) malloc ((3*nFunCorr2 * nValCorr + 1)*sizeof(real ));
     indexFluct2 = (int *) malloc ((nFunCorr2 * nValCorr + 1)*sizeof(int ));

  // VISCOSITY
    indexAcf = (int *)malloc((nBuffAcf+1)*sizeof(int));
    viscAcf = (real **)malloc((nBuffAcf+1)*sizeof(real *));
    for(n = 0 ; n <= nBuffAcf ; n ++)
      viscAcf[n] = (real *)malloc((nValAcf+1)*sizeof(real ));
      viscAcfOrg = (real *)malloc((nBuffAcf+1)*sizeof(real));
      viscAcfAv = (real *)malloc((nValAcf+1)*sizeof(real));

  // RDF
    histRdf = (real **)malloc((sizeHistRdf+1)*sizeof(real*));
    HistRdf = (real *)malloc((sizeHistRdf+1)*sizeof(real));
    for(n = 0; n <= sizeHistRdf; n ++) 
      histRdf[n] = (real *)malloc((sizeHistRdf+1)*sizeof(real));
 
  // Diffusion 
    rDiffuseTrue = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
    rDiffuseOrg = (real **) malloc ((nBuffAcf + 1)*sizeof(real *)); 
    xDiffuse = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
    yDiffuse = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
    rDiffuse = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
   for(n = 0; n <= nBuffAcf; n ++){
     rDiffuseTrue[n] = (real *) malloc ((NDIM*nAtom + 1)*sizeof(real));
     rDiffuseOrg[n] = (real *) malloc ((NDIM*nAtom + 1)*sizeof(real));
     xDiffuse[n] = (real *) malloc ((nValAcf + 1)*sizeof(real));
     yDiffuse[n] = (real *) malloc ((nValAcf + 1)*sizeof(real));
     rDiffuse[n] = (real *) malloc ((nValAcf + 1)*sizeof(real));
    }
     xDiffuseAv = (real *) malloc((nValAcf + 1)*sizeof(real)) ;
     yDiffuseAv = (real *) malloc((nValAcf + 1)*sizeof(real)) ;
     rDiffuseAv = (real *) malloc((nValAcf + 1)*sizeof(real)) ;
     MSDAv = (real *) malloc((nValAcf + 1)*sizeof(real)) ;
     indexDiffuse = (int *) malloc((nBuffAcf +1)*sizeof(int)); 
 
  // Velocity autocorrelation function
     VeloAcfOrg = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
     VeloAcf = (real **) malloc ((nBuffAcf + 1)*sizeof(real *));
     for(n = 0; n <= nBuffAcf; n ++){
       VeloAcfOrg[n] = (real *) malloc ((NDIM*nAtom + 1)*sizeof(real ));
       VeloAcf[n] = (real *) malloc ((nValAcf + 1)*sizeof(real ));
    }
      indexVeloAcf = (int *) malloc((nBuffAcf +1)*sizeof(int));
      VeloAcfAv = (real *) malloc((nValAcf + 1)*sizeof(real)) ;

  // Velocity distribution
     histVel = (real *) malloc((sizeHistVel + 1)*sizeof(real));

  // Configurational Temperature 
     A = (double *)malloc((nAtom+1)*sizeof(double));
     Ax = (double *)malloc((nAtom+1)*sizeof(double));
     Ay = (double *)malloc((nAtom+1)*sizeof(double));
     TValSum = (double *)malloc((nAtom+1)*sizeof(double)); 
}






