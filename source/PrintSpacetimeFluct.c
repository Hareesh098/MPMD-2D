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
void PrintSpacetimeFluct (){
  int j, n, nn;
  real tVal = 0.;
  fprintf(fplongfluctKx, "longitudinal current fluctuation for k = (kx,0)\n");
  fprintf(fplongfluctKx,"NDIM %d\n", NDIM);
  fprintf(fplongfluctKx,"nAtom %d\n", nAtom);
  fprintf(fplongfluctKx,"region %lf\n", region[1]);
  fprintf(fplongfluctKx,"nFunCorr %d\n", nFunCorr);
  fprintf(fplongfluctKx,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fplongfluctKx,"stepCorr %d\n", stepCorr);
  fprintf(fplongfluctKx,"nValCorr %d\n", nValCorr);
  fprintf(fplongfluctKx,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fplongfluctKx, "%e\t", tVal);
   nn = 3*nFunCorr * (n - 1);
   for (j = 1; j <= 3*nFunCorr; j +=3)
    fprintf(fplongfluctKx, "%e\t", spacetimeFluctAv[nn + j]);
    fprintf(fplongfluctKx, "\n");
  }


  fprintf(fptransfluctKx, "transverse current fluctuation for k = (kx,0)\n");
  fprintf(fptransfluctKx,"NDIM %d\n", NDIM);
  fprintf(fptransfluctKx,"nAtom %d\n", nAtom);
  fprintf(fptransfluctKx,"region %lf\n", region[1]);
  fprintf(fptransfluctKx,"nFunCorr %d\n", nFunCorr);
  fprintf(fptransfluctKx,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fptransfluctKx,"stepCorr %d\n", stepCorr);
  fprintf(fptransfluctKx,"nValCorr %d\n", nValCorr);
  fprintf(fptransfluctKx,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fptransfluctKx, "%e\t", tVal);
   nn = 3*nFunCorr * (n - 1);
   for (j = 2; j <= 3*nFunCorr; j += 3)
    fprintf(fptransfluctKx, "%e\t", spacetimeFluctAv[nn + j]);
   fprintf(fptransfluctKx, "\n");
  }

  fprintf(fpdnstyfluctKx, "density fluctuation for k = (kx,0)\n");
  fprintf(fpdnstyfluctKx,"NDIM %d\n", NDIM);
  fprintf(fpdnstyfluctKx,"nAtom %d\n", nAtom);
  fprintf(fpdnstyfluctKx,"region %lf\n", region[1]);
  fprintf(fpdnstyfluctKx,"nFunCorr %d\n", nFunCorr);
  fprintf(fpdnstyfluctKx,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fpdnstyfluctKx,"stepCorr %d\n", stepCorr);
  fprintf(fpdnstyfluctKx,"nValCorr %d\n", nValCorr);
  fprintf(fpdnstyfluctKx,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fpdnstyfluctKx, "%e\t", tVal); 
   nn = 3*nFunCorr * (n - 1);
   for (j = 3; j <= 3*nFunCorr; j += 3)
    fprintf(fpdnstyfluctKx, "%e\t", spacetimeFluctAv[nn + j]);
    fprintf(fpdnstyfluctKx, "\n");
  }



}
