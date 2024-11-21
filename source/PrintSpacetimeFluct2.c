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
void PrintSpacetimeFluct2 (){
  int j, n, nn;
  real tVal = 0.;
  fprintf(fplongfluctKy, "longitudinal current fluctuation for k = (0,ky)\n");
  fprintf(fplongfluctKy,"NDIM %d\n", NDIM);
  fprintf(fplongfluctKy,"nAtom %d\n", nAtom);
  fprintf(fplongfluctKy,"region[2] %lf\n", region[2]);
  fprintf(fplongfluctKy,"nFunCorr2 %d\n", nFunCorr2);
  fprintf(fplongfluctKy,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fplongfluctKy,"stepCorr %d\n", stepCorr);
  fprintf(fplongfluctKy,"nValCorr %d\n", nValCorr);
  fprintf(fplongfluctKy,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fplongfluctKy, "%e\t", tVal);
   nn = 3*nFunCorr2 * (n - 1);
   for (j = 1; j <= 3*nFunCorr2; j +=3)
    fprintf(fplongfluctKy, "%e\t", spacetimeFluctAv2[nn + j]);
    fprintf(fplongfluctKy, "\n");
  }


  fprintf(fptransfluctKy, "transverse current fluctuation for k = (0,ky)\n");
  fprintf(fptransfluctKy,"NDIM %d\n", NDIM);
  fprintf(fptransfluctKy,"nAtom %d\n", nAtom);
  fprintf(fptransfluctKy,"region[2] %lf\n", region[2]);
  fprintf(fptransfluctKy,"nFunCorr2 %d\n", nFunCorr2);
  fprintf(fptransfluctKy,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fptransfluctKy,"stepCorr %d\n", stepCorr);
  fprintf(fptransfluctKy,"nValCorr %d\n", nValCorr);
  fprintf(fptransfluctKy,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fptransfluctKy, "%e\t", tVal);
   nn = 3*nFunCorr2 * (n - 1);
   for (j = 2; j <= 3*nFunCorr2; j += 3)
    fprintf(fptransfluctKy, "%e\t", spacetimeFluctAv2[nn + j]);
   fprintf(fptransfluctKy, "\n");
  }

  fprintf(fpdnstyfluctKy, "density fluctuation for k = (0,ky)\n");
  fprintf(fpdnstyfluctKy,"NDIM %d\n", NDIM);
  fprintf(fpdnstyfluctKy,"nAtom %d\n", nAtom);
  fprintf(fpdnstyfluctKy,"region[2] %lf\n", region[2]);
  fprintf(fpdnstyfluctKy,"nFunCorr2 %d\n", nFunCorr2);
  fprintf(fpdnstyfluctKy,"limitCorrAv %d\n", limitCorrAv);
  fprintf(fpdnstyfluctKy,"stepCorr %d\n", stepCorr);
  fprintf(fpdnstyfluctKy,"nValCorr %d\n", nValCorr);
  fprintf(fpdnstyfluctKy,"deltaT %lf\n", deltaT); 	
  for (n = 1; n <= nValCorr; n ++){
   tVal = (n - 1) * stepCorr * deltaT;
   fprintf(fpdnstyfluctKy, "%e\t", tVal); 
   nn = 3*nFunCorr2 * (n - 1);
   for (j = 3; j <= 3*nFunCorr2; j += 3)
    fprintf(fpdnstyfluctKy, "%e\t", spacetimeFluctAv2[nn + j]);
    fprintf(fpdnstyfluctKy, "\n");
  }



}
