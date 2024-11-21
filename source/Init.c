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
#include<stdlib.h>
#include"global.h"
void Init(){
  char dummy[128];
  FILE *fp;
  fp = fopen("input-data","r");
  fscanf(fp, "%s %lf", dummy, &GAMMA);
  fscanf(fp, "%s %lf", dummy, &rCut);
  fscanf(fp, "%s %lf", dummy, &kappa);
  fscanf(fp, "%s %lf", dummy, &deltaT);
  fscanf(fp, "%s %d", dummy, &stepAvg);
  fscanf(fp, "%s %d", dummy, &stepEquil);
  fscanf(fp, "%s %d", dummy, &stepLimit);
  fscanf(fp, "%s %d", dummy, &stepDump);
  fscanf(fp, "%s %d", dummy, &stepTrajectory);
  fscanf(fp, "%s %c", dummy, &thermo);
  fscanf(fp, "%s %c", dummy, &BC);
  fscanf(fp, "%s %d", dummy, &limitCorrAv);
  fscanf(fp, "%s %d", dummy, &nBuffCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr2);
  fscanf(fp, "%s %d", dummy, &nValCorr);
  fscanf(fp, "%s %d", dummy, &stepCorr);
  fscanf(fp, "%s %d", dummy, &limitAcfAv);
  fscanf(fp, "%s %d", dummy, &nBuffAcf);
  fscanf(fp, "%s %d", dummy, &nValAcf);
  fscanf(fp, "%s %d", dummy, &stepAcf);
  fscanf(fp, "%s %lf", dummy, &rangeRdf);
  fscanf(fp, "%s %d", dummy, &limitRdf);
  fscanf(fp, "%s %d", dummy, &sizeHistRdf);
  fscanf(fp, "%s %d", dummy, &stepRdf);
  fscanf(fp, "%s %d", dummy, &limitVel);
  fscanf(fp, "%s %lf", dummy, &rangeVel);
  fscanf(fp, "%s %d", dummy, &sizeHistVel);
  fscanf(fp, "%s %d", dummy, &stepVel);
  

  fclose(fp);

  FILE *fpSTATE;
  if((fpSTATE = fopen("../STATE","r"))==NULL)
    fprintf(fpresult,"Could not open ../STATE file\n");
  fscanf(fpSTATE, "%s %lf", dummy, &timeNow);
  fscanf(fpSTATE, "%s %d", dummy, &nAtom);
  fscanf(fpSTATE, "%s %lf", dummy, &region[1]);
  fscanf(fpSTATE, "%s %lf", dummy, &region[2]);
  density = nAtom/(region[1]*region[2]);
  cells[1] = region[1] / rCut;
  cells[2] = region[2] / rCut;
  cellList = (int *)malloc((nAtom + cells[1] * cells[2] + 1) * sizeof(int));
  regionH[1] = 0.5*region[1];
  regionH[2] = 0.5*region[2];

  rx = (real *)malloc( (nAtom + 1) * sizeof(real));
  ry = (real *)malloc( (nAtom + 1) * sizeof(real));
  vx = (real *)malloc( (nAtom + 1) * sizeof(real));
  vy = (real *)malloc( (nAtom + 1) * sizeof(real));
  ax = (real *)malloc( (nAtom + 1) * sizeof(real));
  ay = (real *)malloc( (nAtom + 1) * sizeof(real));
  fax = (real *)malloc( (nAtom + 1) * sizeof(real));
  fay = (real *)malloc( (nAtom + 1) * sizeof(real));
  atomType = (int *)malloc((nAtom+1)*sizeof(int));

  int n, idx;
  for(n = 1; n <= nAtom; n ++)
    fscanf(fpSTATE, "%d %lf %lf %lf %lf", &idx, &rx[n], &ry[n], &vx[n], &vy[n]);
  fclose(fpSTATE);

  
  if(rank == master){
    fprintf(fpresult, "------------------------------------\n");
    fprintf(fpresult, "------------ PARAMETERS ------------\n");
    fprintf(fpresult, "------------------------------------\n");
    fprintf(fpresult, "nAtom            %d\n", nAtom);
    fprintf(fpresult, "kappa            %lf\n", kappa);
    fprintf(fpresult, "density          %lf\n", density);
    fprintf(fpresult, "rCut             %lf\n", rCut);
    fprintf(fpresult, "deltaT		%lf\n", deltaT);
    fprintf(fpresult, "stepEquil	%d\n", stepEquil);
    fprintf(fpresult, "stepLimit        %d\n", stepLimit);
    fprintf(fpresult, "region[1]        %lf\n", region[1]);
    fprintf(fpresult, "region[2]        %lf\n", region[2]);
    fprintf(fpresult, "cells[1]         %d\n", cells[1]);
    fprintf(fpresult, "cells[2]         %d\n", cells[2]);
    fprintf(fpresult, "------------------------------------\n");
  }
}
