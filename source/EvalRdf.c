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
void EvalRdf(){    
  real dr[NDIM+1], deltaR, normFac, rr, rrRange, theta, thetaRange, deltaTheta;
  int j1, j2, n, m;
  countRdf ++;
  if(countRdf == 1){
    for(n = 1 ; n <= sizeHistRdf ; n ++){
      for(m = 1 ; m <= sizeHistRdf ; m ++){
        histRdf[n][m] = 0.;
        HistRdf[n] = 0.;
      }
    }
  }
  rrRange = Sqr(rangeRdf);
  deltaR = rangeRdf / sizeHistRdf;
  thetaRange = 2.*M_PI;
  deltaTheta = thetaRange / sizeHistRdf;
  for(j1 = 1 ; j1 <= nAtom  ; j1 ++){
    for(j2 = 1 ; j2 <= nAtom ; j2 ++){

      dr[1] = rx[j1] - rx[j2];
      if(fabs(dr[1]) > regionH[1])
	dr[1] -= SignR(region[1], dr[1]);

      dr[2] = ry[j1] - ry[j2];
      if(fabs(dr[2]) > regionH[2])
	dr[2] -= SignR(region[2], dr[2]);

      rr = Sqr(dr[1]) + Sqr(dr[2]);
     theta = atan(dr[2]/dr[1]);
     if (dr[1] == 0. && dr[2] == 0.) theta = 0.;
     if( dr[1] < 0. && dr[2] < 0.) theta = theta;
     else if( dr[1] > 0. && dr[2] < 0.) theta += M_PI;
     else if( dr[1] > 0. && dr[2] > 0.) theta += M_PI;
     else if( dr[1] < 0. && dr[2] > 0.) theta += 2.*M_PI;
     else if( dr[1] < 0. && dr[2] == 0.) theta = 0.;
     else if( dr[1] == 0.&& dr[2] < 0.) theta = 0.5*M_PI;
     else if( dr[1] > 0. && dr[2] == 0.) theta = M_PI;
     else if( dr[1] == 0.&& dr[2] > 0.) theta = 1.5*M_PI;	
     else if( dr[1] < 0. && dr[2] == 0.) theta = 2.*M_PI; 


     if(rr < rrRange){
	n = (int)(sqrt(rr)/deltaR) + 1;
	m = (int)(theta/deltaTheta) + 1; 
	histRdf[n][m] ++;
      }
    }
  }

//For angular vales of Rdf
  if(countRdf == limitRdf){
    normFac = region[1]*region[2] / (M_PI*Sqr(deltaR)*nAtom*nAtom*countRdf );
    for(n = 1 ; n <= sizeHistRdf ; n ++){
       for(m = 1 ; m <= sizeHistRdf ; m ++){
      histRdf[n][m] *= normFac/(n-0.5);
      }}
//Print the Rdf data
    int n;
     histRdf[1][1] = 0.;
    fprintf(fprdf,"rdf @ timeNow %lf\n", timeNow);
    for(n = 1 ; n <= sizeHistRdf ; n ++){
      rBin = (n - 0.5)*rangeRdf/sizeHistRdf;
      //fprintf(fprdf, "%lf\t", rBin);
      for(m = 1 ; m <= sizeHistRdf ; m ++){
      fprintf(fprdf, "%lf\t",histRdf[n][m]);
       }
     fprintf(fprdf,"\n");
   }
//Summing all the theta values
    
    for(n = 1 ; n <= sizeHistRdf ; n ++){
     for(m = 1 ; m <= sizeHistRdf ; m ++){
       HistRdf[n] += histRdf[n][m];
     }
   }
   for(n = 1 ; n <= sizeHistRdf ; n ++) {
     rBin = (n - 0.5)*rangeRdf/sizeHistRdf;
     fprintf(fprdfAvg, "%lf\t%lf\n",rBin,HistRdf[n]);
    }
  }
}

