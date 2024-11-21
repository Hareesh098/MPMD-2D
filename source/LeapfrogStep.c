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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include"global.h"


void LeapfrogStep(char thermo, gsl_rng * rnd){
  if(stepCount <= stepEquil){
    double gSum, varS, massS;
    temperature = 1./GAMMA;
    
   if(stepCount == 1) varS = 0.;
    double A, S1, S2, T;
    int n;
    S1 = 0.; S2 = 0.; gSum = 0.; massS = 0.1; 

    vvSum = 0.;
    double halfdt = 0.5*deltaT;
    for (n = 1; n <= nAtom; n++){
      T = vx[n] + halfdt * ax[n];
      S1 += T * ax[n];
      S2 += Sqr(T);
      
      T = vy[n] + halfdt * ay[n];
      S1 += T * ay[n];
      S2 += Sqr(T);
     vvSum += (Sqr(vx[n]) +  Sqr(vy[n]));
    }

    A = -S1 / S2;
    S2 = vvSum;
   
    double C = 1 + A*deltaT ;
    double D = deltaT * (1 + 0.5 * A * deltaT);
    
    int i,j;
    real dr[NDIM+1], r, rr, ri, rrCut;
    double vv;

    double uVal, uSum, fcVal, f, AA, AASum;
    double TVal, T_Config;

    double deno, VVSum;	 
    double PE;
    PE=0.;
    deno = 0.;
    VVSum = 0.;
    keConfig = 0.;
    AASum = 0.;
   
  for(n=1;n<=nAtom; n++)
     TValSum[n] = 0.;

  rrCut = Sqr(rCut);

/*****Calculating Configarational temperature*****/
if(thermo == 'C'){
for(i = 1 ; i <= nAtom; i ++){
    for(j = i+1 ; j <= nAtom ; j ++){
      dr[1] = rx[i] - rx[j];
      if(fabs(dr[1]) > regionH[1])
	dr[1] -= SignR(region[1], dr[1]);

     dr[2] = ry[i] - ry[j];
      if(fabs(dr[2]) > regionH[2])
	dr[2] -= SignR(region[2], dr[2]);
  
     rr = Sqr(dr[1]) + Sqr(dr[2]);
     if(rr < rrCut ){
     r = sqrt(rr);
     ri = 1/r;
     uVal = ri*exp(-kappa*r);

     TVal = (1./rr + Sqr(kappa) + kappa/r)*uVal;
     TValSum[i] += TVal;
     TValSum[j] += TVal;
   } }
     AA = Sqr(ax[i]) + Sqr(ay[i]);
     AASum += AA;
     vv = Sqr(vx[i]) + Sqr(vy[i]);
     VVSum += vv; 
     deno += TValSum[i];    
} 
     PE = uSum/nAtom; 
     keConfig = AASum/deno;

     double gSumconfig, varSconfig, massSconfig;
     if(stepCount == 1) varSconfig = 0.;
     gSumconfig = 0.; massSconfig = 2.0; 
   
     gSumconfig = (AASum/temperature - deno)/massSconfig;
     varSconfig += deltaT*gSumconfig; 

      /*****Configarational Nose-Hoover thermostat*****/
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * (vx[n] + varSconfig * ax[n]);
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * (vy[n] + varSconfig * ay[n]);
    }
      /*****Kinetic Nose-Hoover thermostat*****/
  }else if(thermo == 'N'){  
    gSum = (0.5*S2 - (nAtom + 1)*temperature)/massS;
    varS += deltaT*gSum; 
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * (ax[n] - varS *vx[n]);
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * (ay[n] - varS *vy[n]);
      ry[n] += deltaT * vy[n];
   }
      /*****for Gaussian thermostat*****/
 }else if(thermo == 'G'){              
      for (n = 1; n <= nAtom; n++){
      vx[n] = C * vx[n] + D * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] = C * vy[n] + D * ay[n];
      ry[n] += deltaT * vy[n];
     }
   }else if (thermo == 'L'){
  double nu = 0.03066;
  double var = sqrt(2*nu/(GAMMA*deltaT));
  double scale = 1. + nu*deltaT/2.;
  double scale_v = 2./scale - 1.;
  double scale_f = deltaT/scale;
  int n, k;
  for(n = 1 ; n <= nAtom ; n ++){
      vx[n] = scale_v*vx[n] + scale_f*(ax[n] + var*gsl_ran_gaussian(rnd,1));
      rx[n] += deltaT * vx[n];
      vy[n] = scale_v*vy[n] + scale_f*(ay[n] + var*gsl_ran_gaussian(rnd,1));
      ry[n] += deltaT * vy[n];
    }
  }
 }else{
    int n;
    for(n = 1 ; n <= nAtom ; n ++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * vy[n];
    }
  }
}
  
