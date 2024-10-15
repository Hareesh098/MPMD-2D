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
void ZeroDiffusion() ;
void PrintDiffusion() ;
void AccumDiffusion(){
 int j, nb;								
 real xfac,yfac,fac;
 for(nb = 1; nb <= nBuffAcf; nb ++){
 if (indexDiffuse[nb] == nValAcf){
 for (j = 1; j <= nValAcf; j ++ ){
 xDiffuseAv[j] += xDiffuse[nb][j] ;
 yDiffuseAv[j] += yDiffuse[nb][j] ;
 rDiffuseAv[j] += rDiffuse[nb][j] ;
}
 indexDiffuse[nb] = 0 ;
 countDiffuseAv ++ ;
 if ( countDiffuseAv == limitAcfAv){
 xfac = 1./(1. * 2.* nAtom * stepAcf * deltaT * limitAcfAv) ;
 yfac = 1./(1. * 2.* nAtom * stepAcf * deltaT * limitAcfAv) ;
 fac = 1./(NDIM * 2.* nAtom * stepAcf * deltaT * limitAcfAv) ;

 for(j = 1; j <= nValAcf; j++ ) MSDAv[j] = (rDiffuseAv[j] /(NDIM * 2. * nAtom*limitAcfAv)); //Calculating mean-square-displacement(MSD)
 
 for(j = 2; j <= nValAcf; j++){
 rDiffuseAv[j] = rDiffuseAv[j] * fac / (j - 1) ;
 xDiffuseAv[j] = xDiffuseAv[j] * xfac / (j - 1) ;
 yDiffuseAv[j] = yDiffuseAv[j] * yfac / (j - 1) ;
 }
 PrintDiffusion() ;
 ZeroDiffusion() ;
} } }
}
