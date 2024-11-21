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

void AccumDiffusion();
void EvalDiffusion(){
 int k, n, nb, ni;
 for(nb = 1; nb <= nBuffAcf; nb ++){
 indexDiffuse[nb] ++ ;
 if (indexDiffuse[nb] <= 0) continue ;
 if (indexDiffuse[nb] == 1){
 k = 1;
 for(n = 1; n <= nAtom; n++){
 rDiffuseTrue[nb][k] = rx[n] ;
 rDiffuseTrue[nb][k+1] = ry[n] ;
 rDiffuseOrg[nb][k] = rx[n] ;
 rDiffuseOrg[nb][k+1] = ry[n] ;
 k += 2 ;
 }}
 ni = indexDiffuse[nb] ;
 xDiffuse[nb][ni] = 0. ; yDiffuse[nb][ni] = 0. ; rDiffuse[nb][ni] = 0. ;
 k = 1 ;
 for(n = 1; n <= nAtom; n ++){
 rDiffuseTrue[nb][k]   = rx[n] +   rint((rDiffuseTrue[nb][k] - rx[n]) / region[1]) *region[1] ;
 rDiffuseTrue[nb][k+1] = ry[n] + rint((rDiffuseTrue[nb][k+1] - ry[n]) / region[2]) *region[2] ;
 rDiffuse[nb][ni] +=  Sqr(rDiffuseTrue[nb][k]  - rDiffuseOrg[nb][k] ) + Sqr(rDiffuseTrue[nb][k+1]  - rDiffuseOrg[nb][k+1] );
 xDiffuse[nb][ni] +=  Sqr(rDiffuseTrue[nb][k]  - rDiffuseOrg[nb][k] );
 yDiffuse[nb][ni] +=  Sqr(rDiffuseTrue[nb][k+1]  - rDiffuseOrg[nb][k+1] ) ;     
 k += 2 ;
 }}
 AccumDiffusion() ;
}
