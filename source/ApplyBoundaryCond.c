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
#include<math.h>
#include<string.h>
#include"global.h"
void ApplyBoundaryCond(char BC){
  int n, vSign;
  for(n = 1 ; n <= nAtom ; n ++){
    // P.B.C along x and y axis
  if(BC == 'P'){
    rx[n] -= region[1]*rint(rx[n]/region[1]);
    ry[n] -= region[2]*rint(ry[n]/region[2]);
    }
   // R.B.C along y axis
 else if (BC == 'R'){
     rx[n] -= region[1]*rint(rx[n]/region[1]); 
     vSign = 0;
     if (ry[n] >= regionH[2]){
     ry[n] = regionH[2] * 0.999999;  vSign = -1;
     } else if (ry[n] <= -regionH[2]){
     ry[n] = -0.999999*regionH[2];       vSign = 1;  
     }if(vSign) {
     if (vy[n] * vSign < 0.)
     vy[n] = -vy[n];                                             
   }
  } 
 }
}
