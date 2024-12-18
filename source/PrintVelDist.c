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
void PrintVelDist(){
  real vBin;
  int n;
  fprintf(fpveldist,"velocity distribution at timenow = %lf\n", timeNow);
  for(n = 1; n <= sizeHistVel; n++){
   vBin = (n - 0.5) * rangeVel / sizeHistVel;
   fprintf(fpveldist,"%lf\t%lf\n" ,vBin, histVel[n]);
  }
  fprintf(fpveldist,"\n");
}
