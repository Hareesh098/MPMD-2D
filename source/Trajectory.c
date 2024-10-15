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
void Trajectory(){
  int n;
  fprintf(fpxyz, "%d\n", nAtom);
  fprintf(fpxyz, "timeNow %lf region[1] %lf region[2] %lf\n", timeNow, region[1], region[2]);
  for(n = 1 ; n <= nAtom ; n ++)
    fprintf(fpxyz, "%d\t %lf\t %lf\t %lf\t %lf\n", atomType[n], rx[n], ry[n], vx[n], vy[n]);
  fprintf(fpxyz, "\n");
}
