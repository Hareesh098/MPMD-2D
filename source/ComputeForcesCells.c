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
void ComputeForcesCells(){
  double dr[NDIM+1], invWid[NDIM+1], shift[NDIM+1], f, fcVal, rr, rrCut, ri, r, uVal;
  int c, I, J, m1, m1X, m1Y, m2, m2X, m2Y, n, offset;
  int iofX[] = {0, 0, 1, 1, 0, -1, -1, -1, 0, 1},
      iofY[] = {0, 0, 0, 1 ,1, 1, 0, -1, -1, -1};
  
  rrCut = Sqr(rCut);
  invWid[1] = cells[1]/region[1];
  invWid[2] = cells[2]/region[2];
 
  for(n = nAtom+1; n <= nAtom+cells[1]*cells[2] ; n++)
    cellList[n] = 0;
  
  for(n = 1 ; n <= nAtom ; n ++){
    c = ((int)((ry[n] + regionH[2])*invWid[2]))*cells[1] + (int)((rx[n]+regionH[1])*invWid[1]) + nAtom+ 1;
    cellList[n] = cellList[c];
    cellList[c] = n;
  }
  
  for(n = 1 ; n <= nAtom ; n ++){
    ax[n] = 0.;
    ay[n] = 0.;
  }
  
  uSum = 0.0 ;
  virSum = 0.0;
  rfAtom = 0.0;

  int start = 1 + rank*(cells[2]/size);
  int end = (rank+1)*(cells[2]/size);

  for(m1Y = start ; m1Y <= end ; m1Y ++){
    for(m1X = 1 ; m1X <= cells[1] ; m1X ++){
      m1 = (m1Y-1) * cells[1] + m1X + nAtom;
      for(offset = 1 ; offset <= 9 ; offset ++){
	m2X = m1X + iofX[offset]; shift[1] = 0.;
	if(m2X > cells[1]){
	  m2X = 1; shift[1] = region[1];
	}else if(m2X == 0){
	  m2X = cells[1]; shift[1] = -region[1];
	}
	m2Y = m1Y + iofY[offset]; shift[2] = 0.;
	if(m2Y > cells[2]){
	  m2Y = 1; shift[2] = region[2];
	}else if(m2Y == 0){
	  m2Y = cells[2]; shift[2] = -region[2];
	}
	m2 = (m2Y-1)*cells[1] + m2X + nAtom;
	I = cellList[m1];
	while(I > 0){
	  J = cellList[m2];
	  while(J > 0){
	    if(m1 == m2 && J != I){
	      dr[1] = rx[I] - rx[J] - shift[1];
	      dr[2] = ry[I] - ry[J] - shift[2];
	      rr = Sqr(dr[1]) + Sqr(dr[2]);
	      if(rr < rrCut){
		r = sqrt(rr);
		ri = 1/r;
		uVal = ri*exp(-kappa*r);
		fcVal = ri*uVal*(ri+kappa);
		f = fcVal * dr[1];
		ax[I] += f;
		f = fcVal * dr[2];
		ay[I] += f;
		uSum +=  0.5 * uVal;
		virSum += 0.5 * fcVal * rr;
		rfAtom += 0.5 * dr[1] * fcVal * dr[2];
	      }
	    }else if(m1 != m2){
	      dr[1] = rx[I] - rx[J] - shift[1];
	      dr[2] = ry[I] - ry[J] - shift[2];
	      rr = Sqr(dr[1]) + Sqr(dr[2]);
	      if(rr < rrCut){
		r = sqrt(rr);
		ri = 1/r;
		uVal = ri*exp(-kappa*r);
		fcVal = ri*uVal*(ri+kappa);
		f = fcVal * dr[1];
		ax[I] += f;
		f = fcVal * dr[2];
		ay[I] += f;
		uSum +=  0.5 * uVal;
		virSum += 0.5 * fcVal * rr;
		rfAtom += 0.5 * dr[1] * fcVal * dr[2];
	      }
	    }
       	    J = cellList[J];
	  }
	  I = cellList[I];
	}
      }
    }
  }
}
	
