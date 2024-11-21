/*-----------------------------------------------------------------------------

    MPMD-v2.0 : MULTI POTENTIAL MOLECULAR DYNAMICS 
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
    
    See the README file in the top-level MPMD directory.

-----------------------------------------------------------------------------*/



#include<stdio.h>
#include<math.h>
#include"global.h"
void AccumProps(int icode){
  if(icode == 0){
    sTotEnergy = ssTotEnergy = 0.;
    sKinEnergy = ssKinEnergy = 0.;
    sPressure = ssPressure = 0.;
    sPotEnergy = ssPotEnergy = 0.;
    sKEConfig = ssKEConfig = 0.;
    svirSum = 0.;
  }else if(icode == 1){
    sTotEnergy += totEnergy;
    ssTotEnergy += Sqr(totEnergy);
    sKinEnergy += kinEnergy;
    ssKinEnergy += Sqr(kinEnergy);
    sPressure += pressure;
    ssPressure += Sqr(pressure);
    sPotEnergy += potEnergy;
    ssPotEnergy += Sqr(potEnergy);
    sKEConfig += keConfig;
    ssKEConfig += Sqr(keConfig);
    svirSum += virSum;
  }else if(icode == 2){
    sTotEnergy /= stepAvg;
    ssTotEnergy = sqrt(ssTotEnergy/stepAvg - Sqr(sTotEnergy));
    sKinEnergy /= stepAvg;
    ssKinEnergy = sqrt(ssKinEnergy/stepAvg - Sqr(sKinEnergy));
    sPressure /= stepAvg;
    ssPressure = sqrt(ssPressure/stepAvg - Sqr(sPressure));
    sPotEnergy /= stepAvg;
    ssPotEnergy = sqrt(ssPotEnergy/stepAvg - Sqr(ssPotEnergy));
    sKEConfig /= stepAvg;
    ssKEConfig = sqrt(ssKEConfig/stepAvg - Sqr(sKEConfig));
    svirSum /= stepAvg;
  }
}
