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



typedef double real;

#define	NDIM 2
#define Sqr(x) ((x) * (x))
#define SignR(x, y) (((y) >= 0) ? (x) : (- (x)))

real	*rx, *ry, *vx, *vy, *ax, *ay;
real    *A, *Ax, *Ay;
real    *TValSum;

real	region[NDIM+1], regionH[NDIM+1], GAMMA, temperature, 
	deltaT, density, kinEnergy, potEnergy, pressure, rCut, kappa,
	sKinEnergy, sPressure, sTotEnergy, ssKinEnergy, ssPressure, 
	ssTotEnergy, timeNow, totEnergy, uSum, virSum, svirSum,
	sPotEnergy, ssPotEnergy, potEnergy, 
	vSum, vvSum, keConfig, sKEConfig, ssKEConfig;

char    BC, thermo;

int	*atomType;


int	initUcell[NDIM+1], moreCycles, nAtom, stepAvg, 
	stepCount, stepEquil, stepLimit, stepTrajectory, stepDump;

int	*cellList, cells[NDIM+1];

int 	rank, size, master;
real	*fax, *fay, fuSum, fvirSum;

// Spacetime correlations for k = k_x
real	**spacetimeFluct,**complxFuncOrg, *complxFunc,*spacetimeFluctAv;
int	*indexFluct, countFluctAv;
int	nFunCorr;

// Spacetime Correlations for k = k_y
real	**spacetimeFluct2,**complxFuncOrg2, *complxFunc2,*spacetimeFluctAv2;
int	*indexFluct2, countFluctAv2;
int      limitCorrAv, nBuffCorr, nFunCorr2, nValCorr, stepCorr;

// Viscosity
real  	rfAtom, frfAtom;
real  	 **viscAcf, *viscAcfOrg, *viscAcfAv, viscAcfInt;
int     nValAcf, nBuffAcf, stepAcf, countAcfAv, limitAcfAv, *indexAcf;

// Radial distribution function
real 	**histRdf,*HistRdf, rangeRdf, rBin;
int	countRdf, limitRdf, sizeHistRdf, stepRdf;

// Diffusion by Einstein's defination 
real	**rDiffuseTrue, **rDiffuseOrg, **xDiffuse, **yDiffuse, **rDiffuse, *xDiffuseAv,*yDiffuseAv,*rDiffuseAv,*MSDAv;
int	*indexDiffuse, countDiffuseAv ;

// Velocity auto correlation function
real	**VeloAcfOrg, **VeloAcf, *VeloAcfAv, diffuseVeloAcfInt;
int	*indexVeloAcf, countVeloAcfAv ; 

// Velocity distribution
real *histVel, rangeVel;
int countVel, limitVel, sizeHistVel, stepVel;


// Output files prefixes
char    *prefix;
char    dirprefix[250];

char    result[250];
FILE    *fpresult;

char 	xyz[256];
FILE 	*fpxyz;

char 	dnsty[256];
FILE 	*fpdnsty;

char    visc[256];
FILE 	*fpvisc;

char    rdf[256];
FILE    *fprdf;

char    rdfAvg[256];
FILE    *fprdfAvg;

char    diffuse[256];
FILE    *fpdiffuse;

char    veloAcf[256];
FILE    *fpveloAcf;

char    dnstyfluctKx[256];
FILE    *fpdnstyfluctKx;

char    longfluctKx[256];
FILE    *fplongfluctKx;

char    transfluctKx[256];
FILE    *fptransfluctKx;

char    dnstyfluctKy[256];
FILE    *fpdnstyfluctKy;

char    longfluctKy[256];
FILE    *fplongfluctKy;

char    transfluctKy[256];
FILE    *fptransfluctKy;

char    veldist[256];
FILE    *fpveldist;
