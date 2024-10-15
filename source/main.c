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
#include<mpi.h>
#include<time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include"global.h"

const gsl_rng_type * T;
gsl_rng * rnd;

void Init();
void SetupJob();
void EvalSpacetimeCorr();
void Trajectory();
void DumpState();
void ComputeForcesCells();
void LeapfrogStep(char thermo, gsl_rng *);
void ApplyBoundaryCond(char BC);
void EvalProps();
void EvalVacf();
void EvalRdf();
void AccumProps(int icode);
void PrintSummary();
void Close();
void EvalDiffusion();
void EvalVeloAcf();
void EvalSpacetimeFluct ();
void EvalSpacetimeFluct2 ();
void EvalVelDist();

int main(int argc, char **argv){
  time_t t1, t2;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  master = 0;
  
  if(rank == master){
     sprintf(dirprefix,"../output/");
    prefix = strcat(dirprefix, argv[1]);
    sprintf(result, "%s.result", prefix);
    fpresult = fopen (result, "w");
    sprintf (xyz, "%s.xyz", prefix);
    fpxyz = fopen (xyz, "w");
    sprintf (dnsty, "%s.curr-dnsty", prefix);
    fpdnsty = fopen (dnsty, "w");
    sprintf (visc, "%s.viscosity", prefix);
    fpvisc = fopen (visc, "w");
    sprintf (rdf, "%s.rdf", prefix);
    fprdf = fopen (rdf, "w");
    sprintf (diffuse, "%s.diffuse", prefix);
    fpdiffuse = fopen (diffuse, "w");
    sprintf (veloAcf, "%s.veloAcf", prefix);
    fpveloAcf = fopen (veloAcf, "w");
    sprintf (dnstyfluctKx, "%s.dnstyfluct-kx", prefix);
    fpdnstyfluctKx = fopen (dnstyfluctKx, "w");	
    sprintf (longfluctKx, "%s.longfluct-kx", prefix);
    fplongfluctKx = fopen (longfluctKx, "w");
    sprintf (transfluctKx, "%s.transfluct-kx", prefix);
    fptransfluctKx = fopen (transfluctKx, "w");	
    sprintf (dnstyfluctKy, "%s.dnstyfluct-ky", prefix);
    fpdnstyfluctKy = fopen (dnstyfluctKy, "w");	
    sprintf (longfluctKy, "%s.longfluct-ky", prefix);
    fplongfluctKy = fopen (longfluctKy, "w");
    sprintf (transfluctKy, "%s.transfluct-ky", prefix);
    fptransfluctKy = fopen (transfluctKy, "w");		
    sprintf (veldist, "%s.veldist", prefix);
    fpveldist = fopen (veldist, "w");
    sprintf (rdfAvg, "%s.rdfAvg", prefix);
    fprdfAvg = fopen (rdfAvg, "w");

  }

  Init();

  if( (cells[2] % size) != 0 ){
    if(rank == master){
      fprintf(fpresult, "\n size = %d\n", size);
      fprintf(fpresult, "cells[2] = %d\n", cells[2]);
      fprintf(fpresult, "\ncells[2] not exactly divisible by size\n");
      fprintf(fpresult, "Exiting .. \n\n");
      MPI_Finalize();
      exit(1);
    }else{
      MPI_Finalize();
      exit(1);
    }
  }

  SetupJob();

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rnd = gsl_rng_alloc (T);
  gsl_rng_set (rnd, time (0));

  if(rank == master)
    t1 = time(NULL);
  moreCycles = 1;
  while(moreCycles){
    stepCount ++;
    timeNow = stepCount*deltaT;
    
    ComputeForcesCells();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(ax, fax, (nAtom+1), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    MPI_Reduce(ay, fay, (nAtom+1), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    MPI_Reduce(&uSum, &fuSum, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    MPI_Reduce(&virSum, &fvirSum, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    MPI_Reduce(&rfAtom, &frfAtom, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);

    if(rank == master){
      memcpy(ax, fax, (nAtom+1)*sizeof(real));
      memcpy(ay, fay, (nAtom+1)*sizeof(real));
      memcpy(&uSum, &fuSum, sizeof(real));
      memcpy(&virSum, &fvirSum, sizeof(real));
      memcpy(&rfAtom, &frfAtom, sizeof(real));
      LeapfrogStep(thermo,rnd);
      ApplyBoundaryCond(BC);
      EvalProps();
      AccumProps(1);
      if(stepCount % stepAvg == 0){
	AccumProps(2);
	PrintSummary();
  	AccumProps(0);
      }
      if (stepCount >= stepEquil && (stepCount-stepEquil)%stepCorr == 0)
	EvalSpacetimeFluct2 (); 
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepAcf == 0)
	EvalVacf();
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepRdf == 0)
	EvalRdf();
      if(stepCount % stepTrajectory == 0)
	Trajectory();
      if(stepCount % stepDump == 0)
	DumpState();
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepAcf == 0)
	EvalDiffusion();   
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepAcf == 0)
	EvalVeloAcf(); 
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepAcf == 0)
	EvalSpacetimeFluct ();	
      if(stepCount >= stepEquil && (stepCount - stepEquil)%stepVel == 0)
	EvalVelDist(); 	
    }

    MPI_Bcast(rx, (nAtom+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    MPI_Bcast(ry, (nAtom+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

    if(stepCount >= stepLimit)
      moreCycles = 0;
  }

  if(rank == master){
    t2 = time(NULL);
    fprintf(fpresult, "Execution time %lf secs\n", difftime(t2,t1));
  }

  if(rank == master){
    fclose(fpresult);
    fclose(fpxyz);
    fclose(fpdnsty);
    fclose(fpvisc);
    fclose(fprdf);
    fclose(fpdiffuse);	
    fclose(fpveloAcf);
    fclose(fpdnstyfluctKx);
    fclose(fplongfluctKx);
    fclose(fptransfluctKx);
    fclose(fpdnstyfluctKy);
    fclose(fplongfluctKy);
    fclose(fptransfluctKy);
    fclose(fpveldist);
    fclose(fprdfAvg);	 
  }

  Close();
  MPI_Finalize();
  return 0;
}
