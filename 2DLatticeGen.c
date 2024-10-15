#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Sqr(x) ((x)*(x))

#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9
double RandR(int *seed);
double RandR(int *seed){
  *seed = (*seed * IMUL + IADD) & MASK;
  return (*seed * SCALE);
}

int main(int argc, char **argv){

  int NX;
  int NY;
  double density;
  double GAMMA;

  char dummy[25];
  FILE *fp;
  fp = fopen("STATE-PARAMS","r");
  fscanf(fp, "%s %d", dummy, &NX);
  fscanf(fp, "%s %d", dummy, &NY);
  fscanf(fp, "%s %lf", dummy, &density);
  fscanf(fp, "%s %lf", dummy, &GAMMA);
  fclose(fp);

  // SET THE COORDINATES ON A SQUARE LATTICE
  double *rx, *ry;
  rx = (double *)malloc((NX*NY+1)*sizeof(double));
  ry = (double *)malloc((NX*NY+1)*sizeof(double));
  double c[3], gap[3], region[3];
  int n,nx,ny;
  
  region[1] = NX/sqrt(density);
  region[2] = NY/sqrt(density);
  
  gap[1] = region[1]/NX;
  gap[2] = region[2]/NY;
  n = 0;
  for(ny = 1 ; ny <= NY ; ny++){
    c[2] = (ny - 0.5) * gap[2] - 0.5*region[2];
    for(nx = 1 ; nx <= NX ; nx++){
      c[1] = (nx - 0.5) * gap[1] - 0.5*region[1];
      n ++;
      rx[n] = c[1];
      ry[n] = c[2];
    }
  }

  int nAtom = n;
  // SET THE VELOCITIES FOR ALL PARTICLES
  double *vx, *vy;		
  vx = (double *)malloc((nAtom+1)*sizeof(double));
  vy = (double *)malloc((nAtom+1)*sizeof(double));
  double vSum[3], ang;
  int randSeed = 21;
  vSum[1] = vSum[2] = 0.0;
  double vMag = sqrt(2*(1.-1./nAtom)/GAMMA);
  for(n = 1 ; n <= nAtom ; n++){
    ang = 2 * M_PI * RandR(&randSeed);
    vx[n] = vMag * cos(ang);
    vy[n] = vMag * sin(ang);
    vSum[1] += vx[n];
    vSum[2] += vy[n];
  }
  vSum[1] /= nAtom;
  vSum[2] /= nAtom;
  for(n = 1 ; n <= nAtom ; n++){
    vx[n] -= vSum[1];
    vy[n] -= vSum[2];
  }

  // DUMP THE STATE ON THE 'STATE' FILE
  FILE *fpSTATE;
  fpSTATE = fopen("STATE","w");
  fprintf(fpSTATE,"timeNow 0\n");
  fprintf(fpSTATE,"nAtom %d\n", nAtom);
  fprintf(fpSTATE,"region[1] %E\n", region[1]);
  fprintf(fpSTATE,"region[2] %E\n", region[2]);

  for(n = 1 ; n <= nAtom ; n ++)
    fprintf(fpSTATE,"%d %E %E %E %E\n", n, rx[n], ry[n], vx[n], vy[n]);
  fclose(fpSTATE);

  free(rx);
  free(ry);
  free(vx);
  free(vy);
  return 0;
}
 
