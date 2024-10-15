CC	= mpiicc -O3 -g -Wall
INCLUDE = -I/soft/psxe2015/impi/5.0.1.035/intel64/include -I/home/hcharan/Installed-Softwares/gsl-2.3/include/
LIBS 	= -L/home/hcharan/Installed-Softwares/gsl-2.3/lib
LFLAG	= -lm -lgsl -lgslcblas	

OBJECTS	= AccumProps.o AllocArrays.o ApplyBoundaryCond.o Close.o    \
	  ComputeForcesCells.o EvalProps.o Init.o LeapfrogStep.o    \
	  main.o PrintSummary.o SetupJob.o Trajectory.o 	    \
	  EvalVacf.o AccumVacf.o InitVacf.o     \
	  ZeroVacf.o PrintVacf.o Integrate.o DumpState.o EvalRdf.o \
	  EvalDiffuse.o AccumDiffuse.o InitDiffuse.o ZeroDiffuse.o  \
	  PrintDiffuse.o EvalVeloAcf.o AccumVeloAcf.o InitVeloAcf.o \
 	  ZeroVeloAcf.o PrintVeloAcf.o EvalSpacetimeFluct.o AccumSpacetimeFluct.o \
	  InitSpacetimeFluct.o PrintSpacetimeFluct.o ZeroSpacetimeFluct.o \
          EvalSpacetimeFluct2.o AccumSpacetimeFluct2.o \
	  InitSpacetimeFluct2.o PrintSpacetimeFluct2.o ZeroSpacetimeFluct2.o\
	  EvalVelDist.o PrintVelDist.o


main	: $(OBJECTS) 
	$(CC) -o main $(INCLUDE) $(OBJECTS) $(LIBS) $(LFLAG)

%.o: 	%.c
	$(CC) -O -c $< $(INCLUDE) $(LIBS) $(LFLAG)

$(OBJECTS) : global.h

.PHONY : clean

clean :
	rm -rf main $(OBJECTS) 

