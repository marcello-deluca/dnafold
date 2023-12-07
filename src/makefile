CPPC = g++ #clang++
CPPFLAGS = -O3 -freciprocal-math -funsafe-math-optimizations -fno-math-errno  -fassociative-math -fno-signed-zeros -fno-trapping-math -Wall -pedantic -std=c++11 -lpthread #-Werror

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SRCS))
#DBGOBJS = $(patsubst %.cpp, %.dbg.o, $(SRCS))

.PHONY: clean depend all

PROGRNAME = ../dnafold

all: $(PROGRNAME) #$(PROGRNAME)_dbg

#RULE:
$(PROGRNAME): $(OBJS)
	$(CPPC) -o $@ $(CPPFLAGS) $(OBJS)

#$(PROGRNAME)_dbg: $(DBGOBJS)
#	g++ -o $@ -ggdb3 $(DBGOBJS)

clean:
	rm -f $(PROGRNAME) $(PROGRNAME)_dbg *.o *.cpp~ *.hpp~

depend:
	makedepend $(SRCS)

dna_bd.o: importer_function.hpp recordFrame.hpp noOverlapsAtInit.hpp makeConnectivityMatrix.hpp parameters.hpp headers.hpp simInit.hpp resetForces.hpp resetStochasticForces.hpp calculateStochasticForces.hpp calculatePairForcesParallel.hpp calculateBendingForcesParallel.hpp integrateMotion.hpp calculateICSTorsionAndBending.hpp setParameters.hpp readInputFile.hpp parse_json.hpp
importer_function.o: importer_function.hpp

