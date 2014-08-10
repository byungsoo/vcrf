########################
# Mex compile
########################
#CPP = mex
#
#CFLAGS =   $(SFLAGS) -g -v -DA64BITS -DMAT73 -DMEX_COMPILE -largeArrayDims #min added -DA64BITS
#LDFLAGS =  $(SFLAGS) -g -v -lm -Wall -fPIC -DA64BITS -DMAT73 -DMEX_COMPILE -largeArrayDims


########################
# g++ compile
########################
CPP = g++

CFLAGS =   $(SFLAGS) -O3 -fomit-frame-pointer -ffast-math -Wall -fPIC
LDFLAGS =  $(SFLAGS) -O3 -lm -Wall -fPIC


all: callable.o

clean:
	rm -f *.o *.tcov *.d core gmon.out *.stackdump

#example.o: GCoptimization.o example.cpp
#	$(CPP) -c $(CFLAGS) example.cpp
#	$(CPP) -c $(CFLAGS) example.cpp GCoptimization.o

GCoptimization.o: GCoptimization.cpp GCoptimization.h LinkedBlockList.o LinkedBlockList.h maxflow.o graph.o graph.h energy.h QPBO.o
	$(CPP) -c $(CFLAGS) GCoptimization.cpp -o GCoptimization.o
#	$(CPP) -c -DNO_UNARY_DET_COMPILE $(CFLAGS) GCoptimization.cpp -o GCoptimization.o
#	$(CPP) -c $(CFLAGS) GCoptimization.cpp LinkedBlockList.o maxflow.o graph.o -o GCoptimization.o

QPBO.o: QPBO.cpp QPBO.h QPBO_extra.o QPBO_maxflow.o QPBO_postprocessing.o
	$(CPP) -c $(CFLAGS) QPBO.cpp -o QPBO.o

QPBO_extra.o: QPBO_extra.cpp
	$(CPP) -c $(CFLAGS) QPBO_extra.cpp -o QPBO_extra.o

QPBO_maxflow.o: QPBO_maxflow.cpp
	$(CPP) -c $(CFLAGS) QPBO_maxflow.cpp -o QPBO_maxflow.o

QPBO_postprocessing.o: QPBO_postprocessing.cpp
	$(CPP) -c $(CFLAGS) QPBO_postprocessing.cpp -o QPBO_postprocessing.o

LinkedBlockList.o: LinkedBlockList.cpp LinkedBlockList.h QPBO.o
	$(CPP) -c $(CFLAGS) LinkedBlockList.cpp -o LinkedBlockList.o

maxflow.o: maxflow.cpp graph.h
	$(CPP) -c $(CFLAGS) maxflow.cpp -o maxflow.o


graph.o: graph.cpp graph.h block.h
	$(CPP) -c $(CFLAGS) graph.cpp -o graph.o



callable.o: callable.cpp GCoptimization.o QPBO.o
	$(CPP) -c $(CFLAGS)  callable.cpp -o callable.o
