CC=g++
OPTIM=-O4
UTIL=util
CFLAGS=$(OPTIM) -c -I. -I$(UTIL)/include $(OPTIONS)
CSPICE_LIBS=$(UTIL)/lib/cspice.a $(UTIL)/lib/csupport.a
GSL_LIB=$(UTIL)/lib/libgsl.a $(UTIL)/lib/libgslcblas.a 
LFLAGS=-lm $(CSPICE_LIBS) $(GSL_LIB)

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@

%.o:%.cpp util.cpp
	$(CC) $(CFLAGS) $< -o $@

cleanout:
	rm -rf *.out *.o *.exe

clean:cleanout
	rm -rf *.log *~ *.png *.dat \#*\#
	find . -name *~ -exec rm -rf {} \;
	rm -rf animation/*

compile:
	make comet-cloud.out

run:
	./comet-cloud.out
