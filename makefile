CC=g++
OPTIM=-O4
UTIL=util
CFLAGS=$(OPTIM) -c -I. -I$(UTIL)/include $(OPTIONS)
CSPICE_LIBS=$(UTIL)/lib/cspice.a $(UTIL)/lib/csupport.a
GSL_LIB=$(UTIL)/lib/libgsl.a $(UTIL)/lib/libgslcblas.a 
LFLAGS=-lm $(CSPICE_LIBS) $(GSL_LIB)
PROJECT_FILES=\
makefile README.md TODO \
cometsyn.cpp comet-simulation.cpp \
comet-analysis.py analysis.cfg comet-animation.sh \
plot-orbit.gpl plot-orbit-fragments.gpl plot-trajectories.gpl plot-fragments.gpl fragments.gph

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@

%.o:%.cpp cometsyn.cpp
	$(CC) $(CFLAGS) $< -o $@

cleanout:
	rm -rf *.out *.o *.exe

clean:cleanout
	rm -rf *.log *~ *.png *.dat \#*\#
	find . -name *~ -exec rm -rf {} \;
	rm -rf animation/*.png

edit:
	emacs -nw $(PROJECT_FILES)

compile:
	make comet-simulation.out

run:
	./comet-simulation.out

force:
	./comet-simulation.out yes
