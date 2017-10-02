include compiler.in

BRANCH=$(shell bash .getbranch)
UTIL=util
CFLAGS=$(OPTIM) -c -I. -w -I$(UTIL)/include $(OPTIONS)
LFLAGS=-lm $(UTIL)/lib/$(ARCH)/cspice.a $(UTIL)/lib/$(ARCH)/csupport.a $(UTIL)/lib/$(ARCH)/libgsl.a $(UTIL)/lib/$(ARCH)/libgslcblas.a

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

commit:
	@echo "Commiting changes..."
	@touch .htaccess
	@-git commit -am "Commit"
	@git push origin $(BRANCH)

pull:
	@echo "Pulling from repository..."
	@git reset --hard HEAD	
	@git pull origin $(BRANCH)

show:
	@echo "Branch: $(BRANCH)"

