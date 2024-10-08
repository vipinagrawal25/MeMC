include hosts/quest#
opt=-O3
# opt=-pg
ifeq ($(debug), y)
	opt = -g3  -Wall -Wpedantic
endif

link = $(opt) -lm -std=c++17 -lhdf5 -Isrc#
sources =  src/vector.cpp src/metropolis.cpp src/random_gen.cpp src/bending.cpp
sources += src/stretching.cpp src/hdf5_io.cpp src/misc.cpp
sources += src/electrostatics.cpp
# sources += src/misc.cpp
#
object = obj/vector.o obj/metropolis.o obj/random_gen.o obj/bending.o
object += obj/stretching.o obj/hdf5_io.o obj/misc.o
object += obj/multicomp.o
object += obj/electrostatics.o
#
# Find directories that match the pattern (00000, 00001, 00002, etc.)
DIRS := $(shell find . -maxdepth 1 -type d -name '[0-9][0-9][0-9][0-9][0-9]')

all : memc
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; fi
#
obj/readnml.o: src/read_namelist.f90
	gfortran -c -Jobj src/read_namelist.f90 -o obj/readnml.o
#
memc: $(object) obj/main.o obj/readnml.o
	$(CC) $(object) obj/main.o obj/readnml.o  $(link) -lgfortran -o exe_memc
#
obj/main.o: main.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)

object : $(object)
obj/%.o : src/%.cpp $(includes)
	@mkdir -p $(@D)
	$(info Compiling $<)
	$(CC) -Iobj -c $< -o $@ $(link)
#
clean:
	@rm -rf bin $(object) exe_*
	@echo "all obj bin cleared"

distclean:
	@rm -rf bin $(object) 
	@echo "all data cleared"
	@echo "Deleting .out files in directories named like 00000, 00001, etc..."
	@for dir in $(DIRS); do find "$$dir" -maxdepth 1 -name "$(OUT_PATTERN)" -exec rm -f {} \; -print; done
	@echo "Deleting files matching $(SNAP_PATTERN) in directories named like 00000, 00001, etc..."
	@for dir in $(DIRS); do find "$$dir" -maxdepth 1 -name "$(SNAP_PATTERN)" -exec rm -f {} \; -print; done
	@echo "Deleting mc_log and restartindex.txt in directories named like 00000, 00001, etc..."
	@for dir in $(DIRS); do find "$$dir" -maxdepth 1 -name "mc_log" -exec rm -f {} \; -print; done
	@for dir in $(DIRS); do find "$$dir" -maxdepth 1 -name "restartindex.txt" -exec rm -f {} \; -print; done

curv: $(object) obj/curv.o
	echo $(CC) $(object) $(link)
	@$(CC) $(object) obj/readnml.o obj/curv.o $(link) -lgfortran -o exe_curv
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe_curv $(bindir)

obj/curv.o: utils/getcurv.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)
##
gradient : $(object) obj/gradient.o
	$(CC) $(object) obj/gradient.o obj/readnml.o  $(link) -lgfortran -o exe_grad

obj/gradient.o: tests/gradient.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)

# Define the patterns for the files to delete
OUT_PATTERN = *.out
SNAP_PATTERN = snap_*.h5
OTHER_FILES = mc_log restartindex.txt