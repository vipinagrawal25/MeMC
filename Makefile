HOST=norlx65
include hosts/$(HOST)
# CC = mpic++
#
opt=-O3
# opt=-pg
# ifeq ($(debug), y)
# 	opt = -g3  -Wall -pedantic
# endif

link = $(opt) -lm -std=c++17 -lhdf5 -Iincludes # 
sources = src/forces_lj.cpp src/forces_surf.cpp src/Metropolis.cpp
sources += src/init.cpp  src/hdf5_io.cpp
sources += src/cubic_solve.cpp
sources += src/misc.cpp src/random_gen.cpp
#
object =  obj/forces_lj.o obj/init.o obj/forces_surf.o obj/Metropolis.o
object += obj/hdf5_io.o 
object += obj/cubic_solve.o obj/misc.o obj/vector.o obj/random_gen.o
#
includes += includes/global.h includes/subroutine.h includes/Vector.h 
bindir = ./bin
#
#
all : start memc 
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe* $(bindir)/

obj/readnml.o: src/read_namelist.f90
	gfortran -c src/read_namelist.f90 -o obj/readnml.o

start: $(object) obj/start.o obj/readnml.o
	$(CC) $(object) obj/start.o obj/readnml.o $(link) -lgfortran -o exe_start

memc: $(object) obj/memc.o obj/readnml.o
	$(CC) $(object) obj/readnml.o obj/memc.o  $(link) -lgfortran -o exe_memc
#

obj/memc.o: mains/memc.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)
# 

obj/start.o: mains/start.cpp $(includes)
	$(CC) -Jobj -c $< -o $@ $(link)

object : $(object) obj/readnml.o
obj/%.o : src/%.cpp $(includes) obj/readnml.o
	@mkdir -p $(@D)
	$(info Compiling $<)
	$(CC) -Iobj -c $< obj/readnml.o -o $@ $(link)
#
clean:
	@rm -rf bin $(object) exe_*
	@echo "all obj bin cleared"

distclean:
	@rm -rf bin $(object) 
	@echo "all data cleared"
	@rm -rf $(dir $(wildcard */mc_log))
	@echo "Deleted all the output data"
