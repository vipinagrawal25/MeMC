HOST=su
include hosts/$(HOST)
#
#opt=-O3
opt=-pg
# ifeq ($(debug), y)
# 	opt = -g3  -Wall -pedantic
# endif

link = $(opt) -lm -std=c++17 -lhdf5 -Iincludes
# 
sources = src/forces_lj.cpp src/forces_surf.cpp src/Metropolis.cpp
sources += src/init.cpp src/hdf5_io.cpp
sources += src/cubic_solve.cpp
sources += src/misc.cpp
sources += src/visit_io.cpp
#
object =  obj/forces_lj.o obj/init.o obj/forces_surf.o obj/Metropolis.o
object += obj/hdf5_io.o 
object += obj/cubic_solve.o obj/misc.o
object += obj/visit_io.o
#
includes += includes/global.h includes/subroutine.h includes/Vector.h includes/misc.h
bindir = ./bin
#
#
all : start memc
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe* $(bindir)/

start: $(object) obj/start.o
	echo $(CC) $(object) $(link)
	@$(CC) $(object) obj/start.o $(link) -o exe_start

memc: $(object) obj/memc.o
	echo $(CC) $(object) $(link)
	@$(CC) $(object) obj/memc.o $(link) -o exe_memc

energy: $(object) obj/energy.o
	echo $(CC) $(object) $(link)
	@$(CC) $(object) obj/energy.o $(link) -o exe_energy
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe_energy $(bindir)

obj/memc.o: mains/memc.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)
#
obj/start.o: mains/start.cpp $(includes)
	$(CC) -Jobj -c $< -o $@ $(link)
#
obj/energy.o: utils/energy.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)

object : $(object)
obj/%.o : src/%.cpp $(includes)
	@mkdir -p $(@D)
	$(info Compiling $<)
	$(CC) -Iobj -c $< -o $@ $(link)
#
clean:
	@rm -rf bin $(object)
	@echo "all obj bin cleared"

distclean:
	@rm -rf bin $(object) 
	@echo "all data cleared"
	@rm -rf $(dir $(wildcard */mc_log))
	@echo "Deleted all the output data"