CC = g++ 
#
opt = -O3
ifeq ($(debug), y)
	opt = -g3  -Wall -pedantic
endif

link = $(opt) -lm -std=c++17 -lhdf5 -Iincludes
# 
sources = src/forces_lj.cpp src/forces_surf.cpp src/Metropolis.cpp
sources += src/init.cpp  src/vtk_io.cpp src/hdf5_io.cpp
sources += src/cubic_solve.cpp
sources += src/misc.cpp
#
object =  obj/forces_lj.o obj/init.o obj/forces_surf.o obj/Metropolis.o
object += obj/vtk_io.o obj/hdf5_io.o 
object += obj/cubic_solve.o obj/misc.o
#
includes += includes/global.h includes/subroutine.h includes/Vector.h 
bindir = ./bin
#
#
all : exe_start 
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe* $(bindir)/

exe_start: $(object) obj/start.o
	echo $(CC) $(object) $(link)
	$(CC) $(object) obj/start.o $(link) -o exe_start


exe_memc: $(object) obj/memc.o
	echo $(CC) $(object) $(link)
	$(CC) $(object) $(link) -o run

obj/memc.o: mains/main.cpp $(includes)
	$(CC) -Jobj -c $< -o $@ $(link)
# 

obj/start.o: mains/start.cpp $(includes)
	$(CC) -Jobj -c $< -o $@ $(link)

object : $(object)
obj/%.o : src/%.cpp $(includes)
	@mkdir -p $(@D)
	$(info Compiling $<)
	$(CC) -Iobj -c $< -o $@ $(link)
#
clean:
	@rm $(object) run
	@echo "all object file and executables removed"
#
distclean:
	@rm -rf data_0* $(object) run
	@echo "all data cleared"
	@rm -rf $(dir $(wildcard */mc_log))
	@echo "Deleted all the output data"
