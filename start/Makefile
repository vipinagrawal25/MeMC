CC=/software/mpi/openmpi-4.1.4-gcc-11.2.0/bin/mpic++ -O3 -w  -L/software/hdf5/1.12.0/lib/ -I/software/hdf5/1.12.0/include/
opt=-O3
ifeq ($(debug), y)
	opt = -g3  -Wall -pedantic
endif

link = $(opt) -lm -std=c++17 -lhdf5 -Iincludes # 
sources = forces_lj.cpp Metropolis.cpp
sources += init.cpp  hdf5_io.cpp
#
object =  obj/forces_lj.o obj/init.o  obj/Metropolis.o
object += obj/hdf5_io.o
#
includes += global.h subroutine.h 
bindir = ./bin
#
#
all : objdir start 

objdir:
	@mkdir -p obj

start: $(object) obj/start.o 
	$(CC) $(object) obj/start.o  $(link)  -o exe_start


obj/start.o: start.cpp $(includes)
	$(CC) -Jobj -c $< -o $@ $(link)

object : $(object) 
obj/%.o : %.cpp $(includes) 
	$(info Compiling $<)
	$(CC) -Iobj -c $<  -o $@ $(link)
#
clean:
	@rm -rf bin $(object) exe_*
	@echo "all obj bin cleared"

distclean:
	@rm -rf bin $(object) 
	@echo "all data cleared"
	@rm -rf $(dir $(wildcard */mc_log))
	@echo "Deleted all the output data"
