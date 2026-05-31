CC = mpic++
#
opt=-O3
# opt=-pg
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
