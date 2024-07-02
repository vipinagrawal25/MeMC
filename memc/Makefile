# HOST=norlx65
# include hosts/$(HOST)
# CC = mpic++
CC=mpic++ -O3 -w  -L/home/vipin-agrawal/opt/hdf5/HDF_Group/HDF5/1.12.3/lib/ \
	-I/home/vipin-agrawal/opt/hdf5/HDF_Group/HDF5/1.12.3/include/
#
opt=-O3
# opt=-pg
ifeq ($(debug), y)
	opt = -g3  -Wall -Wpedantic
endif

link = $(opt) -lm -std=c++14 -lhdf5 -Isrc # 
sources =  src/vector.cpp src/metropolis.cpp src/random_gen.cpp src/bending.cpp
sources += src/stretching.cpp  src/hdf5_io.cpp  src/misc.cpp
# sources += src/misc.cpp
#
object = obj/vector.o obj/metropolis.o obj/random_gen.o obj/bending.o
object += obj/stretching.o  obj/hdf5_io.o obj/misc.o
bindir = ./bin
#
#
all : memc 
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe* $(bindir)/
#
obj/readnml.o: src/read_namelist.f90
	gfortran -c -Jobj src/read_namelist.f90 -o obj/readnml.o
#
memc: $(object) obj/main.o obj/readnml.o
	$(CC) $(object) obj/main.o obj/readnml.o  $(link) -lgfortran -o exe_memc
#
obj/main.o: main.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)
#
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
	@rm -rf $(dir $(wildcard */mc_log))
	@echo "Deleted all the output data"

curv: $(object) obj/curv.o
	echo $(CC) $(object) $(link)
	@$(CC) $(object) obj/readnml.o obj/curv.o $(link) -lgfortran -o exe_curv
	@if [ ! -d $(bindir) ] ; then echo "directory bin does not exist creating it" ; mkdir $(bindir) ; fi
	mv exe_curv $(bindir)

obj/curv.o: utils/getcurv.cpp $(includes)
	@$(CC) -Jobj -c $< -o $@ $(link)
