#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <stdbool.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <iostream>

#include "metropolis.hpp"
#include "bending.hpp"
#include "stretching.hpp"
#include "random_gen.hpp"
#include "hdf5_io.hpp"
#include "misc.hpp"
#include "multicomp.hpp"
#include "electrostatics.hpp"
#include "selfavoidance.hpp"
#include "linetension.hpp"

template<typename T>
string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}
/*----------------------------------------------------------*/
void start_simulation(MESH_p &mesh, McP mcobj, STE &stretchobj, string outfolder, double radius, int &residx){

    double Pole_zcoord;
    double ave_bond_len;
    int fnumber, titer;
    string resfile;

    stretchobj.init_eval_lij_t0(mesh, mcobj.isfluid());
    if(!mcobj.isrestart()) { residx = 0; }
    else{
        fstream restartfile(outfolder+"/restartindex.txt", ios::in);
        if (restartfile.is_open()) {
            restartfile >> titer >> fnumber;
            restartfile.close();
        } else {
            cerr << "Error: restartindex.txt does not exist or could not be opened." 
            << "\nPlease change the para_file to start code from begining." << std::endl;
            exit(EXIT_FAILURE);
        }
        residx = titer;
        resfile=outfolder+"/snap_"+ZeroPadNumber(fnumber)+".h5";
        hdf5_io_read_double( (double *)mesh.pos,  resfile, "pos");
        hdf5_io_read_mesh((int *) mesh.numnbr, (int *) mesh.node_nbr_list, resfile);
    }
}
/*----------------------------------------------------------*/
int main(int argc, char *argv[]){
    int mpi_err, mpi_rank, residx, world_size;

    mpi_err = MPI_Init(0x0, 0x0);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    mpi_err =  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    pid_t pid = getpid();
    uint32_t seed_v;
    int iter, start, num_moves, num_bond_change, recaliter, num_exchange=0;
    int fnumber=0;
    double av_bond_len, Etot;
    clock_t timer;
    timer = clock(); // Initialize timer
    string outfolder, para_file, outfile, filename;

    start=0;
    outfolder = ZeroPadNumber(mpi_rank+start)+"/";
    fstream fileptr(outfolder+"/mc_log", ios::app);

    // Check if the file opened successfully
    seed_v = (uint32_t) (mpi_rank + time(0));
    RandomGenerator::init(seed_v);

    MESH_p mesh(outfolder);

    BE bendobj(mesh, outfolder);
    STE stretchobj(mesh, outfolder);
    MulCom lipidobj(mesh, outfolder);
    ESP chargeobj(mesh, outfolder);
    SelfAvoid repulsiveobj(mesh, outfolder);
    LTN lineobj(mesh, outfolder);

    McP mcobj(bendobj, stretchobj, lipidobj, chargeobj, repulsiveobj, lineobj);
    mcobj.initMC(mesh, outfolder);
    //
    if(!mcobj.exchange()) mesh.ncomp=1;
    if (mcobj.isfluid()) recaliter=mcobj.fluidizeevery();
    else recaliter=1;

    ostream* terminal;
    ofstream out_file;

    if (world_size == 1) terminal = &std::cout;
    else{
        out_file.open(outfolder+"/terminal.out", std::ios_base::app);
        terminal = &out_file;
    }
    
    (*terminal) << "# Simulation started on " << ctime(&timer) << endl;
    start_simulation(mesh, mcobj, stretchobj, outfolder, mesh.radius, residx);
    if(!mcobj.isrestart()) mcobj.wHeader(mesh,fileptr);
    if(mcobj.isrestart()) fileptr << "# Restart index " << residx << endl;

    (*terminal) << "# The seed value is " << seed_v << endl;
    Etot = mcobj.evalEnergy(mesh);
    mcobj.write_energy(fileptr, iter, mesh);

    //
    for(int cycle = 0; cycle<5; cycle++){
        mcobj.startcycle(cycle);
        for(int anneal=0; anneal < 6; anneal++){
            mcobj.updateparam(anneal, outfolder);
            for(iter=residx; iter < mcobj.totaliter(); iter++){
                if(iter%mcobj.dumpskip() == 0){
                    outfile=outfolder+"/snap_"+ZeroPadNumber(fnumber)+".h5";
                    hdf5_io_delete(outfile);
                    hdf5_io_write((double*) mesh.pos, 3*mesh.N, outfile, "pos");
                    hdf5_io_write_mesh(mesh.numnbr, mesh.node_nbr_list, mesh.N,
                                    mesh.nghst, outfile);
                    if (mesh.ncomp>1) hdf5_io_write(mesh.compA, mesh.N, outfile, "lip");
                    fstream restartfile(outfolder+"/restartindex.txt", ios::out);
                    restartfile << iter << " " << fnumber << endl;
                    restartfile.close();
                    fnumber++;
                    (*terminal) << "Snapshot written to " << outfile << endl;
                }
                if(repulsiveobj.isSelfRepulsive()) repulsiveobj.buildCellList(mesh);
                num_moves = mcobj.monte_carlo_3d(mesh.pos, mesh);
                if (mcobj.exchange()) num_exchange = mcobj.monte_carlo_lipid(mesh.pos, mesh);
                if (mcobj.isfluid() && !(iter % mcobj.fluidizeevery())){
                    num_bond_change = mcobj.monte_carlo_fluid(mesh.pos, mesh);
                    (*terminal) << "fluid stats " << num_bond_change
                                << " bonds flipped" << endl;
                }
                if(!(iter % recaliter)){
                    Etot = mcobj.evalEnergy(mesh);
                    (*terminal) << "iter = " << iter << 
                    "; Accepted Moves = " << (double)num_moves*100/mcobj.onemciter() << " %;"
                    "; Exchanged Moves = " << (double)num_exchange * 100 / mcobj.onemciter() << " %;"
                    << " totalener = " << Etot << "; volume = " << mcobj.getvolume() << endl;
                }
                mcobj.write_energy(fileptr, iter, mesh);
                residx=iter;
            }
        }
    }
    // Final snapshot
    (*terminal) << "Total time taken = " << (clock()-timer)/CLOCKS_PER_SEC << "s" << endl;
    if (world_size > 1) out_file.close();
    fileptr.close();
    mesh.free();
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_err = MPI_Finalize();
    return 0;
}