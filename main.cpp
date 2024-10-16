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

template<typename T>
string ZeroPadNumber(T num){
    ostringstream ss;
    ss << setw( 5 ) << setfill( '0' ) << (int)num;
    return ss.str();
}

void scale_pos(Vec3d *pos, double R, int N){
  for(int i = 0; i<N; i++) pos[i] = pos[i]*R;
}

bool isPlaner(Vec3d *pos, int Np){
    for(int i=0;i<Np;i++){if(pos[i].z!=0){return false;}}
    return true;
}

pair<double, double> get_box_dim(MESH_p mesh){
    Vec3d *Pos = mesh.pos;
    double xmin=1e7,xmax=-1e7,ymin=1e7,ymax=-1e7;
    for (int i = 0; i < mesh.N; ++i){
        if (Pos[i].x<xmin)  xmin=Pos[i].x;
        if (Pos[i].x>xmax)  xmax=Pos[i].x;
        if (Pos[i].y<ymin)  ymin=Pos[i].y;
        if (Pos[i].y>ymax)  ymax=Pos[i].y;
    }
    double xlen = xmax-xmin;
    double ylen = ymax-ymin;
    return {xlen, ylen};
}

double start_simulation(MESH_p &mesh, McP mcobj, STE &stretchobj, 
    string outfolder, double radius, int &residx){

    double Pole_zcoord;
    double ave_bond_len;
    int tdumpskip, titer;
    string resfile;

    hdf5_io_read_double( (double *)mesh.pos,  outfolder+"/input.h5", "pos" );
    scale_pos(mesh.pos, radius, mesh.N);
    hdf5_io_read_mesh((int *) mesh.numnbr, (int *) mesh.node_nbr_list, 
        outfolder+"/input.h5");

    if (isPlaner(mesh.pos,mesh.N)){
        mesh.sphere=false;
        mesh.edge = get_nstart(mesh.N, 1);
        // To make sure that the edges do not coincide after pbc.
        mesh.boxlen=get_box_dim(mesh).first*(1+1/sqrt(mesh.N));
    }
    else{
        mesh.sphere=true;
        mesh.edge = -1;
        mesh.boxlen=0;
        mesh.bdry_type=2;   // Sphere always has a pbc.
    }
    ave_bond_len = stretchobj.init_eval_lij_t0(mesh, mcobj.isfluid());
    // cout << ave_bond_len << endl;
    if(!mcobj.isrestart()) { residx = 0; }
    else{
        fstream restartfile(outfolder+"/restartindex.txt", ios::in);
        restartfile >> titer >>  tdumpskip;
        restartfile.close();
        residx = titer;
        // if(area_para.do_area)init_area_t0(Pos,mesh,mbrane_para,area_para);
        // init_spcurv(spcurv_para, Pos, mbrane_para.N);
        // if(stick_para.do_stick)
        // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, 
        // mbrane_para.N);
        // max(&mesh.nPole,&Pole_zcoord,Pos,mbrane_para.N);
        // min(&mesh.sPole,&Pole_zcoord,Pos,mbrane_para.N);
        // identify_attractive_part(Pos, stick_para.is_attractive, stick_para.theta, 
        //  mbrane_para.N);
        resfile=outfolder+"/snap_"+ZeroPadNumber(titer/tdumpskip)+".h5";
        hdf5_io_read_double( (double *)mesh.pos,  resfile, "pos");
        hdf5_io_read_mesh((int *) mesh.numnbr, (int *) mesh.node_nbr_list, resfile);
    }

    // init_activity(act_para, mbrane_para.N);
    return ave_bond_len;
}
/*----------------------------------------------------------*/
void diag_wHeader(BE bendobj, STE steobj, ESP chargeobj, MESH_p mesh,
                std::fstream &fid ){
    std::string log_headers = "#iter acceptedmoves bend_e stretch_e ";
    if (chargeobj.calculate()) log_headers+="electroe ";
    if(steobj.dopressure()) {log_headers+=" Pressure_e ";}
    if(steobj.dovol()) {log_headers+=" Volume_e ";}
    log_headers+="total_e ";
    if (mesh.sphere){log_headers+="volume";}
    fid << log_headers << endl;
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
    double av_bond_len, Etot;

    // Vec3d *Pos;
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

    // makemcobj();
    McP mcobj(bendobj, stretchobj, lipidobj, chargeobj, repulsiveobj);
    mcobj.initMC(mesh, outfolder);

    mesh.av_bond_len = start_simulation(mesh, mcobj, stretchobj, outfolder,
                        mesh.radius, residx);

    // How often do you want to compute the total energy?
    if (mcobj.isfluid()) recaliter=mcobj.fluidizeevery();
    else recaliter=10;

    ostream* terminal;
    ofstream out_file;

    if (world_size == 1) terminal = &std::cout;
    else{
        out_file.open(outfolder+"/terminal.out", std::ios_base::app);
        terminal = &out_file;
    }

    // ofstream terminal(outfolder+"/terminal.out", ios::app);
    (*terminal) << "# The seed value is " << seed_v << endl;
    if(!mcobj.isrestart()) diag_wHeader(bendobj, stretchobj, chargeobj, mesh, fileptr);
    if(mcobj.isrestart()) fileptr << "# Restart index " << residx << endl;

    clock_t timer;
    Etot = mcobj.evalEnergy(mesh);
    mcobj.write_energy(fileptr, iter, mesh);
    for(iter=residx; iter < mcobj.totaliter(); iter++){
        if(iter%mcobj.dumpskip() == 0){
            outfile=outfolder+"/snap_"+ZeroPadNumber(iter/mcobj.dumpskip())+".h5";
            hdf5_io_delete(outfile);
            hdf5_io_write((double*) mesh.pos, 3*mesh.N, outfile, "pos");
            hdf5_io_write_mesh(mesh.numnbr, mesh.node_nbr_list, mesh.N,
                                mesh.nghst, outfile);
            
            if (mesh.ncomp>1) hdf5_io_write(mesh.compA, mesh.N, outfile, "lip");
            fstream restartfile(outfolder+"/restartindex.txt", ios::out);
            restartfile << iter << " " << mcobj.dumpskip() << endl;
            restartfile.close();
        }
        //
        if(repulsiveobj.isSelfRepulsive()) repulsiveobj.buildCellList(mesh);
        num_moves = mcobj.monte_carlo_3d(mesh.pos, mesh);
        if (mesh.ncomp>1) num_exchange = mcobj.monte_carlo_lipid(mesh.pos, mesh);
        if (mcobj.isfluid() && !(iter % mcobj.fluidizeevery())){
            num_bond_change = mcobj.monte_carlo_fluid(mesh.pos, mesh);
            (*terminal) << "fluid stats " << num_bond_change << " bonds flipped" << endl;
        }
        //
        if(!(iter % recaliter)){
            Etot = mcobj.evalEnergy(mesh);
            (*terminal) << "iter = " << iter << 
            "; Accepted Moves = " << (double)num_moves*100/mcobj.onemciter() 
            << " %;"
            "; Exchanged Moves = " << (double)num_exchange * 100 / mcobj.onemciter()
            << " %;"
            << " totalener = " << Etot << "; volume = " << mcobj.getvolume() << endl;
        }
        mcobj.write_energy(fileptr, iter, mesh);
    }
    (*terminal) << "Total time taken = " << (clock()-timer)/CLOCKS_PER_SEC << "s" << endl;
    if (world_size > 1) out_file.close();
    
    fileptr.close();
    mesh.free();
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_err = MPI_Finalize();
    return 0;
}