#include "activity.hpp"
#include <fstream>

extern "C"  void  Activity_listread(char *, bool *, double *, double *, char *);


 
int ACT::initACT(int N, std::string fname){
  char tmp_fname1[128], tmp_fname2[128];
  string parafile, outfile;
  string whichactivity;
  double minA, maxA;
  int i;

  parafile = fname+"/para_file.in";
  sprintf(tmp_fname1, "%s", parafile.c_str() );
  Activity_listread(tmp_fname2, &doactivity, &maxA, &minA, tmp_fname1);

  whichactivity = tmp_fname2;
  ofstream out_;

  out_.open( fname+"/activitypara.out");
  out_<< "# =========== activity parameters ==========" << endl
      << " N " << N << endl
      << " doactivity = " << doactivity << endl
      << " whichactivity = " << whichactivity << endl
      << " minA " << minA << endl
      << " maxA " << maxA << endl;
  out_.close();
 
    for(i = 0; i < N; i++){
      Actarr.push_back(0.0);
    }



  if(whichactivity == "random"){
        for(i=0;i<N;i++) Actarr[i] = RandomGenerator::generateUniform(minA, maxA);
    }

    if(whichactivity == "constant"){
        for(i=0;i<N;i++) Actarr[i] = maxA;
    }

  return 0; 
}

bool ACT::is_active(){ return doactivity;}

double ACT::getActivityIdx(int idx){return Actarr[idx];}
