/*
PDKS是屈崑在新加坡国立大学攻读博士学位期间研究可压缩
LBM模型时开发的程序LBM可压缩流模拟程序。基本特点为：
1 支持二维单块贴体网格
2 支持多个可压缩LBM模型
3 支持I方向网格分块MPI并行
4 支持http://www.math.tu-cottbus.de/~berti/gral中的control-device
5 我忘了……

具体理论部分请参考我的论文
“Alternative method to construct equilibrium distribution functions in 
lattice-Boltzmann method simulation of inviscid compressible flows at 
high Mach number”, Physical Review E, 2007, Vol 75 (3), 036706, Kun Qu, 
Chang Shu and Yong Tian Chew
http://serve.me.nus.edu.sg/shuchang/Publications/Latest%20papers%20for%20web/Qu_Shu_Chew%20(PRE)%202007.pdf

和学位论文
https://scholarbank.nus.edu.sg/handle/10635/15823



PDKS (Parallel Discrete Kinetic Solver), written by Qu Kun in 2006, is a " 
lattice Boltzmann equation solver for two dimentional compressible inviscid/viscous 
flows.
Copyright (C) <2007>  <Kun Qu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include <iostream> 
#include <string> 
#include <sstream> 
 
 
#include "IO/mutator.h" 
#include "IO/control-device.h" 
#include "dks.h" 
#include "explicit.h" 
#include "parallel.h" 
 
double Q[IDIM+3][JDIM+3][DIR];  // current values at step n 
double DQ[IDIM+3][JDIM+3][DIR]; // new values at step n+1 
double M[IDIM+3][JDIM+3][NVAR];       // rho, u, v, e; 
 
int read_par(std::string& case_f, std::string& case_name, 
              double& dt, int& max_step, int& print_step, int& save_step, 
              int& check_step, int& c_count, double& c_time); 
 
int read_case(std::string& case_f, 
              std::string& case_name, 
              double& dt, 
              int& max_step, 
              int& print_step,  
              int& save_step, 
              int& check_step, 
              int& c_count, 
              double& c_time); 
 
 
 
 
int main(int argc, char* argv[]) 
{ 
///////////////////////////////////////////////// 
   
 
  /* Initialize MPI */ 
  MPI_Init(&argc, &argv); 
  /* Find out my identity in the default communicator */ 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  MPI_Comm_size(MPI_COMM_WORLD, &cpusize); 
  if(myrank==0) std::cout<<"I have "<<cpusize<<" brothers"<<std::endl;
  
  if(!partition()) 
  {
     
     std::cout<<"Partition Error!"<<std::endl; 
     MPI_Finalize();
     return 0;
  }
////////////////////////////////////////////////   
 
	 
    double t=0.0E0, dt; 
    int maxstep, display, save,check; 
     
    if(myrank==0) 
    std::cout<<"PDKS (Parallel Discrete Kinetic Solver), written by Qu Kun in 2006, is a " 
             <<"lattice Boltzmann equation solver for two dimentional compressible" 
             <<"inviscid/viscous flows.\n" 
             <<"Usage: Pdks casefile logfile\n\n" 
             <<"Note NI and NJ are "<<IDIM<<" and "<<JDIM<<std::endl; 
     
    if(argv[1]==NULL && myrank==0)  
    { 
        std::cerr<<"You should specify the name of the case file."<<std::endl; 
        return 0; 
    } 
     
    std::string casefile(argv[1]); 
    std::string logname(argv[2]); 
    //std::cout<<casefile<<"\n\n\n"<<logname<<std::endl;     
        	 
  
    std::ofstream logfile;
    if(myrank == 0)
         logfile.open(logname.c_str(), ios::app); 
    std::string dataname(""); 
    int count=1; 
    read_par(casefile, dataname, dt, maxstep, display, save, check, count, t); 
    
    double time_st = MPI_Wtime();
    double rr,rU,re,rL,rt;
    get_ref(rr,rU,re,rL);
    rt = rL/rU;

    for(;count<=maxstep;count++) 
    { 
        step(dt,1.0E0, 1.0E0);  
        synchronize(0); 
        t += dt; 
        bc(t); 
        
        get_global_nfo(dataname, logfile, t, dt, count, save, check, display, 0); 
        if(count%save==0)
      	{        		
            if(myrank==0)  
            std::cout<<t*rt<<" Speed = "<<(MPI_Wtime()-time_st)/save
                     <<" Second/Step"<<std::endl;
       	    time_st = MPI_Wtime();
       	}
    }
    
    //wr_check_2_0(dataname,t,maxstep,0); 
    if(maxstep%save)  
        get_global_nfo(dataname, logfile, t, dt, count, save, check, display, 2);
     
////////////////////////////////////// 
  MPI_Finalize(); 
///////////////////////////////////// 
    return 1;    
} 
 
 
 
/* 
Read parameter from a casefile 
*/ 
int read_par(std::string& case_f, std::string& case_name, 
              double& dt, int& max_step, int& print_step, int& save_step, 
              int& check_step, int& c_count, double& c_time) 
{ 
   ControlDevice Ctrl( GetFileControlDevice(case_f.c_str(), "DKS") ); 
 
  string grid_name; 
  Ctrl.add("case_name", GetMutator(case_name)); 
  Ctrl.add("grid_name", GetMutator(grid_name)); 
 
  int vis_mode; 
  double tgamma,tRe,tPr,rfrho,rfe,rfL,Re_Rho,Re_U,Re_L; 
  Ctrl.add("ref_rho", GetMutator(rfrho)); 
  Ctrl.add("ref_e", GetMutator(rfe)); 
  Ctrl.add("ref_L", GetMutator(rfL)); 
  Ctrl.add("Re_Rho",GetMutator(Re_Rho)); 
  Ctrl.add("Re_U",GetMutator(Re_U)); 
  Ctrl.add("Re_L",GetMutator(Re_L)); 
 
   
 
  ControlDevice Ctrl1 = Ctrl.getSubDevice("lbm_par"); 
  Ctrl1.add("vis_mode", GetMutator(vis_mode)); 
  Ctrl1.add("Gamma",  GetMutator(tgamma)); 
  Ctrl1.add("Re_tau", GetMutator(tRe)); 
  Ctrl1.add("Pr", GetMutator(tPr)); 
 
  int t_nbc; 
  std::string bc_file; 
  Ctrl.add("nbc", GetMutator(t_nbc)); 
  Ctrl.add("bc_file", GetMutator(bc_file)); 
 
  ControlDevice Ctrl3 = Ctrl.getSubDevice("init_condition"); 
  int init_mode; 
  double tr0,tu0,tv0,te0; 
  std::string init_file; 
  Ctrl3.add("init_mode", GetMutator(init_mode)); 
  Ctrl3.add("init_file", GetMutator(init_file)); 
  Ctrl3.add("irho", GetMutator(tr0)); 
  Ctrl3.add("iu", GetMutator(tu0)); 
  Ctrl3.add("iv", GetMutator(tv0)); 
  Ctrl3.add("ie", GetMutator(te0)); 
 
  ControlDevice Ctrl4 = Ctrl.getSubDevice("time_control"); 
  Ctrl4.add("time_step", GetMutator(dt)); 
  Ctrl4.add("max_step", GetMutator(max_step)); 
  Ctrl4.add("print_step", GetMutator(print_step)); 
  Ctrl4.add("save_step", GetMutator(save_step)); 
  Ctrl4.add("check_step", GetMutator(check_step)); 
 
  Ctrl.update(); 
  //Ctrl.print_values(cout); 
 
  // Grid 
  read_grid(grid_name.c_str()); 
  extend_grid(); 
  metric(); 
  cell_center(); 
  //std::cout<<myrank<<" Grid done"<<std::endl; 
 
  // LBM model 
  /* 
  if(vis_mode==0) 
      std::cout<<"Inviscid flow"<<std::endl; 
  else if(vis_mode>=1) 
      std::cout<<"Laminar flow"<<std::endl; 
  */ 
   
  set_ref(rfe, rfrho, rfL); 
/*   
  #ifdef D2Q15L2 
  set_lbm(Re_Rho, Re_U, Re_L,  rfe, tRe, tPr, tgamma); 
  #endif 
*/
  #ifdef D2Q17L2 
  set_lbm(Re_Rho, Re_U, Re_L,  rfe, tRe, tPr, tgamma); 
  #endif   
  #ifdef D2Q13L2 
  set_lbm(rfe, tRe, 0.0, tgamma); 
  #endif 
 
  // init 
  if(init_mode==0) // using irho,iu,iv,ie as initial values 
        init(tr0,tu0,tv0,te0); 
  else if(init_mode==1)  // reading rho/u/v/e from a file 
        usr_init();// 
  else if(init_mode==2) 
  	{
  		std::cout<<myrank<<" CPU Reading Checkpoint file as initial condition"<<std::endl;
        wr_check_2_0(init_file, c_time, c_count,1 );// 
    }
  else 
      { 
            std::cerr<<"Wrong init mode"<<std::endl; 
            return -1; 
      }   
 
  // read bc data 
  read_bc(t_nbc, bc_file); 
  bc(0.0E0); 
 
  return 1; 
} 
 
