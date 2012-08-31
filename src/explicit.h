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


#ifndef EXPLICIT_H
#define EXPLICIT_H

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "dks.h"
#include "grid_2d.h"
#include "bc.h"
#include "parallel.h"

inline double sign(const double var)
{
    return var>0.0E0?1.0E0:-1.0E0;
}

inline double smooth(const double v1, const double v2)
{
const double epsq = 1.0E-12, eps=1.0E-6;
return (2.0E0*v1*v2)/(v1*v1+v2*v2+epsq);
}



extern int flag[IDIM+3][JDIM+3];

extern double DQP[MAXDIM+3],DQM[MAXDIM+3],QP[MAXDIM+3],QM[MAXDIM+3];

void clear_space(void);

// MUSCL scheme
void limiter(  const double rkap,
               const double vm, const double vp,
               double &dp, double &dm, const int index);

void muscli(const int cj, const int cdir,
           const int bcdat[][8],
           const int nbcs
           );

void musclj(const int ci, const int cdir,
           const int bcdat[][8],
           const int nbcs
           );

void fluxi(const int cj, const int cdir);
void fluxj(const int ci, const int cdir);

void convect(const double dt, const double c_alpha);

//
void step(const double dt, const double c_alpha, const double c_beta);
void source(const double dt, const double c_beta);
void update(void);

void usr_init(void);
void init(const double c_rho, 
          const double c_u,
          const double c_v,
          const double c_e);
int read_init(const std::string& init_f);

double get_inf_residual(const double dt);
double get_nsq_residual(const double dt);
double get_inf_residual_P(const double dt);
double get_nsq_residual_P(const double dt);

// IO
int write_plt(std::string& prefix, const double c_time,
              const int c_count, const int field);
// write or read a checkpoint file (restart file)
int wr_check_2_0(std::string& prefix_or_filename,  
                      double& c_time,
                      int& c_count,
                      const int write_or_read
                      );
                      
int compute_force(std::string& prefix, 
	                const int ibc, 
                  const double c_time, 
                  const int c_count,
                  const int write_file,
                  double* Cx,
                  double* Cy);

int get_global_nfo(std::string& dataname, 
	                 std::ofstream& logfile,
	                 double& t, 
	                 const double dt,
	                 int& iter, 
	                 const int save_int, 
	                 const int check_int, 
	                 const int displ_int,
                         const char anyway);

void get_anM(const double ct, const double cx, const double cy, 
             double* cruve);
#endif
