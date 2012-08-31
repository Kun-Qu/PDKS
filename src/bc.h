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

#ifndef BC_H
#define BC_H


#include "dks.h"
#include "lbm.h"
#include "grid_2d.h"

#include <iostream>
#include <string>


const int ITAN  = 101,    // inviscid wall
          IWALL = 102,   // viscous wall
          IFAR  = 103,    // farfield
          ITAN2 = 104,
          // 1xx  are face center boundaries
          
          ISPIN  = 201,   // supersonic inflow
          ISUIN  = 202,   // subsonic inflow
          ISPOUT = 203,  // supersonic outflow
          ISUOUT = 204,  // subsonic outflow
          ISYM   = 205,  // symmetric boundary
          // 2xx are ghost cell center boundaroes
          
          IPERD  = 301,    // period boudanry
          ICGRID = 302,    // for C-type grid
          // 3xx are cut BC
          
          IUSRDEF = 401;   // user defined BC
          
          
          
const int MAXBT=12;
const int MAXBC=50;


extern int nbcs;
/**
bcdat[ibc][0]: BC type
          [1]: BC direction, 0=i, 1=j
          [2]: i(1 or IDIM)
          [3]: istrt
          [4]: iend
          [5]: jstrt
          [6]: jend
          [7]: profile index
*/
extern int bcdat[MAXBC][8];
extern double bcm[MAXBC][4];  // boundary value in rho, u, v, e


void read_bc(const int n_bcs, const std::string& bc_f);
void bc(double time);
void get_bc_macro(double east[JDIM+3][NVAR],
                  double west[JDIM+3][NVAR],
                  double south[IDIM+3][NVAR],
                  double north[IDIM+3][NVAR]);
int setup_cut(const int ibc);


int bc_cgrid(const int ibc);
int bc_tan(const int ibc );
int bc_tan2(const int ibc );
int bc_wall(const int ibc );
int bc_spout(const int ibc );
int bc_spin(const int ibc );
int bc_suin(const int ibc);
int bc_suout(const int ibc);
int bc_far(const int ibc );
int bc_perd(const int ibc);
int bc_usrdef(const int ibc, double time);

void reimann_farfield(double rho[3], double u[3], double v[3], double e[3], 
                double nx, double ny);
                
void metbc(void);
int set_bcdat(const int ibc, const int ibcpar[8], const double dbcpar[4]);

void get_anm(const double ct, const double cx, const double cy, double* cruve);
#endif
