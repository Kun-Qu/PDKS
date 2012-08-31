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

#ifndef D2Q17L2_H
#define D2Q17L2_H

#include <iostream>
#include <cmath>

const double dc[34][2] = 
{
//  0
{0.0E0, 0.0E0},
{0.0E0, 0.0E0},

// 1..4
{2.0E0/3.0E0,  2.0E0/3.0E0},
{2.0E0/3.0E0,  2.0E0/3.0E0},

{-2.0E0/3.0E0,  2.0E0/3.0E0},
{-2.0E0/3.0E0,  2.0E0/3.0E0},

{-2.0E0/3.0E0,  -2.0E0/3.0E0},
{-2.0E0/3.0E0,  -2.0E0/3.0E0},

{2.0E0/3.0E0,  -2.0E0/3.0E0},
{2.0E0/3.0E0,  -2.0E0/3.0E0},

// 5..6
{2.0E0/3.0E0, 4.0E0/3.0E0},
{2.0E0/3.0E0, 4.0E0/3.0E0},

{4.0E0/3.0E0, 2.0E0/3.0E0},
{4.0E0/3.0E0, 2.0E0/3.0E0},

// 7..8
{-2.0E0/3.0E0, 4.0E0/3.0E0},
{-2.0E0/3.0E0, 4.0E0/3.0E0},

{-4.0E0/3.0E0, 2.0E0/3.0E0},
{-4.0E0/3.0E0, 2.0E0/3.0E0},

// 9..10
{-2.0E0/3.0E0, -4.0E0/3.0E0},
{-2.0E0/3.0E0, -4.0E0/3.0E0},

{-4.0E0/3.0E0, -2.0E0/3.0E0},
{-4.0E0/3.0E0, -2.0E0/3.0E0},

// 11..12
{2.0E0/3.0E0, -4.0E0/3.0E0},
{2.0E0/3.0E0, -4.0E0/3.0E0},

{4.0E0/3.0E0, -2.0E0/3.0E0},
{4.0E0/3.0E0, -2.0E0/3.0E0},

// 13..17
{2.0E0, 0.0E0},
{2.0E0, 0.0E0},

{0.0E0, 2.0E0},
{0.0E0, 2.0E0},

{-2.0E0, 0.0E0},
{-2.0E0, 0.0E0},

{0.0E0, -2.0E0},
{0.0E0, -2.0E0}


};


// (c_i*c_i)/2
const double ei[17] =
{
0.0E0,

4.0E0/9.0E0,4.0E0/9.0E0,4.0E0/9.0E0,4.0E0/9.0E0,

1.0E1/9.0E0,1.0E1/9.0E0,1.0E1/9.0E0,1.0E1/9.0E0,
1.0E1/9.0E0,1.0E1/9.0E0,1.0E1/9.0E0,1.0E1/9.0E0,

2.0E0,2.0E0,2.0E0,2.0E0
};




// Prandtl number, theta (for adjusting Pr), Mu of molecule, heat conduction
extern double Pr, theta, Mum;  

inline double get_tau(const double c_p)
{
    return Mum/c_p;
};

void get_feq (const double rho, const double u, const double v, const double e,
              double feq[34]
             );
	 
void mirror( const double cnx, const double cny,   // inside normal vector of the mirror
             const double cf[34],              // the original f
             double rf[34],                     // the f for output
             const int mode=0
           );

void get_weight(const double cx, const double cy, double weight[17]);

void collision(double* MV,
               double c_bf[34],
               double c_df[34],
               const double dt,
               const double c_beta);
               
void macro(const double c_f[34], 
           double& rr, double& ru, double& rv, double& re);

inline void set_lbm(
                const double pRe_Rho,
                const double pRe_U,
                const double pRe_L,
                const double e0, 
                const double c_Re, 
                const double cPr=0.72E0,    // for air
                const double Gamma0=1.4E0  // for air
                )
{
    e_level[0] = 0.0E0;
    e_level[1] = 1.0E0;
    
#ifdef LAMINAR
    Re_L = pRe_L;
    Re_rho = pRe_Rho;
    Re_U = pRe_U;
    
    //viscosity, is used to compute Tau
    Mum = (Re_U/ref_U)*(Re_rho/ref_rho)*(Re_L/ref_L)/c_Re;
    double Rtau = Mum/(ref_rho*ref_U*ref_U);
    //std::cout<<"Dimesionless vsicosity is "<<Mum<<" Characteristic Tau is "<<Rtau<<std::endl;
#endif

#ifdef INVISCID
   lbm_omega = 1.0E0/c_Re;
   //std::cout<<"Dimesionless tau is "<< c_Re <<std::endl;
#endif

    Pr=cPr;
    
    Gamma = Gamma0;
    Gamma_1 = Gamma - 1.0E0;
            
    theta = Gamma/Pr;
    
    return;
}
void get_vstress(double ruve[4], double f[34], double stress[3]);
#endif
