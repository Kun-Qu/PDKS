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


#ifndef LBM_H
#define LBM_H
#include <cmath>
#include <iostream>

extern double ref_e,ref_U,ref_rho,ref_L, ref_p;     
extern double Re_L,Re_rho,Re_U;   
extern double e_level[2], lbm_omega, Gamma, Gamma_1;


inline void reflect(const double cnx,const double cny,
                    const double ox, const double oy,
                    double& rx,      double& ry,
                    const int mode=0)
{
    double on,onx,ony;
    if(mode==0)  // specular reflection
    {
        on = ox*cnx + oy*cny;
        onx = cnx*on;
        ony = cny*on;
        onx += onx;
        ony += ony;
        rx = ox - onx;
        ry = oy - ony;
    }
    else if(mode==1)  // bounce back reflection
    {
        rx = -ox;
        ry = -oy;
    }
    return;
}

inline double get_pressure(const double crho, const double ce)
{
    return Gamma_1 * crho * ce;
}

inline double get_omega(void) {return lbm_omega;}


inline void set_ref(const double e0, const double rho0, const double L0)
{
    ref_e = e0;
    ref_U = sqrt(ref_e);
    ref_rho = rho0;
    ref_p = ref_rho * ref_e;
    ref_L = L0;
}



inline void get_ref(double& vr, double& vU, double& ve, double& vL)
{
    vr = ref_rho;
    vU = ref_U;
    ve = ref_e;
    vL = ref_L;
    
    return;
}

inline void get_ReV(double& Vr, double& VL, double& VU)
{
	Vr = Re_rho;
	VL = Re_L;
	VU = Re_U;
	
	return;
	
}


#ifdef INVISCID
    #include "d2q13l2.h"
#endif

#ifdef LAMINAR
    #include "d2q17l2.h"
#endif 


#endif


