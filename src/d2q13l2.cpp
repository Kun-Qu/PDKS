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

#include "d2q13l2.h"


double e_level[2], lbm_omega, Gamma, Gamma_1;



void macro(const double c_f[26], 
           double& rr, double& ru, double& rv, double& re)
{
    rr = 0.0E0;
    ru = 0.0E0;
    rv = 0.0E0;
    re = 0.0E0;
    
    double mi;
    for(int idir=0;idir<26;idir++)
    {
        mi =  c_f[idir];
        rr += mi;
        ru += mi*dc[idir][0];
        rv += mi*dc[idir][1];
        re += 5.0E-1 * mi*dc2[idir];
        if(idir%2==1) re += mi * e_level[1]; 
    }
    ru /= rr;
    rv /= rr;
    re /= rr;
    re -= 5.0E-1 * (ru*ru + rv*rv); 
    
    
    return;
}


void mirror( const double cnx, const double cny,   // inside unit normal vector of the mirror
             const double cf[26],              // the original f
             double rf[26]                     // the f for output
           )
{
    for(int idir=0;idir<26;idir++) 
        rf[idir] = 0.0E0;
    
    
    double rx,ry,weight[26];
    for(int dir=0;dir<26;dir++)
    {
        // compute reflected location (rx,ry) of node[dir]
        reflect(cnx,cny, dc[dir][0],dc[dir][1], rx,ry);
        //std::cout<<dir<<std::endl;
        double rz;
        rz = dir%2;
        get_weight(rx,ry,rz,weight);
       
        // assign f(rx,ry) onto new nodes
        for(int idir=0;idir<26;idir++)
            rf[idir] += weight[idir]*cf[dir]; 
    }
       

    
    return;
}




void collision(double ruve[4],
               double c_bf[26],
               double c_df[26],
               const double dt,
               const double c_beta)
{
    double c_fq[26], tdf;
    
    get_feq(ruve[0],ruve[1],ruve[2],ruve[3],c_fq);
    
    for(int dir=0;dir<26;dir++)
    {
        tdf = ( c_fq[dir] - c_bf[dir] )*lbm_omega*dt*c_beta;
        c_bf[dir] += tdf;
        c_df[dir] += tdf;
    }
}



void get_weight(const double cx, const double cy, const double cz,double weight[26])
{
    double t1 ,t3 ,t5 ,t7 ,t8 ,t11,t12,t13,t14,t15,t16,t17,t18,t20,t21,t22,
           t23,t24,t25,t26,t30,t35,t36,t37,t38,t40,t41,t42,t43;
    double xyw[13];
    
      t1  = cx*cx;
      t3  = cy*cy;
      t5  = t1*t1;
      t7  = t1*t3;
      t8  = t3*t3;
      t11 = 2.0/3.0*cx;
      t12 = 2.0/3.0*t1;
      t13 = t1*cx;
      t14 = t13/6.0;
      t15 = cx*t3;
      t16 = t15/2.0;
      t17 = t5/6.0;
      t18 = t7/2.0;
      t20 = 2.0/3.0*cy;
      t21 = 2.0/3.0*t3;
      t22 = t1*cy;
      t23 = t22/2.0;
      t24 = t3*cy;
      t25 = t24/6.0;
      t26 = t8/6.0;
      t30 = cy*cx;
      t35 = cx/12.0;
      t36 = t1/24.0;
      t37 = t13/12.0;
      t38 = t5/24.0;
      t40 = cy/12.0;
      t41 = t3/24.0;
      t42 = t24/12.0;
      t43 = t8/24.0;
      
      xyw[0] = 1.0-5.0/4.0*t1-5.0/4.0*t3+t5/4.0+t7+t8/4.0;
      xyw[1] = t11+t12-t14-t16-t17-t18;
      xyw[2] = t20+t21-t23-t25-t18-t26;
      xyw[3] = -t11+t12+t14+t16-t17-t18;
      xyw[4] = -t20+t21+t23+t25-t18-t26;
      xyw[5] = t30/4.0+t22/4.0+t15/4.0+t7/4.0;
      xyw[6] = -t30/4.0+t22/4.0-t15/4.0+t7/4.0;
      xyw[7] = t30/4.0-t22/4.0-t15/4.0+t7/4.0;
      xyw[8] = -t30/4.0-t22/4.0+t15/4.0+t7/4.0;
      xyw[9] = -t35-t36+t37+t38;
      xyw[10] = -t40-t41+t42+t43;
      xyw[11] = t35-t36-t37+t38;
      xyw[12] = t40-t41-t42+t43;
      
      for(int dir=0;dir<13;dir++)
      {
        weight[dir+dir]     = xyw[dir] * (1.0E0 - cz);
        weight[dir+dir+1]  =  xyw[dir] * cz;
      }

    return;
}




/////////////////////////////////////////////////////////////////

void
get_feq (const double rho, const double u, const double v, const double e,
	  double feq[26])
{
  double c = (Gamma - 1) * e;
  double ep = e - c;
  c = sqrt (2.0E0 * c);

  double t1 = v * v;
  double t2 = t1 * t1;
  double t3 = 4.0E0 * t2;
  double t4 = c * c;
  double t6 = u * u;
  double t7 = t6 * t1;
  double t9 = t4 * t4;
  double t11 = t1 * t4;
  double t13 = t6 * t4;
  double t17 = t6 * t6;
  double t18 = 4.0E0 * t17;
  double t20 =
    rho * (t3 - 20.0E0 * t4 + 16.0E0 * t7 + 16.0E0 + 5.0E0 * t9 +
	   20.0E0 * t11 + 20.0E0 * t13 - 20.0E0 * t1 - 20.0E0 * t6 + t18);
  double t21 = e_level[1] - ep;
  double t23 = 1.E0 / (e_level[1] - e_level[0]);
  double t24 = t21 * t23;
  double t27 = ep - e_level[0];
  double t28 = t27 * t23;
  double t31 = 18.0E0 * t13;
  double t32 = 12.0E0 * t7;
  double t33 = t4 * u;
  double t34 = 12.0E0 * t33;
  double t35 = 6.0E0 * t11;
  double t36 = 16.0E0 * u;
  double t37 = t6 * u;
  double t38 = 4.0E0 * t37;
  double t39 = u * t1;
  double t40 = 12.0E0 * t39;
  double t41 = 3.0E0 * t9;
  double t42 = 16.0E0 * t6;
  double t43 = 8.0E0 * t4;
  double t44 =
    t31 + t32 + t34 + t35 - t36 + t38 + t40 + t41 + t18 - t42 - t43;
  double t45 = rho * t44;
  double t50 = 16.0E0 * v;
  double t51 = t6 * v;
  double t52 = 12.0E0 * t51;
  double t53 = 6.0E0 * t13;
  double t54 = 18.0E0 * t11;
  double t55 = t4 * v;
  double t56 = 12.0E0 * t55;
  double t57 = t1 * v;
  double t58 = 4.0E0 * t57;
  double t59 = 16.0E0 * t1;
  double t60 =
    -t50 - t43 + t52 + t32 + t53 + t54 + t56 + t58 + t41 - t59 + t3;
  double t61 = rho * t60;
  double t66 =
    t31 + t32 - t34 + t35 + t36 - t38 - t40 + t41 + t18 - t42 - t43;
  double t67 = rho * t66;
  double t72 =
    -t52 - t58 + t50 - t56 + t32 + t41 + t3 - t59 + t53 + t54 - t43;
  double t73 = rho * t72;
  double t78 = 8.0E0 * t51;
  double t79 = 4.0E0 * t13;
  double t80 = 4.0E0 * t55;
  double t82 = 8.0E0 * u * v;
  double t83 = 8.0E0 * t7;
  double t84 = 8.0E0 * t39;
  double t85 = 4.0E0 * t11;
  double t86 = 4.0E0 * t33;
  double t88 = rho * (t78 + t79 + t80 + t82 + t83 + t9 + t84 + t85 + t86);
  double t94 = rho * (t78 + t79 + t80 - t82 + t83 + t9 - t84 + t85 - t86);
  double t100 = rho * (-t78 + t79 - t80 + t82 + t83 + t9 - t84 + t85 - t86);
  double t106 = rho * (-t78 + t79 - t80 - t82 + t83 + t9 + t84 + t85 + t86);
  double t112 = rho * u / 12.0E0;
  double t114 = rho * t17 / 24.0E0;
  double t115 = rho * t4;
  double t116 = t115 / 48.0E0;
  double t117 = rho * t6;
  double t118 = t117 / 24.0E0;
  double t120 = t117 * t4 / 8.0E0;
  double t122 = rho * t9 / 64.0E0;
  double t124 = rho * t37 / 12.0E0;
  double t126 = t115 * u / 8.0E0;
  double t127 = -t112 + t114 - t116 - t118 + t120 + t122 + t124 + t126;
  double t132 = 24.0E0 * t55;
  double t133 = 8.0E0 * t1;
  double t134 = 8.0E0 * t2;
  double t135 = 4.0E0 * t4;
  double t136 = 24.0E0 * t11;
  double t137 = 16.0E0 * t57;
  double t139 = rho * (-t50 + t132 - t133 + t134 - t135 + t136 + t137 + t41);
  double t144 = t112 + t114 - t116 - t118 + t120 + t122 - t124 - t126;
  double t150 = rho * (t50 - t132 - t133 + t134 - t135 + t136 - t137 + t41);

  feq[0] = t20 * t24 / 16.0E0;
  feq[1] = t20 * t28 / 16.0E0;
  
  feq[2] = -t45 * t24 / 24.0E0;
  feq[3] = -t45 * t28 / 24.0E0;
  
  feq[4] = -t61 * t24 / 24.0E0;
  feq[5] = -t61 * t28 / 24.0E0;
  
  feq[6] = -t67 * t24 / 24.0E0;
  feq[7] = -t67 * t28 / 24.0E0;
  
  feq[8] = -t73 * t24 / 24.0E0;
  feq[9] = -t73 * t28 / 24.0E0;
  
  feq[10] = t88 * t24 / 32.0E0;
  feq[11] = t88 * t28 / 32.0E0;
  
  feq[12] = t94 * t24 / 32.0E0;
  feq[13] = t94 * t28 / 32.0E0;
  
  feq[14] = t100 * t24 / 32.0E0;
  feq[15] = t100 * t28 / 32.0E0;
  
  feq[16] = t106 * t24 / 32.0E0;
  feq[17] = t106 * t28 / 32.0E0;
  
  feq[18] = t127 * t21 * t23;
  feq[19] = t127 * t27 * t23;
  
  feq[20] = t139 * t24 / 192.0E0;
  feq[21] = t139 * t28 / 192.0E0;
  
  feq[22] = t144 * t21 * t23;
  feq[23] = t144 * t27 * t23;
  
  feq[24] = t150 * t24 / 192.0E0;
  feq[25] = t150 * t28 / 192.0E0;
  
  return;
}

void get_vstress(double ruve[4], double f[26], double stress[3])
{
    stress[0]=0.0E0;
    stress[1]=0.0E0;
    stress[2]=0.0E0;
    
    return;
    
}
