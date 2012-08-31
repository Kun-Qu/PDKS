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

#include "d2q15l2.h"
double e_level[2], Gamma, lbm_omega, Gamma_1,
        Pr,theta, Mum, Kam, Rgas, Cv, Cp;



void macro(const double c_f[30], 
           double& rr, double& ru, double& rv, double& re)
{
    rr = 0.0E0;
    ru = 0.0E0;
    rv = 0.0E0;
    re = 0.0E0;
    
    double mi;
    for(int idir=0;idir<15;idir++)
    {
        mi =  c_f[idir+idir+1] + c_f[idir+idir];
        rr += mi;
        ru += mi*dc[idir+idir][0];
        rv += mi*dc[idir+idir][1];
        re +=  mi * ei[idir] + c_f[idir+idir+1];
    }
    
    ru /= rr;
    rv /= rr;
    re /= rr;
    re -= 5.0E-1 * (ru*ru + rv*rv); 
    
    
    return;
}


void mirror( const double cnx, const double cny,   // inside unit normal vector of the mirror
             const double cf[30],              // the original f
             double rf[30],                  // the f for output
             const int mode
           )
{
    for(int idir=0;idir<30;idir++) 
        rf[idir] = 0.0E0;
    
    
    double rx,ry,weight[15];
    for(int dir=0;dir<15;dir++)
    {
        // compute reflected location (rx,ry) of node[dir]
        reflect(cnx,cny, dc[dir+dir][0],dc[dir+dir][1],  rx,ry, mode);
        get_weight(rx,ry,weight);       
        // assign f(rx,ry) onto new nodes
        for(int idir=0;idir<15;idir++)
        {
            rf[idir+idir] += weight[idir]*cf[dir+dir]; 
            rf[idir+idir+1] += weight[idir]*cf[dir+dir+1]; 
        }
    }

    return;
}



////////////////////////////////////////////////////////////////////////////////
void collision(double* MV,           
               double c_bf[30],
               double c_df[30],
               const double dt,
               const double c_beta)
{
    double c_fq[30], tdf;
    
    get_feq(MV[0],MV[1],MV[2],MV[3],c_fq);

#ifdef LAMINAR
    // compute Tau according to pressure and molecule viscosity
    lbm_omega = MV[0]*MV[3]*Gamma_1 / Mum;
    MV[4]=lbm_omega*dt; 
#endif

#ifdef INVISCID
    
#endif

    for(int dir=0;dir<30;dir++)
    {
        tdf = ( c_fq[dir] - c_bf[dir] )*lbm_omega*dt*c_beta;
        c_bf[dir] += tdf;
        c_df[dir] += tdf;
    }
}


void get_weight(const double x, const double y, 
        double weight[15])
{
double      t1 = sqrt(3.0E0);
double      t2 = t1*y;
double      t4 = y*y;
double      t6 = t4*y;
double      t7 = t1*t6;
double      t9 = t4*t4;
double      t12 = 4.0E0/81.0E0*x;
double      t13 = 20.0E0/243.0E0*t2;
double      t14 = t1*x;
double      t15 = t14*y;
double      t16 = 2.0E0/27.0E0*t15;
double      t17 = 10.0E0/27.0E0*t4;
double      t18 = x*t4;
double      t19 = 2.0E0/9.0E0*t18;
double      t20 = 10.0E0/81.0E0*t7;
double      t21 = t14*t6;
double      t22 = 4.0E0/27.0E0*t21;
double      t23 = 4.0E0/27.0E0*t9;
double      t26 = 5.0E0/27.0E0*x;
double      t27 = 5.0E0/81.0E0*t2;
double      t28 = x*x;
double      t29 = t28/9.0E0;
double      t30 = 7.0E0/18.0E0*t15;
double      t31 = 5.0E0/18.0E0*t4;
double      t33 = t1*t28*y;
double      t34 = 5.0E0/18.0E0*t33;
double      t35 = 5.0E0/54.0E0*t7;
double      t36 = t28*t4;
double      t37 = t36/3.0E0;
double      t38 = 2.0E0/9.0E0*t21;
double      t39 = t9/9.0E0;
double      t40 = 4.0E0/81.0E0-t26+t27+t29-t30-t31+t34-t35+t37+t38+t39;
double      t42 = 2.0E0/9.0E0*t28;
double      t43 = 2.0E0/9.0E0*t4;
double      t49 = 4.0E0/81.0E0+t26+t27+t29+t30-t31+t34-t35+t37-t38+t39;
double      t50 = 4.0E0/27.0E0*x;
double      t51 = 4.0E0/243.0E0*t2;
double      t52 = 2.0E0/9.0E0*t15;
double      t53 = 4.0E0/27.0E0*t4;
double      t54 = t28*x;
double      t55 = 2.0E0/9.0E0*t54;
double      t56 = t33/9.0E0;
double      t57 = t7/81.0E0;
double      t59 = t1*t54*y;
double      t60 = t59/9.0E0;
double      t61 = t21/9.0E0;
double      t62 = t9/27.0E0;
double      t63 = -16.0E0/243.0E0+t50+t51+t42+t52+t53-t55-t56+t57-t60-t37-t61-t62;
double      t64 = 32.0E0/27.0E0*x;
double      t65 = 32.0E0/81.0E0*t2;
double      t66 = 2.0E0/3.0E0*t54;
double      t67 = 2.0E0/3.0E0*t18;
double      t68 = 7.0E0/27.0E0*t7;
double      t69 = t59/3.0E0;
double      t70 = 32.0E0/81.0E0-t64-t65-t42-t43+t66+t56+t67+t68+t69+t37-t61-t39;
double      t71 = 32.0E0/81.0E0+t64-t65-t42-t43-t66+t56-t67+t68-t69+t37+t61-t39;
double      t72 = -16.0E0/243.0E0-t50+t51+t42-t52+t53+t55-t56+t57+t60-t37+t61-t62;
double      t73 = 11.0E0/324.0E0*x;
double      t74 = 11.0E0/972.0E0*t2;
double      t75 = 7.0E0/72.0E0*t28;
double      t76 = 7.0E0/108.0E0*t15;
double      t77 = 7.0E0/216.0E0*t4;
double      t78 = t54/36.0E0;
double      t79 = t33/36.0E0;
double      t80 = t18/36.0E0;
double      t81 = t7/324.0E0;
double      t82 = t28*t28;
double      t83 = t82/24.0E0;
double      t84 = t59/18.0E0;
double      t85 = t36/12.0E0;
double      t86 = t21/54.0E0;
double      t87 = t9/216.0E0;
double      t88 = 5.0E0/243.0E0-t73-t74-t75-t76-t77+t78+t79+t80+t81+t83+t84+t85+t86+t87;
double      t89 = 8.0E0/81.0E0*x;
double      t90 = 16.0E0/243.0E0*t2;
double      t91 = t28/3.0E0;
double      t92 = 4.0E0/27.0E0*t15;
double      t93 = t4/27.0E0;
double      t94 = t54/18.0E0;
double      t95 = t33/6.0E0;
double      t96 = 5.0E0/18.0E0*t18;
double      t97 = 7.0E0/162.0E0*t7;
double      t98 = t82/6.0E0;
double      t99 = t21/27.0E0;
double      t100 = t9/54.0E0;
double      t101 = -16.0E0/243.0E0+t89+t90+t91+t92+t93-t94-t95-t96-t97-t98-t60+t99+t100;
double      t109 = -16.0E0/243.0E0-t89+t90+t91-t92+t93+t94-t95+t96-t97-t98+t60-t99+t100;
double      t110 = 5.0E0/243.0E0+t73-t74-t75+t76-t77-t78+t79-t80+t81+t83-t84+t85-t86+t87;

      weight[0] = 5.0E0/243.0E0+11.0E0/486.0E0*t2-7.0E0/54.0E0*t4-2.0E0/81.0E0*t7+2.0E0/27.0E0*
t9;
      weight[1] = -16.0E0/243.0E0+t12-t13+t16+t17-t19+t20-t22-t23;
      weight[2] = -16.0E0/243.0E0-t12-t13-t16+t17+t19+t20+t22-t23;
      weight[3] = t40;
      weight[4] = 32.0E0/81.0E0+64.0E0/81.0E0*t2-t42-t43-5.0E0/9.0E0*t33-11.0E0/27.0E0*t7
-2.0E0/3.0E0*t36+2.0E0/9.0E0*t9;
      weight[5] = t49;
      weight[6] = t63;
      weight[7] = t70;
      weight[8] = t71;
      weight[9] = t72;
      weight[10] = t88;
      weight[11] = t101;
      weight[12] = 4.0E0/81.0E0-10.0E0/81.0E0*t2-17.0E0/36.0E0*t28+11.0E0/36.0E0*t4+t34-t35+
t82/4.0E0-t36/6.0E0+t9/36.0E0;
      weight[13] = t109;
      weight[14] = t110;

return;
}

////////////////////////////////////////////////////////////////////////////////
void get_feq(const double rho, const double u, const double v, const double e,
            double feq[30])
{
double C, ep;
C = sqrt(2.0E0*Gamma_1*e);
ep = (2.0E0-Gamma)*e;

double      t1 = v*v;
double      t3 = sqrt(3.0E0);
double      t4 = t3*v;
double      t5 = 22.0E0*t4;
double      t7 = C*C;
double      t8 = t1*theta*t7;
double      t9 = 216.0E0*t8;
double      t10 = theta*t7;
double      t12 = theta*theta;
double      t13 = t7*t7;
double      t14 = t12*t13;
double      t15 = 27.0E0*t14;
double      t16 = t1*t1;
double      t17 = 72.0E0*t16;
double      t18 = t1*v;
double      t19 = t3*t18;
double      t21 = t4*t10;
double      t22 = 36.0E0*t21;
double      t24 = rho*(-126.0E0*t1+20.0E0+t5+t9-63.0E0*t10+t15+t17-24.0E0*t19-t22);
double      t25 = 1/theta;
double      t28 = t25*(1.0E0-theta*ep);
double      t32 = rho*(theta-1.0E0);
double      t43 = u*t1;
double      t44 = 108.0E0*t43;
double      t46 = u*theta*t7;
double      t47 = 54.0E0*t46;
double      t48 = 40.0E0*t4;
double      t49 = 60.0E0*t19;
double      t50 = 24.0E0*u;
double      t51 = 180.0E0*t1;
double      t52 = 90.0E0*t10;
double      t53 = t3*u;
double      t54 = t53*v;
double      t55 = 36.0E0*t54;
double      t56 = t53*t18;
double      t57 = 72.0E0*t56;
double      t58 = 90.0E0*t21;
double      t61 = t53*v*theta*t7;
double      t62 = 108.0E0*t61;
double      t63 = t44+t47+t48-t49+t9-t50+t17+t15-t51-t52-t55+t57-t58+32.0E0+t62;
double      t64 = rho*t63;
double      t67 = 4.0E0/81.0E0*u;
double      t68 = 20.0E0/243.0E0*t4;
double      t69 = 2.0E0/27.0E0*t54;
double      t70 = 10.0E0/27.0E0*t1;
double      t71 = 2.0E0/9.0E0*t43;
double      t72 = 10.0E0/81.0E0*t19;
double      t73 = 4.0E0/27.0E0*t56;
double      t74 = 4.0E0/27.0E0*t16;
double      t81 = -t9+t51-t15+t49-t17-t48+t58+t57+t62+t47+t44-t50-t55+t52-32.0E0;
double      t82 = rho*t81;
double      t91 = 162.0E0*t8;
double      t92 = 60.0E0*u;
double      t93 = u*u;
double      t95 = t93*theta*t7;
double      t96 = 54.0E0*t95;
double      t97 = t93*t1;
double      t98 = 108.0E0*t97;
double      t99 = 36.0E0*t93;
double      t100 = 36.0E0*t16;
double      t101 = 90.0E0*t1;
double      t102 = 20.0E0*t4;
double      t103 = 30.0E0*t19;
double      t104 = 126.0E0*t54;
double      t106 = t3*t93*v;
double      t107 = 90.0E0*t106;
double      t108 = 27.0E0*t10;
double      t109 = t91+16.0E0-t92+t96+t98+t99+t100+t15-t101+t102-t103+t62+t57-t104+t107
-t108;

double      t110 = rho*t109;
double      t113 = 5.0E0/27.0E0*u;
double      t114 = 5.0E0/81.0E0*t4;
double      t115 = t93/9.0E0;
double      t116 = 7.0E0/18.0E0*t54;
double      t117 = 5.0E0/18.0E0*t1;
double      t118 = 5.0E0/18.0E0*t106;
double      t119 = 5.0E0/54.0E0*t19;
double      t120 = t97/3.0E0;
double      t121 = 2.0E0/9.0E0*t56;
double      t122 = t16/9.0E0;
double      t123 = 4.0E0/81.0E0-t113+t114+t115-t116-t117+t118-t119+t120+t121+t122;
double      t130 = 27.0E0*t95;
double      t131 = 72.0E0*t21;
double      t133 = 18.0E0*t16;
double      t134 = 18.0E0*t1;
double      t136 = 27.0E0*t8;
double      t138 = 64.0E0*t4;
double      t139 = 54.0E0*t97;
double      t140 = 45.0E0*t106+t130+t131+18.0E0*t93-t133+t134+18.0E0*t10-t136+33.0E0*t19-t138
+t139-32.0E0;

double      t141 = rho*t140;
double      t145 = 2.0E0/9.0E0*t93;
double      t146 = 2.0E0/9.0E0*t1;
double      t157 = 16.0E0+t107-t108+t100+t92-t62-t57+t104+t15+t99+t98+t96+t91+t102-t103
-t101;

double      t158 = rho*t157;
double      t161 = 4.0E0/81.0E0+t113+t114+t115+t116-t117+t118-t119+t120-t121+t122;
double      t167 = 72.0E0*t1;
double      t168 = 162.0E0*t61;
double      t169 = 6.0E0*t19;
double      t170 = 135.0E0*t8;
double      t171 = 162.0E0*t97;
double      t172 = 81.0E0*t95;
double      t173 = 108.0E0*t54;
double      t174 = 54.0E0*t56;
double      t175 = t93*u;
double      t177 = t3*t175*v;
double      t178 = 54.0E0*t177;
double      t179 = 18.0E0*t21;
double      t180 = 54.0E0*t106;
double      t181 = 162.0E0*t46;
double      t182 = 8.0E0*t4;
double      t183 = 72.0E0*u;
double      t184 = 108.0E0*t93;
double      t185 = 108.0E0*t175;
double      t186 = -t167-t52+t168-t169+t170+t171+t172-t173+t174+t178+t179+t180+t181-
t182+32.0E0-t183-t184+t185+t133+t15;

double      t187 = rho*t186;
double      t190 = 4.0E0/27.0E0*u;
double      t191 = 4.0E0/243.0E0*t4;
double      t192 = 2.0E0/9.0E0*t54;
double      t193 = 4.0E0/27.0E0*t1;
double      t194 = 2.0E0/9.0E0*t175;
double      t195 = t106/9.0E0;
double      t196 = t19/81.0E0;
double      t197 = t177/9.0E0;
double      t198 = t56/9.0E0;
double      t199 = t16/27.0E0;
double      t200 = -16.0E0/243.0E0+t190+t191+t145+t192+t193-t194-t195+t196-t197-t120-t198
-t199;

double      t206 = 216.0E0*t46;
double      t207 = 42.0E0*t19;
double      t208 = 36.0E0*t1;
double      t209 = 36.0E0*t10;
double      t210 = 54.0E0*t61;
double      t211 = 18.0E0*t56;
double      t212 = 18.0E0*t106;
double      t213 = 192.0E0*u;
double      t214 = t44+t206-t138+t207-t136+t139+t130-t99+t185-t133-t208-t209+t210-
t211+t178+t131+t212+64.0E0-t213;

double      t215 = rho*t214;
double      t218 = 32.0E0/27.0E0*u;
double      t219 = 32.0E0/81.0E0*t4;
double      t220 = 2.0E0/3.0E0*t175;
double      t221 = 2.0E0/3.0E0*t43;
double      t222 = 7.0E0/27.0E0*t19;
double      t223 = t177/3.0E0;
double      t224 = 32.0E0/81.0E0-t218-t219-t145-t146+t220+t195+t221+t222+t223+t120-t198-
t122;

double      t230 = -64.0E0+t208+t209-t213+t99+t185+t133-t207+t44-t130+t136+t138+t178-
t211+t210-t131+t206-t212-t139;

double      t231 = rho*t230;
double      t234 = 32.0E0/81.0E0+t218-t219-t145-t146-t220+t195-t221+t222-t223+t120+t198-
t122;

double      t240 = t184-t15-t133+t167-32.0E0+t52-t171-t172-t170-t179-t180+t185+t182+
t169+t178+t168-t173+t174+t181-t183;

double      t241 = rho*t240;
double      t244 = -16.0E0/243.0E0-t190+t191+t145-t192+t193+t194-t195+t196+t197-t120+t198
-t199;

double      t250 = 36.0E0*t56;
double      t251 = 108.0E0*t177;
double      t252 = 216.0E0*t61;
double      t253 = 54.0E0*t43;
double      t254 = 108.0E0*t46;
double      t255 = 108.0E0*t8;
double      t256 = -t104+t250+t251+t22+t180+t252+t253+t254-t5+t169+t255;
double      t257 = 324.0E0*t95;
double      t258 = 66.0E0*u;
double      t259 = 189.0E0*t93;
double      t260 = 54.0E0*t175;
double      t261 = 9.0E0*t16;
double      t262 = 54.0E0*t14;
double      t263 = 63.0E0*t1;
double      t264 = t93*t93;
double      t265 = 81.0E0*t264;
double      t266 = 126.0E0*t10;
double      t267 = 40.0E0+t171+t257-t258-t259+t260+t261+t262-t263+t265-t266;
double      t269 = rho*(t256+t267);
double      t272 = 11.0E0/324.0E0*u;
double      t273 = 11.0E0/972.0E0*t4;
double      t274 = 7.0E0/72.0E0*t93;
double      t275 = 7.0E0/108.0E0*t54;
double      t276 = 7.0E0/216.0E0*t1;
double      t277 = t175/36.0E0;
double      t278 = t106/36.0E0;
double      t279 = t43/36.0E0;
double      t280 = t19/324.0E0;
double      t281 = t264/24.0E0;
double      t282 = t177/18.0E0;
double      t283 = t97/12.0E0;
double      t284 = t56/54.0E0;
double      t285 = t16/216.0E0;
double      t286 = 5.0E0/243.0E0-t272-t273-t274-t275-t276+t277+t278+t279+t280+t281+t282+
t283+t284+t285;

double      t292 = 48.0E0*u;
double      t293 = 162.0E0*t93;
double      t294 = 27.0E0*t175;
double      t296 = 243.0E0*t95;
double      t297 = 135.0E0*t43;
double      t298 = 32.0E0*t4;
double      t299 = 21.0E0*t19;
double      t300 = 72.0E0*t54;
double      t301 = 81.0E0*t106;
double      t302 = -t211+t178+t131+t296+t297+t254-t298+t299-t136-t300+t301;
double      t304 = rho*(32.0E0-t261-t134+t265-t52+t15-t292-t293+t294+t210+t302);
double      t307 = 8.0E0/81.0E0*u;
double      t308 = 16.0E0/243.0E0*t4;
double      t309 = t93/3.0E0;
double      t310 = 4.0E0/27.0E0*t54;
double      t311 = t1/27.0E0;
double      t312 = t175/18.0E0;
double      t313 = t106/6.0E0;
double      t314 = 5.0E0/18.0E0*t43;
double      t315 = 7.0E0/162.0E0*t19;
double      t316 = t264/6.0E0;
double      t317 = t56/27.0E0;
double      t318 = t16/54.0E0;
double      t319 = -16.0E0/243.0E0+t307+t308+t309+t310+t311-t312-t313-t314-t315-t316-t197
+t317+t318;

double      t328 = 16.0E0+t107-t103-t139+216.0E0*t95-t48+99.0E0*t1+t15+t261-153.0E0*t93-t108+
t265;

double      t329 = rho*t328;
double      t345 = 32.0E0-t261+t300+t211+t292+t15-t134-t52+t299-t136+t296;
double      t347 = rho*(t265+t131+t301-t298-t178-t293-t210-t294-t297-t254+t345);
double      t350 = -16.0E0/243.0E0-t307+t308+t309-t310+t311+t312-t313+t314-t315-t316+t197
-t317+t318;

double      t356 = -t252-t253-t254-t5+t169+t255+t171+t257+t22+t180+t104;
double      t357 = 40.0E0-t250-t251+t258-t259-t260+t261+t262-t263+t265-t266;
double      t359 = rho*(t356+t357);
double      t362 = 5.0E0/243.0E0+t272-t273-t274+t275-t276-t277+t278-t279+t280+t281-t282+
t283-t284+t285;


      feq[0] = t24*t28/972.0E0+t32*t25*(5.0E0/243.0E0+11.0E0/486.0E0*t4-7.0E0/54.0E0*t1-2.0E0/
81.0E0*t19+2.0E0/27.0E0*t16);
      feq[1] = t24*ep/972.0E0;
      feq[2] = -t64*t28/486.0E0+t32*t25*(-16.0E0/243.0E0+t67-t68+t69+t70-t71+t72-t73-
t74);
      feq[3] = -t64*ep/486.0E0;
      feq[4] = t82*t28/486.0E0+t32*t25*(-16.0E0/243.0E0-t67-t68-t69+t70+t71+t72+t73-
t74);
      feq[5] = t82*ep/486.0E0;
      feq[6] = t110*t28/324.0E0+t32*t25*t123;
      feq[7] = t110*ep/324.0E0;
      feq[8] = -t141*t28/81.0E0+t32*t25*(32.0E0/81.0E0+64.0E0/81.0E0*t4-t145-t146-5.0E0/9.0E0
*t106-11.0E0/27.0E0*t19-2.0E0/3.0E0*t97+2.0E0/9.0E0*t16);
      feq[9] = -t141*ep/81.0E0;
      feq[10] = t158*t28/324.0E0+t32*t25*t161;
      feq[11] = t158*ep/324.0E0;
      feq[12] = -t187*t28/486.0E0+t32*t25*t200;
      feq[13] = -t187*ep/486.0E0;
      feq[14] = t215*t28/162.0E0+t32*t25*t224;
      feq[15] = t215*ep/162.0E0;
      feq[16] = -t231*t28/162.0E0+t32*t25*t234;
      feq[17] = -t231*ep/162.0E0;
      feq[18] = t241*t28/486.0E0+t32*t25*t244;
      feq[19] = t241*ep/486.0E0;
      feq[20] = t269*t28/1944.0E0+t32*t25*t286;
      feq[21] = t269*ep/1944.0E0;
      feq[22] = -t304*t28/486.0E0+t32*t25*t319;
      feq[23] = -t304*ep/486.0E0;
      feq[24] = t329*t28/324.0E0+t32*t25*(4.0E0/81.0E0-10.0E0/81.0E0*t4-17.0E0/36.0E0*t93+
11.0E0/36.0E0*t1+t118-t119+t264/4.0E0-t97/6.0E0+t16/36.0E0);
      feq[25] = t329*ep/324.0E0;
      feq[26] = -t347*t28/486.0E0+t32*t25*t350;
      feq[27] = -t347*ep/486.0E0;
      feq[28] = t359*t28/1944.0E0+t32*t25*t362;
      feq[29] = t359*ep/1944.0E0;


return;
}


void get_vstress(double ruve[4], double f[30], double stress[3])
{
    
    double tfeq[30];
    get_feq(ruve[0],ruve[1],ruve[2],ruve[3],tfeq);
    
    stress[0]=0.0E0;
    stress[1]=0.0E0;
    stress[2]=0.0E0;
    
    for(int dir=0;dir<30;dir++) 
    {
        tfeq[dir] = f[dir] - tfeq[dir];
        stress[0] += tfeq[dir]*dc[dir][0]*dc[dir][0];
        stress[1] += tfeq[dir]*dc[dir][0]*dc[dir][1];
        stress[2] += tfeq[dir]*dc[dir][1]*dc[dir][1];
    }
    
    return;
    
}
