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

#include "d2q17l2.h"
double e_level[2], Gamma, lbm_omega, Gamma_1,
        Pr,theta, Mum, Kam, Rgas, Cv, Cp;



void macro(const double c_f[34], 
           double& rr, double& ru, double& rv, double& re)
{
    rr = 0.0E0;
    ru = 0.0E0;
    rv = 0.0E0;
    re = 0.0E0;
    
    double mi;
    for(int idir=0;idir<17;idir++)
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
             const double cf[34],              // the original f
             double rf[34],                  // the f for output
             const int mode
           )
{
    for(int idir=0;idir<34;idir++) 
        rf[idir] = 0.0E0;
    
    
    double rx,ry,weight[15];
    for(int dir=0;dir<17;dir++)
    {
        // compute reflected location (rx,ry) of node[dir]
        reflect(cnx,cny, dc[dir+dir][0],dc[dir+dir][1],  rx,ry, mode);
        get_weight(rx,ry,weight);       
        // assign f(rx,ry) onto new nodes
        for(int idir=0;idir<17;idir++)
        {
            rf[idir+idir] += weight[idir]*cf[dir+dir]; 
            rf[idir+idir+1] += weight[idir]*cf[dir+dir+1]; 
        }
    }

    return;
}



////////////////////////////////////////////////////////////////////////////////
void collision(double* MV,           
               double c_bf[34],
               double c_df[34],
               const double dt,
               const double c_beta)
{
    double c_fq[34], tdf;
    
    get_feq(MV[0],MV[1],MV[2],MV[3],c_fq);

#ifdef LAMINAR
    // compute Tau according to pressure and molecule viscosity
    lbm_omega = MV[0]*MV[3]*Gamma_1 / Mum;
    MV[4]=lbm_omega*dt; 
#endif

#ifdef INVISCID
    
#endif

    for(int dir=0;dir<34;dir++)
    {
        tdf = ( c_fq[dir] - c_bf[dir] )*lbm_omega*dt*c_beta;
        c_bf[dir] += tdf;
        c_df[dir] += tdf;
    }
}


void get_weight(const double x, const double y, 
        double weight[17])
{
double t1 = x*x;
double t3 = y*y;
double t5 = t1*t1;
double t7 = t1*t3;
double t9 = t3*t3;
double t12 = 27.0/40.0*x;
double t13 = 27.0/40.0*y;
double t14 = 9.0/16.0*t1;
double t15 = x*y;
double t16 = 15.0/16.0*t15;
double t17 = 9.0/16.0*t3;
double t18 = t1*x;
double t19 = 27.0/64.0*t18;
double t21 = 9.0/32.0*t1*y;
double t23 = 9.0/32.0*x*t3;
double t24 = t3*y;
double t25 = 27.0/64.0*t24;
double t26 = 9.0/64.0*t5;
double t27 = t18*y;
double t28 = 27.0/64.0*t27;
double t29 = 63.0/64.0*t7;
double t30 = x*t24;
double t31 = 27.0/64.0*t30;
double t32 = 9.0/64.0*t9;
double t33 = t5*x;
double t34 = 81.0/1280.0*t33;
double t35 = t9*y;
double t36 = 81.0/1280.0*t35;
double t37 = t12+t13+t14+t16+t17-t19-t21-t23-t25-t26-t28-t29-t31-t32+t34+t36;
double t38 = -t12+t13+t14-t16+t17+t19-t21+t23-t25-t26+t28-t29+t31-t32-t34+t36;
double t39 = -t12-t13+t14+t16+t17+t19+t21+t23+t25-t26-t28-t29-t31-t32-t34-t36;
double t40 = t12-t13+t14-t16+t17-t19+t21-t23+t25-t26+t28-t29+t31-t32+t34-t36;
double t41 = 9.0/80.0*x;
double t42 = 9.0/80.0*y;
double t43 = 9.0/32.0*t1;
double t44 = 3.0/32.0*t15;
double t45 = 9.0/64.0*t3;
double t46 = 9.0/256.0*t18;
double t47 = 9.0/32.0*t24;
double t48 = 9.0/128.0*t5;
double t49 = 9.0/32.0*t7;
double t50 = 27.0/128.0*t30;
double t51 = 9.0/256.0*t9;
double t52 = 81.0/5120.0*t33;
double t53 = -t41-t42-t43-t44+t45-t46+t23+t47+t48+t49+t50-t51+t52-t36;
double t54 = 9.0/64.0*t1;
double t55 = 9.0/32.0*t3;
double t56 = 9.0/32.0*t18;
double t57 = 9.0/256.0*t24;
double t58 = 9.0/256.0*t5;
double t59 = 27.0/128.0*t27;
double t60 = 9.0/128.0*t9;
double t61 = 81.0/5120.0*t35;
double t62 = -t41-t42+t54-t44-t55+t56+t21-t57-t58+t59+t49+t60-t34+t61;
double t63 = t41-t42-t43+t44+t45+t46-t23+t47+t48+t49-t50-t51-t52-t36;
double t64 = t41-t42+t54+t44-t55-t56+t21-t57-t58-t59+t49+t60+t34+t61;
double t65 = t41+t42-t43-t44+t45+t46-t23-t47+t48+t49+t50-t51-t52+t36;
double t66 = t41+t42+t54-t44-t55-t56-t21+t57-t58+t59+t49+t60+t34-t61;
double t67 = -t41+t42-t43+t44+t45-t46+t23-t47+t48+t49-t50-t51+t52+t36;
double t68 = -t41+t42+t54+t44-t55+t56-t21+t57-t58-t59+t49+t60-t34-t61;
double t69 = x/40.0;
double t70 = t1/16.0;
double t71 = 3.0/32.0*t3;
double t72 = 9.0/128.0*t18;
double t73 = 3.0/64.0*t5;
double t74 = 3.0/32.0*t7;
double t75 = 3.0/128.0*t9;
double t76 = 81.0/2560.0*t33;
double t78 = y/40.0;
double t79 = 3.0/32.0*t1;
double t80 = t3/16.0;
double t81 = 9.0/128.0*t24;
double t82 = 3.0/128.0*t5;
double t83 = 3.0/64.0*t9;
double t84 = 81.0/2560.0*t35;


      weight[0] = 1.0-7.0/4.0*t1-7.0/4.0*t3+3.0/8.0*t5+33.0/16.0*t7+3.0/8.0*
9;
      weight[1] = t37;
      weight[2] = t38;
      weight[3] = t39;
      weight[4] = t40;
      weight[5] = t53;
      weight[6] = t62;
      weight[7] = t63;
      weight[8] = t64;
      weight[9] = t65;
      weight[10] = t66;
      weight[11] = t67;
      weight[12] = t68;
      weight[13] = t69-t70+t71-t72+t73-t74-t75+t76;
      weight[14] = t78+t79-t80-t81-t82-t74+t83+t84;
      weight[15] = -t69-t70+t71+t72+t73-t74-t75-t76;
      weight[16] = -t78+t79-t80+t81-t82-t74+t83-t84;


return;
}

////////////////////////////////////////////////////////////////////////////////
void get_feq(const double rho, const double u, const double v, const double e,
            double feq[34])
{
double C, ep;
C = sqrt(2.0E0*Gamma_1*e);
ep = (2.0E0-Gamma)*e;

double      t1 = theta*theta;
double      t2 = C*C;
double      t3 = t2*t2;
double      t4 = t1*t3;
double      t6 = u*u;
double      t7 = v*v;
double      t8 = t6*t7;
double      t11 = t7*t7;
double      t14 = theta*t2;
double      t15 = t14*t7;
double      t17 = t6*t6;
double      t20 = t14*t6;
double      t23 = rho*(69.0*t4+128.0+264.0*t8-224.0*t7+48.0*t11-224.0*t6+276.0
*t15+48.0*t17-224.0*t14+276.0*t20);
double      t24 = 1/theta;
double      t27 = t24*(1.0-theta*ep);
double      t31 = rho*(theta-1.0);
double      t43 = 780.0*t4;
double      t44 = u*v;
double      t46 = 4320.0*t14*t44;
double      t47 = 1920.0*t14;
double      t48 = u*t7;
double      t49 = 960.0*t48;
double      t50 = 3200.0*t44;
double      t51 = t6*u;
double      t52 = t51*v;
double      t53 = 1440.0*t52;
double      t54 = t11*v;
double      t55 = 216.0*t54;
double      t56 = t17*u;
double      t57 = 216.0*t56;
double      t58 = t7*v;
double      t59 = u*t58;
double      t60 = 1440.0*t59;
double      t61 = 1440.0*t51;
double      t62 = 1440.0*t58;
double      t63 = t4*u;
double      t64 = 405.0*t63;
double      t65 = t14*u;
double      t66 = 2640.0*t65;
double      t67 = -t43-t46+t47-t49+t50-t53+t55+t57-t60-t61-t62+t64-t66;
double      t69 = v*theta*t2;
double      t70 = 2640.0*t69;
double      t71 = t14*t58;
double      t72 = 1080.0*t71;
double      t73 = t4*v;
double      t74 = 405.0*t73;
double      t75 = t14*t51;
double      t76 = 1080.0*t75;
double      t77 = 1920.0*t7;
double      t78 = 480.0*t11;
double      t79 = 480.0*t17;
double      t80 = 2304.0*u;
double      t81 = 2304.0*v;
double      t82 = 3360.0*t8;
double      t83 = 1920.0*t6;
double      t84 = 3120.0*t15;
double      t85 = 3120.0*t20;
double      t86 = t6*v;
double      t87 = 960.0*t86;
double      t88 = -t70+t72+t74+t76+t77-t78-t79+t80+t81-t82+t83-t84-t85-t87;
double      t90 = rho*(t67+t88);
double      t93 = 27.0/40.0*u;
double      t94 = 27.0/40.0*v;
double      t95 = 9.0/16.0*t6;
double      t96 = 15.0/16.0*t44;
double      t97 = 9.0/16.0*t7;
double      t98 = 27.0/64.0*t51;
double      t99 = 9.0/32.0*t86;
double      t100 = 9.0/32.0*t48;
double      t101 = 27.0/64.0*t58;
double      t102 = 9.0/64.0*t17;
double      t103 = 27.0/64.0*t52;
double      t104 = 63.0/64.0*t8;
double      t105 = 27.0/64.0*t59;
double      t106 = 9.0/64.0*t11;
double      t107 = 81.0/1280.0*t56;
double      t108 = 81.0/1280.0*t54;
double      t109 = t93+t94+t95+t96+t97-t98-t99-t100-t101-t102-t103-t104-t105
-t106+t107+t108;
double      t115 = -t43+t46+t47+t49-t50+t53+t55-t57+t60+t61-t62-t64+t66;
double      t116 = -t70+t72+t74-t76+t77-t78-t79-t80+t81-t82+t83-t84-t85-t87;
double      t118 = rho*(t115+t116);
double      t121 = -t93+t94+t95-t96+t97+t98-t99+t100-t101-t102+t103-t104+t105
-t106-t107+t108;
double      t127 = t43+t46-t47-t49-t50+t53+t55+t57+t60-t61-t62+t64-t66;
double      t128 = -t70+t72+t74+t76-t77+t78+t79+t80+t81+t82-t83+t84+t85-t87;
double      t130 = rho*(t127+t128);
double      t133 = -t93-t94+t95+t96+t97+t98+t99+t100+t101-t102-t103-t104-t105
-t106-t107-t108;
double      t139 = t43-t46-t47+t49+t50-t53+t55-t57-t60+t61-t62-t64+t66;
double      t140 = -t70+t72+t74-t76-t77+t78+t79-t80+t81+t82-t83+t84+t85-t87;
double      t142 = rho*(t139+t140);
double      t145 = t93-t94+t95-t96+t97-t98+t99-t100+t101-t102+t103-t104+t105
-t106+t107-t108;
double      t151 = 660.0*t4;
double      t152 = 960.0*t14;
double      t153 = 3840.0*t48;
double      t154 = 1280.0*t44;
double      t155 = 864.0*t54;
double      t156 = 2880.0*t59;
double      t157 = 480.0*t51;
double      t158 = 3840.0*t58;
double      t159 = 1200.0*t65;
double      t160 = -t151-t46+t152-t153+t154+t155-t57-t156+t157-t158-t64-t159;
double      t161 = 5760.0*t69;
double      t162 = 4320.0*t71;
double      t163 = 1620.0*t73;
double      t164 = 960.0*t17;
double      t165 = 1536.0*u;
double      t166 = 1536.0*v;
double      t167 = 3840.0*t8;
double      t168 = 3840.0*t6;
double      t169 = 480.0*t15;
double      t170 = 4800.0*t20;
double      t171 = -t161+t162+t163-t76-t77+t78-t164+t165+t166-t167+t168-t169
-t170;
double      t173 = rho*(t160+t171);
double      t176 = 9.0/80.0*u;
double      t177 = 9.0/80.0*v;
double      t178 = 9.0/32.0*t6;
double      t179 = 3.0/32.0*t44;
double      t180 = 9.0/64.0*t7;
double      t181 = 9.0/256.0*t51;
double      t182 = 9.0/32.0*t58;
double      t183 = 9.0/128.0*t17;
double      t184 = 9.0/32.0*t8;
double      t185 = 27.0/128.0*t59;
double      t186 = 9.0/256.0*t11;
double      t187 = 81.0/5120.0*t56;
double      t188 = -t176-t177-t178-t179+t180-t181+t100+t182+t183+t184+t185-t186
+t187-t108;
double      t194 = 2880.0*t52;
double      t195 = 864.0*t56;
double      t196 = 3840.0*t51;
double      t197 = 480.0*t58;
double      t198 = 1620.0*t63;
double      t199 = 5760.0*t65;
double      t200 = 1200.0*t69;
double      t201 = t151+t46-t152-t154+t194+t55-t195+t196-t197-t198+t199+t200;
double      t202 = 4320.0*t75;
double      t203 = 3840.0*t7;
double      t204 = 960.0*t11;
double      t205 = 4800.0*t15;
double      t206 = 480.0*t20;
double      t207 = 3840.0*t86;
double      t208 = t72+t74-t202-t203+t204-t79-t165-t166+t167+t83+t205+t206+t207;
double      t210 = rho*(t201+t208);
double      t213 = 9.0/64.0*t6;
double      t214 = 9.0/32.0*t7;
double      t215 = 9.0/32.0*t51;
double      t216 = 9.0/256.0*t58;
double      t217 = 9.0/256.0*t17;
double      t218 = 27.0/128.0*t52;
double      t219 = 9.0/128.0*t11;
double      t220 = 81.0/5120.0*t54;
double      t221 = -t176-t177+t213-t179-t214+t215+t99-t216-t217+t218+t184+t219-
t107+t220;
double      t227 = -t151+t46+t152+t153-t154+t155+t57+t156-t157-t158+t64+t159;
double      t228 = -t161+t162+t163+t76-t77+t78-t164-t165+t166-t167+t168-t169
-t170;
double      t230 = rho*(t227+t228);
double      t233 = t176-t177-t178+t179+t180+t181-t100+t182+t183+t184-t185-t186
-t187-t108;
double      t239 = t151-t46-t152+t154-t194+t55+t195-t196-t197+t198-t199+t200;
double      t240 = t72+t74+t202-t203+t204-t79+t165-t166+t167+t83+t205+t206+t207;
double      t242 = rho*(t239+t240);
double      t245 = t176-t177+t213+t179-t214-t215+t99-t216-t217-t218+t184+t219
+t107+t220;
double      t251 = t151+t46-t152-t153-t154+t155-t57+t156+t157-t158-t64-t159;
double      t252 = -t161+t162+t163-t76+t77-t78+t164+t165+t166+t167-t168+t169
+t170;
double      t254 = rho*(t251+t252);
double      t257 = t176+t177-t178-t179+t180+t181-t100-t182+t183+t184+t185-t186
-t187+t108;
double      t263 = -t151-t46+t152+t154-t194+t55-t195+t196-t197-t198+t199+t200;
double      t264 = t72+t74-t202+t203-t204+t79-t165-t166-t167-t83-t205-t206+t207;
double      t266 = rho*(t263+t264);
double      t269 = t176+t177+t213-t179-t214-t215-t99+t216-t217+t218+t184+t219
+t107-t220;
double      t275 = t151-t46-t152+t153+t154+t155+t57-t156-t157-t158+t64+t159;
double      t276 = -t161+t162+t163+t76+t77-t78+t164-t165+t166+t167-t168+t169
+t170;
double      t278 = rho*(t275+t276);
double      t281 = -t176+t177-t178+t179+t180-t181+t100-t182+t183+t184-t185-t186
+t187+t108;
double      t287 = -t151+t46+t152-t154+t194+t55+t195-t196-t197+t198-t199+t200;
double      t288 = t72+t74+t202+t203-t204+t79+t165-t166-t167-t83-t205-t206+t207;
double      t290 = rho*(t287+t288);
double      t293 = -t176+t177+t213+t179-t214+t215-t99+t216-t217-t218+t184+t219
-t107-t220;
double      t299 = 1920.0*t8;
double      t300 = 2400.0*t15;
double      t301 = 1920.0*t20;
double      t302 = 1280.0*t6;
double      t303 = 320.0*t14;
double      t304 = 60.0*t4;
double      t305 = 648.0*t56;
double      t306 = 512.0*u;
double      t307 = 1215.0*t63;
double      t308 = 3240.0*t75;
double      t309 = 2160.0*t65;
double      t310 = t299+t300-t301+t78+t302-t164-t303+t304-t77-t305-t306-t307
-t308+t309+t61;
double      t311 = rho*t310;
double      t314 = u/40.0;
double      t315 = t6/16.0;
double      t316 = 3.0/32.0*t7;
double      t317 = 9.0/128.0*t51;
double      t318 = 3.0/64.0*t17;
double      t319 = 3.0/32.0*t8;
double      t320 = 3.0/128.0*t11;
double      t321 = 81.0/2560.0*t56;
double      t328 = 1280.0*t7;
double      t329 = 512.0*v;
double      t330 = 648.0*t54;
double      t331 = 1215.0*t73;
double      t332 = 2160.0*t69;
double      t333 = 1920.0*t15;
double      t334 = 2400.0*t20;
double      t335 = 3240.0*t71;
double      t336 = t204+t83-t79+t303-t304-t328+t329-t62+t330+t331-t332-t299+t333
-t334+t335;
double      t337 = rho*t336;
double      t340 = v/40.0;
double      t341 = 3.0/32.0*t6;
double      t342 = t7/16.0;
double      t343 = 9.0/128.0*t58;
double      t344 = 3.0/128.0*t17;
double      t345 = 3.0/64.0*t11;
double      t346 = 81.0/2560.0*t54;
double      t353 = t299+t300-t301+t78+t302-t164-t303+t304-t77+t305+t306+t307
+t308-t309-t61;
double      t354 = rho*t353;
double      t363 = -t204-t83+t79-t303+t304+t328+t329-t62+t330+t299-t333+t334
+t335+t331-t332;
double      t364 = rho*t363;
      
      
      feq[0] = t23*t27/128.0+t31*t24*(1.0-7.0/4.0*t6-7.0/4.0*t7+3.0/8.0*t17
      +33.0/16.0*t8+3.0/8.0*t11);
      feq[1] = t23*ep/128.0;
      feq[2] = 3.0/10240.0*t90*t27+t31*t24*t109;
      feq[3] = 3.0/10240.0*t90*ep;
      feq[4] = 3.0/10240.0*t118*t27+t31*t24*t121;
      feq[5] = 3.0/10240.0*t118*ep;
      feq[6] = -3.0/10240.0*t130*t27+t31*t24*t133;
      feq[7] = -3.0/10240.0*t130*ep;
      feq[8] = -3.0/10240.0*t142*t27+t31*t24*t145;
      feq[9] = -3.0/10240.0*t142*ep;
      feq[10] = -3.0/40960.0*t173*t27+t31*t24*t188;
      feq[11] = -3.0/40960.0*t173*ep;
      feq[12] = 3.0/40960.0*t210*t27+t31*t24*t221;
      feq[13] = 3.0/40960.0*t210*ep;
      feq[14] = -3.0/40960.0*t230*t27+t31*t24*t233;
      feq[15] = -3.0/40960.0*t230*ep;
      feq[16] = 3.0/40960.0*t242*t27+t31*t24*t245;
      feq[17] = 3.0/40960.0*t242*ep;
      feq[18] = 3.0/40960.0*t254*t27+t31*t24*t257;
      feq[19] = 3.0/40960.0*t254*ep;
      feq[20] = -3.0/40960.0*t266*t27+t31*t24*t269;
      feq[21] = -3.0/40960.0*t266*ep;
      feq[22] = 3.0/40960.0*t278*t27+t31*t24*t281;
      feq[23] = 3.0/40960.0*t278*ep;
      feq[24] = -3.0/40960.0*t290*t27+t31*t24*t293;
      feq[25] = -3.0/40960.0*t290*ep;
      feq[26] = -t311*t27/20480.0+t31*t24*(t314-t315+t316-t317+t318-t319-t320
      +t321);
      feq[27] = -t311*ep/20480.0;
      feq[28] = t337*t27/20480.0+t31*t24*(t340+t341-t342-t343-t344-t319+t345
      +t346);
      feq[29] = t337*ep/20480.0;
      feq[30] = -t354*t27/20480.0+t31*t24*(-t314-t315+t316+t317+t318-t319-t320
      -t321);
      feq[31] = -t354*ep/20480.0;
      feq[32] = -t364*t27/20480.0+t31*t24*(-t340+t341-t342+t343-t344-t319+t345
      -t346);
      feq[33] = -t364*ep/20480.0;


return;
}


void get_vstress(double ruve[4], double f[34], double stress[3])
{
    
    double tfeq[34];
    get_feq(ruve[0],ruve[1],ruve[2],ruve[3],tfeq);
    
    stress[0]=0.0E0;
    stress[1]=0.0E0;
    stress[2]=0.0E0;
    
    for(int dir=0;dir<34;dir++) 
    {
        tfeq[dir] = f[dir] - tfeq[dir];
        stress[0] += tfeq[dir]*dc[dir][0]*dc[dir][0];
        stress[1] += tfeq[dir]*dc[dir][0]*dc[dir][1];
        stress[2] += tfeq[dir]*dc[dir][1]*dc[dir][1];
    }
    
    return;
    
}
