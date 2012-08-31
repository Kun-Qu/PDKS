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

#include "explicit.h" 
#include "bc.h" 
 
double ref_e,ref_U,ref_rho,ref_p, ref_L; 
double Re_L,Re_rho,Re_U;    
int flag[IDIM+3][JDIM+3]; 
double DQP[MAXDIM+3],DQM[MAXDIM+3],QP[MAXDIM+3],QM[MAXDIM+3];//,FLUX[MAXDIM+3]; 
 
 
 
 
/* ////////////////////////////////////////////////// 
Euler forward integration 
*/ 
void step( const double dt, const double c_alpha, const double c_beta) 
{ 
    //std::cout<<myrank<<" SYNCH = "<<get_SYNCH()<<std::endl; 
    if(get_SYNCH()==0) return; 
	    
    MPI_Barrier(MPI_COMM_WORLD); 
    //std::cout<<myrank<<"  MPI_Barrier(MPI_COMM_WORLD);"<<std::endl; 
	  	 
	  	 
    convect(dt, c_alpha); 
    //std::cout<<myrank<<"  convect(dt, c_alpha);"<<std::endl; 
    	 
    update(); 
    //std::cout<<myrank<<"  update();"<<std::endl; 
    	 
    source(dt,  c_beta); 
    //std::cout<<myrank<<"  source(dt,  c_beta);"<<std::endl; 
     
    clr_SYNCH(); 
    clr_GATHR(); 
    clr_ALL_UPDATED(); 
    //std::cout<<myrank<<"  Clear flags;"<<std::endl; 
    	 
    return; 
} 
 
 
 
 
//////////////////////////////////////////////////////////////////////////////// 
/* 
Just change Q 
*/ 
void source(const double dt, const double c_beta) 
{ 
	int iis,iit,jjs,jjt; 
	my_range(myrank,iis,iit,jjs,jjt); 
     
    for(int i=iis;i<=iit;i++) 
    for(int j=jjs;j<=jjt;j++) 
    { 
        collision(M[i][j],   Q[i][j], DQ[i][j],   dt, c_beta); 
    } 
     
    return; 
} 
 
 
/* ///////////////////////////////////////////////////////////////// 
 
 
*/ 
int get_global_nfo(std::string& dataname,  
	                 std::ofstream& logfile, 
	                 double& t,
	                 const double dt, 
	                 int& iter,  
	                 const int save_int,  
	                 const int check_int,  
	                 const int displ_int,
                         const char anyway) 
{ 
const char DISPLAY_ANYWAY= 0x1;
const char SAVE_ANYWAY   = 0x2;
const char CHECK_ANYWAY  = 0x4;

   if(iter%displ_int==0 || (anyway & DISPLAY_ANYWAY) )
   { 
       double res; 
       res = get_nsq_residual_P(dt);  
       //compute_force(dataname,2,t,iter, 0, &Cx, &Cy); 
       if(myrank==0)
       {
           double rr,rU,re,rL,rt;
           get_ref(rr,rU,re,rL);
           rt = rL/rU;
           logfile<<iter<<" "<<res<<std::endl; 
            std::cout<<"\nIterating "<<iter<<"  Time "<<t*rt<<" residual = "<<res<<std::endl; 
       }
   } 
    
    
   if(iter%save_int==0 || (anyway & SAVE_ANYWAY) )
   { 
   	gather_array();
   	//std::cout<<myrank<<", I wanto to save"<<std::endl;
   	if(myrank==0)
   	{
        for(int i=2;i<=IDIM;i++)
        for(int j=2;j<=JDIM;j++)
           macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
   	    bc(t);
   	    set_ALL_UPDATED();
    }
    
   	
   	write_plt(dataname,t,iter,true);
        if(myrank==0)
        {
   	      double Cx, Cy;
              int writeforce=1;
   	      compute_force(dataname,4,t,iter, writeforce, &Cx, &Cy);
   	      std::cout<<"\nIterating "<<iter<<"  Cx = "<<Cx<<" Cy = "<<Cy<<std::endl;
        }
   } 
 
   if(iter%check_int==0 || (anyway & CHECK_ANYWAY) )
   { 
        gather_array();
        //std::cout<<myrank<<"I wanto to write checkpoint file"<<std::endl;
        if(myrank==0) wr_check_2_0(dataname,t,iter,0); 
   } 
   	    
   return 1; 
	 
} 
 
/* ////////////////////////////////////
Compute residual after gathering
*/ 
 
double get_inf_residual(const double dt) 
{ 
    if(myrank!=0)
    {
      std::cerr<<"Only available for Rank 0"<<std::endl;
      return 1.0E0;
     }
    
    double tmax=0.0E0, tabs; 
     
    for(int i=2;i<=IDIM;i++) 
    for(int j=2;j<=JDIM;j++) 
    { 
        for(int dir=0;dir<DIR;dir++) 
        { 
            tabs = fabs(DQ[i][j][dir]); 
            tmax = tmax>tabs? tmax: tabs; 
        } 
    } 
 
    return tmax/dt; 
     
} 

double get_nsq_residual(const double dt) 
{ 
    if(myrank!=0)
    {
      std::cerr<<"Only available for Rank 0"<<std::endl;
      return 1.0E0;
     }

    
    double sum_sq=0.0E0; 
    for(int i=2;i<=IDIM;i++) 
    for(int j=2;j<=JDIM;j++) 
    { 
 
        for(int dir=0;dir<DIR;dir++) 
            sum_sq += DQ[i][j][dir] * DQ[i][j][dir]; 
    } 
 
    sum_sq = sqrt(sum_sq/(IDIM*JDIM))/dt; 
    return sum_sq; 
  
}

/* ////////////////////////////////////
Compute residual distributely
*/ 
double get_inf_residual_P(const double dt) 
{ 
    int iis,iit,jjs,jjt; 
    my_range(myrank,iis,iit,jjs,jjt); 
	
    double tmax=0.0E0, tabs; 
     
    for(int i=iis;i<=iit;i++) 
    for(int j=jjs;j<=jjt;j++) 
    { 
        for(int dir=0;dir<DIR;dir++) 
        { 
            tabs = fabs(DQ[i][j][dir]); 
            tmax = tmax>tabs? tmax: tabs; 
        } 
    } 
    //std::cout<<myrank<<" tmax = "<<tmax<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank==0)
       for(int rank=1;rank<cpusize;rank++)
       {
           double tt;
           MPI_Status status;
           MPI_Recv( &tt, 1, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &status );
           if(tt>tmax) tmax = tt;
       }
    else
         MPI_Send( &tmax, 1, MPI_DOUBLE, 0,    myrank, MPI_COMM_WORLD);
     
    if(myrank==0)
      return tmax/dt;
    else
      return 1.0E0; 
     
} 

 
double get_nsq_residual_P(const double dt) 
{ 
    int iis,iit,jjs,jjt; 
    my_range(myrank,iis,iit,jjs,jjt);
    
    double sum_sq=0.0E0; 
    for(int i=iis;i<=iit;i++) 
    for(int j=jjs;j<=jjt;j++) 
    { 
 
        for(int dir=0;dir<DIR;dir++) 
            sum_sq += DQ[i][j][dir] * DQ[i][j][dir]; 
    } 

    if(myrank==0)
       for(int rank=1;rank<cpusize;rank++)
       {
           double tt;
           MPI_Status status;
           MPI_Recv( &tt, 1, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &status );
           sum_sq += tt;
       }
    else
           MPI_Send( &sum_sq, 1, MPI_DOUBLE, 0,    myrank, MPI_COMM_WORLD);
           
             
    if(myrank==0) 
    {
       sum_sq = sqrt(sum_sq/(IDIM*JDIM))/dt; 
       return sum_sq;
    }
    else
       return 1.0E0; 
  
} 
//////////////////////////////////////////////////////////////////////////////// 
/* 
Compute DQ 
*/ 
void convect(const double dt, const double c_alpha) 
{ 
	int iis,iit,jjs,jjt; 
	my_range(myrank,iis,iit,jjs,jjt); 
	 
    //  reset DQ 
    for(int i=iis;i<=iit;i++) 
    for(int j=jjs;j<=jjt;j++) 
    for(int dir=0;dir<DIR;dir++) 
        DQ[i][j][dir] = 0.0E0; 
         
    for(int dir=1;dir<DIR;dir++)   // no dir=0 
    {         
        //1 j direction 
        for(int j=jjs;j<=jjt;j++) 
        { 
            clear_space(); 
            muscli(j,dir,bcdat,nbcs);  // compute slope and limited interface values 
            fluxi(j,dir);              // compute interface flux              
             
        } 
    } 
     
     
    for(int dir=1;dir<DIR;dir++)   // include rest particles 
    { 
        //2 j direction 
        for(int i=iis;i<=iit;i++) 
        { 
            clear_space(); 
            musclj(i,dir,bcdat,nbcs); 
            fluxj(i,dir); 
        } 
    } 
     
    for(int dir=1;dir<DIR;dir++)          
        //3 volunm average 
        for(int i=iis;i<=iit;i++) 
        for(int j=jjs;j<=jjt;j++) 
            DQ[i][j][dir] *= (dt/VOL[i][j]);  // df = flux*dt/V 
         
         
 
         
    return; 
} 
 
//////////////////////////////////////////////////////////////////////////////// 
/* 
Sum of Q and DQ 
reset DQ 
Update M 
*/ 
void update(void) 
{ 
    int ne=0; 
    double TQ[DIR]; 
    int iis,iit,jjs,jjt; 
	  my_range(myrank,iis,iit,jjs,jjt); 
     
    for(int i=iis;i<=iit;i++) 
    for(int j=jjs;j<=jjt;j++) 
    { 
        for(int dir=0;dir<DIR;dir++) 
        { 
             TQ[dir] = Q[i][j][dir] + DQ[i][j][dir]; 
             //DQ[i][j][dir] = 0.0E0;   // reset DQ 
        } 
	    flag[i][j]=0; 
         
        macro(TQ, M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
         
        if(M[i][j][3]<0.0E0) //neg energy, don't altering Q and restore M 
        { 
             //std::cout<<"Negtive energy at ["<<i<<", "<<j<<"]"<<std::endl; 
             ne++; 
             flag[i][j] = 1; 
             macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
        } 
        else  // positive enrgy, OK, update Q and conserve M 
        { 
            for(int dir=0;dir<DIR;dir++) 
                Q[i][j][dir] = TQ[dir]; 
        } 
         
    } 
    //if(ne>0) std::cout<<ne<<" cells were limited because of negtive internal energy" 
    //                  <<std::endl; 
    return; 
} 
 
 
 
void clear_space(void) 
{ 
    for(int ii=0;ii<MAXDIM+3;ii++) 
    { 
        DQP[ii] = 0.0E0; 
        DQM[ii] = 0.0E0; 
         QP[ii] = 0.0E0; 
         QM[ii] = 0.0E0; 
    } 
     
    return; 
} 
 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
/* 
Compute DQ 
Limit DQ to compute DQP, DQM 
Extrapolate QM, QP 
 
bcdat[ibc][0]: BC type 
          [1]: BC direction 
          [2]: i(1 or IDIM) 
          [3]: istrt 
          [4]: iend 
          [5]: jstrt 
          [6]: jend 
          [7]: profile index 
 
*/ 
void muscli(const int cj, const int cdir, 
           const int bcdat[][8], 
           const int nbcs 
           ) 
{ 
int iis,iit,jjs,jjt; 
my_range(myrank,iis,iit,jjs,jjt); 
 
/** 
1.a  Compute DQ and store it in DQP[0..IDIM+1] 
*/     
for(int i=iis-2;i<=iit+1;i++) 
        DQP[i] = Q[i+1][cj][cdir] - Q[i][cj][cdir]; 
 
/** 
1.b  set delta(q) to 2*delta(q) to adjust for the wall boundary being on the 
 wall rather than at a cell center 
*/ 
for(int ibc=0;ibc<nbcs;ibc++) 
{ 
    int bctype = bcdat[ibc][0]; 
    int idir = bcdat[ibc][1]; 
    if(idir==0 && (bctype<200 && bctype>100)) 
    { 
        int in_out = bcdat[ibc][2]; 
        int jstrt = bcdat[ibc][5]; 
        int jend  = bcdat[ibc][6]; 
        if(cj>=jstrt && cj<=jend)  
             DQP[in_out] += DQP[in_out];         
    }     
} 
 
/** 
2.a limit delta(q) 
*/  
for(int i=iis-1;i<=iit+1;i++) 
{ 
    double dp,dm; 
    limiter(1.0E0/3.0E0, DQP[i-1], DQP[i], dp, dm, cj); 
    DQP[i-1] = dp; 
    DQM[i  ] = dm; 
} 
 
/** 
3.a construct Q(+,-) 
*/ 
for(int i=iis-1;i<=iit;i++) 
{ 
    QP[i] = Q[i+1][cj][cdir] + DQP[i]; 
    QM[i] = Q[i  ][cj][cdir] + DQM[i]; 
} 
 
/** 
3.b Correct wall value of Q(+,-) 
*/ 
for(int ibc=0;ibc<nbcs;ibc++) 
{ 
    int bctype = bcdat[ibc][0]; 
    int idir = bcdat[ibc][1]; 
    if(idir==0 && (bctype<200 && bctype>100)) 
    { 
        int in_out = bcdat[ibc][2]; 
        int jstrt = bcdat[ibc][5]; 
        int jend  = bcdat[ibc][6]; 
        if(cj>=jstrt && cj<=jend)  
        { 
           int ibcpt; 
           if(in_out==1)  
              ibcpt = 1; 
           else  
              ibcpt = IDIM+1; 
           // i=1,IDIM+1 is wall cell 
           QP[in_out] = Q[ibcpt][cj][cdir];  
           QM[in_out] = Q[ibcpt][cj][cdir]; 
        } 
    }     
} 
// Finished with Q(+,-) caculation for muscl in i 
 
 
return;     
} 
 
/** 
Compute flux on all interfaces along I direction and store them in FLUX 
DQP and DQM should be ready 
 
Exact Riemann solver is used to compute flux because the convection term is  
linear. 
*/ 
void fluxi(const int cj, const int cdir) 
{ 
    double UBAR,FLUX; 
    int iis,iit,jjs,jjt; 
	  my_range(myrank,iis,iit,jjs,jjt); 
     
    for(int i=iis-1;i<=iit;i++) // loop all interfaces 
    { 
        UBAR = dc[cdir][0]*S[i][cj][0][0] + dc[cdir][1]*S[i][cj][0][1]; 
         
        FLUX = (UBAR>0.0E0?QM[i]:QP[i]) * UBAR * S[i][cj][0][2]; 
        DQ[i  ][cj][cdir] -= FLUX; 
        DQ[i+1][cj][cdir] += FLUX; 
    } 
     
    return; 
} 
 
//////////////////////////////////////////////////////////////////////////////// 
/* 
Compute DQ 
Limit DQ to compute DQP, DQM 
Extrapolate QM, QP 
*/ 
void musclj(const int ci, const int cdir, 
           const int bcdat[][8], 
           const int nbcs 
           ) 
{ 
int iis,iit,jjs,jjt; 
my_range(myrank,iis,iit,jjs,jjt); 
	 
/** 
1.a  Compute DQ and store it in DQP[] 
*/     
for(int j=jjs-2;j<=jjt+1;j++) 
        DQP[j] = Q[ci][j+1][cdir] - Q[ci][j][cdir]; 
 
/** 
1.b  set delta(q) to 2*delta(q) to adjust for face-center boundaries being on the 
 wall rather than at a cell center 
*/ 
for(int ibc=0;ibc<nbcs;ibc++) 
{ 
    int bctype = bcdat[ibc][0]; 
    int idir = bcdat[ibc][1]; 
    if(idir==1 && (bctype<200 && bctype>100)) 
    { 
        int in_out = bcdat[ibc][2]; 
        int istrt = bcdat[ibc][3]; 
        int iend  = bcdat[ibc][4]; 
        if(ci>=istrt && ci<=iend)  
             DQP[in_out] += DQP[in_out];         
    }     
} 
 
 
 
/** 
2.a limit delta(q) 
*/ 
for(int j=jjs-1;j<=jjt+1;j++) 
{ 
    double dp,dm; 
    limiter(1.0E0/3.0E0, DQP[j-1], DQP[j], dp, dm, ci); 
    DQP[j-1] = dp; 
    DQM[j  ] = dm; 
} 
 
/** 
3.a construct Q(+,-) at face 
*/ 
for(int j=jjs-1;j<=jjt;j++) 
{ 
    QP[j] = Q[ci][j+1][cdir] + DQP[j]; 
    QM[j] = Q[ci][j  ][cdir] + DQM[j]; 
} 
 
/** 
3.b Correct wall value of Q(+,-) 
*/ 
for(int ibc=0;ibc<nbcs;ibc++) 
{ 
    int bctype = bcdat[ibc][0]; 
    int idir = bcdat[ibc][1]; 
    if(idir==1 && (bctype<200 && bctype>100)) 
    { 
        int in_out = bcdat[ibc][2]; 
        int istrt = bcdat[ibc][3]; 
        int iend  = bcdat[ibc][4]; 
        if(ci>=istrt && ci<=iend)  
        { 
           int ibcpt; 
           if(in_out==1)  
              ibcpt = 1; 
           else  
              ibcpt = JDIM+1; 
           QP[in_out] = Q[ci][ibcpt][cdir]; 
           QM[in_out] = Q[ci][ibcpt][cdir]; 
        } 
    }     
} 
// Finished with Q(+,-) caculation for muscl in i 
 
 
return;     
} 
 
/** 
Limiters of MUSCL shceme. Only smooth limiter was available 
*/ 
void limiter(  const double rkap, 
               const double vm, const double vp, 
               double &dp, double &dm, const int index) 
{ 
    double sm; 
    sm  = smooth(vm,vp); 
     
    dp =  -0.25E0 * sm *  
        ( (1.0E0-sm*rkap)*vp + 
          (1.0E0+sm*rkap)*vm ); 
    dm =   0.25E0 * sm *  
        ( (1.0E0-sm*rkap)*vm + 
          (1.0E0+sm*rkap)*vp ); 
           
    return ; 
} 
 

 
/** 
Compute flux on a colum of interfaces along J direction and store them in FLUX 
DQP and DQM should be ready 
 
Exact Riemann solver is used to compute flux because the convection term is  
linear. 
*/ 
void fluxj(const int ci, const int cdir) 
{ 
    double VBAR, FLUX; 
    int iis,iit,jjs,jjt; 
    my_range(myrank,iis,iit,jjs,jjt); 
     
    for(int j=jjs-1;j<=jjt;j++) // loop all interfaces 
    { 
        VBAR = dc[cdir][0]*S[ci][j][1][0] + dc[cdir][1]*S[ci][j][1][1]; 
 
        FLUX = (VBAR>0.0E0?QM[j]:QP[j])*VBAR * S[ci][j][1][2]; 
        DQ[ci][j  ][cdir] -= FLUX; 
        DQ[ci][j+1][cdir] += FLUX; 
    }     
    return; 
} 
 
 
void init(const double c_rho,  
          const double c_u, 
          const double c_v, 
          const double c_e) // Uniform initialization 
{ 
    double irho, iu, iv, ie; 
    double rr,rU,re,rp;     
    get_ref(rr,rU,re,rp);  
     
    irho = c_rho/rr; 
    iu = c_u/rU; 
    iv = c_v/rU; 
    ie = c_e/re; 
     
    double c_feq[DIR]; 
     
    get_feq(irho,iu,iv,ie,c_feq); 
     
    for(int i=2;i<=IDIM;i++) 
    for(int j=2;j<=JDIM;j++) 
    { 
        flag[i][j] = 0; 
        M[i][j][0] = irho; 
        M[i][j][1] = iu; 
        M[i][j][2] = iv; 
        M[i][j][3] = ie; 
         
        for(int dir=0;dir<DIR;dir++) 
        { 
            Q[i][j][dir] = c_feq[dir]; 
            DQ[i][j][dir] = 0.0E0; 
        } 
    } 
    return;     
} 
 
void usr_init(void) 
{ 
    double rr,rU,re,rp;     
    get_ref(rr,rU,re,rp); 
    
  double dx,dy,dis;
  
	for(int j=2;j<=JDIM;j++) 
	for(int i=2;i<=IDIM;i++) 
	{ 
/*
  // circular shock
		      M[i][j][0]=1.0E0/rr;
		      M[i][j][1]=0.0E0/rU;
		      M[i][j][2]=0.0E0;
		      M[i][j][3]=(1.0E0/(0.4E0*1.0E0))/re;

		dx = RC[i][j][0] - 0.3E0;
		dy = RC[i][j][1] - 0.7E0;
		dis= sqrt(dx*dx+dy*dy);
				
		if(dis<0.2E0)
			{
				M[i][j][0]=5.0E0/rr;
				M[i][j][1]=0.0E0;
				M[i][j][2]=0.0E0;
				M[i][j][3]=(5.0E0/(0.4E0*5.0E0))/re;
			}
*/
/*
// JCP shock
                     M[i][j][0]=1.0E0/rr;
                     M[i][j][1]=0.0E0;
                     M[i][j][2]=0.0E0;
                     M[i][j][3]=(1.0E0/(0.4E0*1.0E0))/re;

                     if(RC[i][j][1]+RC[i][j][0]-0.15E0<0.0E0)
                     {
                            M[i][j][0]= 0.125E0/rr;
                            M[i][j][1]=0.0E0;
                            M[i][j][2]=0.0E0;
                            M[i][j][3]=(0.14E0/(0.4E0*0.125E0))/re;
                     }
//end JCP shock
*/

         get_anM(0.0E0, RC[i][j][0], RC[i][j][1], M[i][j]);
         M[i][j][0] /= rr;
         M[i][j][1] /= rU;
         M[i][j][2] /= rU;
         M[i][j][3] /= re;


	  get_feq(M[i][j][0],M[i][j][1],M[i][j][2],M[i][j][3], Q[i][j]); 
	  for(int dir=0;dir<DIR;dir++)  
		    DQ[i][j][dir] = 0.0E0; 
	 
	} 
 
	return;  
} 



void get_anM(const double ct, const double cx, const double cy, double* cruve)
{
double cr,cd,crho;
cr = sqrt(cx*cx+cy*cy);
cruve[0] = 1.0E0 + ct/cr;
cruve[1] = -cx/cr;
cruve[2] = -cy/cr;
cruve[3] = 1.0E-6/(Gamma_1*cruve[0]);

return;
}

 
/* 
const int field: output mode 
0   Write macro vars with cell-center format 
1   Write macro vars with grid-node format 
*/ 
int write_plt(std::string& prefix, const double c_time,const int c_count, const int field) 
{ 
    // only CPU can write PLT file when the whole Q is available on it 
    if(myrank!=0) 
    {
    	//std::cout<<"CPU "<<myrank<<": not my bussiness"<<std::endl;
    	return 0; 
    }
    
    if(get_GATHR()==0) 
    {
    	//std::cout<<"Kao, not gather yet"<<std::endl;
    		return 0;
    }
    	 
    // Compute all macro-states 
    if(get_ALL_UPDATED()==0) 
    { 
       std::cout<<"Kao, I have to compute macro"<<std::endl;
       for(int i=1;i<=IDIM+1;i++) 
       for(int j=1;j<=JDIM+1;j++)
       { 
          macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
          //if(j==1) std::cout<<M[i][j][1]<<std::endl;
       }
    	 
    	 set_ALL_UPDATED(); 
    } 
     
    std::ostringstream  plt_name; 
    plt_name<<prefix<<"_"<<std::setw(6)<<std::setfill('0')<<c_count 
            <<".plt"; 
    std::string filename= plt_name.str(); 
    //std::cout<<"PLT file name is "<<filename.c_str()<<std::endl; 
    std::ofstream plt_file(filename.c_str(),std::ios_base::out); 
    if(!plt_file) 
    { 
        std::cerr<<"Can not create a new plt file"<<std::endl; 
        return -1; 
    } 
     
     
     
     
    plt_file.precision(11); 
    if(field==0)  // Write macro vars with cell-center format 
    { 
        std::cout<<"Output tecplot file in cell center format"<<std::endl;
        // extraplate values on the cornner 
        for(int iv=0;iv<4;iv++) 
        { 
            // i=1,j=1 
            M[1][1][iv] =  0.5E0*(M[2][1][iv]*1.5E0-0.5E0*M[3][1][iv])   // row 
                         + 0.5E0*(M[1][2][iv]*1.5E0-0.5E0*M[1][3][iv]);  // col 
            // i=1, j=JDIM+1 
            M[1][JDIM+1][iv] =  0.5E0*(M[2][JDIM+1][iv]*1.5E0-0.5E0*M[3][JDIM+1][iv])  // row 
                         + 0.5E0*(M[1][JDIM  ][iv]*1.5E0-0.5E0*M[1][JDIM-1][iv]);   // col 
            // i=IDIM+1, j=1 
            M[IDIM+1][1][iv] =  0.5E0*(M[IDIM  ][1][iv]*1.5E0-0.5E0*M[IDIM-1][1][iv])   // row 
                              + 0.5E0*(M[IDIM+1][2][iv]*1.5E0-0.5E0*M[IDIM+1][3][iv]);  // col 
            // i=IDIM+1, j=JDIM+1 
            M[IDIM+1][JDIM+1][iv] =  0.5E0*(M[IDIM  ][JDIM+1][iv]*1.5E0-0.5E0*M[IDIM-1][JDIM+1][iv])   // row 
                                   + 0.5E0*(M[IDIM+1][JDIM  ][iv]*1.5E0-0.5E0*M[IDIM+1][JDIM-1][iv]);  // col 
        } 
 
        plt_file<<"variables=\"x\",\"y\",\"Rho\",\"U\",\"V\",\"e\",\"Flag\" \n" 
        <<" ZONE T=\""<<c_time*ref_L/ref_U<<"\", I="<<IDIM+1<<", J="<<JDIM+1<<", F=POINT"<<std::endl; 
        for(int j=1;j<=JDIM+1;j++) 
        for(int i=1;i<=IDIM+1;i++) 
        { 
            double cM[5]; 
             
            if(i==1 && j==1) 
                for(int iv=0;iv<NVAR;iv++) 
                cM[iv] =  0.5E0*(M[2][1][iv]*1.5E0 - 0.5E0*M[3][1][iv])   // row 
                        + 0.5E0*(M[1][2][iv]*1.5E0 - 0.5E0*M[1][3][iv]);  // col 
                         
            else if(i==1 && j==JDIM+1) 
                for(int iv=0;iv<NVAR;iv++) 
                cM[iv] =0.5E0*(M[2][JDIM+1][iv]*1.5E0 - 0.5E0*M[3][JDIM+1][iv])  // row 
                        + 0.5E0*(M[1][JDIM  ][iv]*1.5E0 - 0.5E0*M[1][JDIM-1][iv]);   // col 
                 
            else if(i==IDIM+1 && j==1) 
                for(int iv=0;iv<NVAR;iv++) 
                cM[iv] =  0.5E0*(M[IDIM  ][1][iv]*1.5E0 - 0.5E0*M[IDIM-1][1][iv])   // row 
                        + 0.5E0*(M[IDIM+1][2][iv]*1.5E0 - 0.5E0*M[IDIM+1][3][iv]);  // col 
                         
            else if(i==IDIM+1 && j==JDIM+1) 
                for(int iv=0;iv<NVAR;iv++) 
                cM[iv] =  0.5E0*(M[IDIM  ][JDIM+1][iv]*1.5E0 - 0.5E0*M[IDIM-1][JDIM+1][iv])   // row 
                        + 0.5E0*(M[IDIM+1][JDIM  ][iv]*1.5E0 - 0.5E0*M[IDIM+1][JDIM-1][iv]);  // col 
                         
            else 
                for(int iv=0;iv<NVAR;iv++) 
                    cM[iv] = M[i][j][iv]; 
 
        plt_file<<RC[i][j][0]<<" "<<RC[i][j][1]<<" " 
                <<cM[0]*ref_rho<<"   " 
                <<cM[1]*ref_U<<"   " 
                <<cM[2]*ref_U<<"   " 
                <<cM[3]*ref_e<<"   " 
                <<flag[i][j]<<"    "<<std::endl;; 
        //if(NVAR==5) plt_file<<cM[4]<<std::endl; 
	      //else  
	      //plt_file<<std::endl; 
        } 
    } 
    else if(field==1)   // grid format 
    { 
         
        plt_file<<"variables=\"x\",\"y\",\"Rho\",\"U\",\"V\",\"e\"\n" 
        <<" ZONE T=\""<<c_time<<"s\", I="<<IDIM<<", J="<<JDIM<<", F=POINT"<<std::endl; 
        // compute all nodes on boundary edges except corners 
        double east[JDIM+3][NVAR], west[JDIM+3][NVAR],south[IDIM+3][NVAR], north[IDIM+3][NVAR]; 
        get_bc_macro(east,west,south,north); 
 
        for(int j=1;j<=JDIM;j++) 
        for(int i=1;i<=IDIM;i++) 
        { 
            double cM[4]; 
            if(i==1) 
               for(int iv=0;iv<NVAR;iv++)  cM[iv] = east[j][iv]; 
            else if(j==1) 
               for(int iv=0;iv<NVAR;iv++)  cM[iv] = south[i][iv]; 
            else if(i==IDIM) 
               for(int iv=0;iv<NVAR;iv++)  cM[iv] = west[j][iv]; 
            else if(j==JDIM) 
               for(int iv=0;iv<NVAR;iv++)  cM[iv] = north[i][iv]; 
            else 
                for(int iv=0;iv<4;iv++)  
                cM[iv] = 0.25E0*(M[i][j][iv] + M[i+1][j][iv] + M[i][j+1][iv] + M[i+1][j+1][iv]); 
             
        plt_file<<R[i][j][0]<<" "<<R[i][j][1]<<" " 
                <<cM[0]*ref_rho<<"   " 
                <<cM[1]*ref_U<<"   " 
                <<cM[2]*ref_U<<"   " 
                <<cM[3]*ref_e<<"   " 
                <<std::endl;             
                        
        } 
    } 
 
    plt_file.close(); 
     
#ifndef WIN32 
    std::ostringstream  tar_cmd, rm_cmd; 
    tar_cmd<<"tar -czf "<<plt_name.str()<<".tar.gz "<<plt_name.str(); 
    filename = tar_cmd.str(); 
    rm_cmd<<"rm -f "<<plt_name.str(); 
    system(filename.c_str()); 
    filename = rm_cmd.str(); 
    system(filename.c_str()); 
#endif 
   
     
    return 1; 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
int read_init(const std::string& init_f) 
{ 
    std::ifstream init_file(init_f.c_str(),std::ios::in); 
    if(!init_file) 
    { 
        std::cerr<<"Can not open the init file"<<std::endl; 
        return -1; 
    } 
     
    // Assume that only rho, u, v and e in the file and the datas are arranged 
    // as block format 
    int idim,jdim,nv; 
    init_file>>idim>>jdim>>nv; 
    if(idim!=IDIM || jdim!=JDIM || nv!=4) 
    { 
        std::cerr<<"In init file, idim!=IDIM || jdim!=JDIM || nv!=4"<<std::endl; 
        return -2; 
    }  
     
    // read values on grid nodes 
    for(int k=0;k<4;k++) 
    for(int j=1;j<=JDIM;j++) 
    for(int i=1;i<=IDIM;i++) 
        init_file>>Q[i][j][k];  //  Q array as temp 
     
    // compute values in cells 
    for(int k=0;k<4;k++) 
    for(int j=1;j<=JDIM-1;j++) 
    for(int i=1;i<=IDIM-1;i++) 
        M[i+1][j+1][k] = (Q[i][j][k] + Q[i+1][j][k] +  
                          Q[i][j+1][k] + Q[i+1][j+1][k])*0.25E0; 
     
    // compute feq as initial values of f 
    for(int j=2;j<=JDIM;j++) 
    for(int i=2;i<=IDIM;i++) 
        get_feq(M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3], Q[i][j]); 
         
    return 1; 
         
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
int wr_check_2_0(std::string& prefix_or_filename, 
                      double& c_time, 
                      int& c_count, 
                      const int write_or_read 
                      ) 
{ 
    if(write_or_read==0)     // write 
    { 
        // only CPU 0 can write checkoint file when the whole Q is available on it 
	      if(myrank!=0) return 0; 
	      if(get_GATHR()==0) 
	      	{
	      		 //std::cerr<<"Kao! You have to gather first"<<std::endl;
	      		 	return 0;
	        }
	  	 
        std::ostringstream  rst_name; 
        rst_name<<prefix_or_filename<<"_"<<std::setw(6)<<std::setfill('0')<<c_count 
                <<".rst"; 
        std::string filename= rst_name.str(); 
        //	std::cout<<"Checkpoint file name is "<<filename<<std::endl;
        std::ofstream rst_file(filename.c_str(), 
                               std::ios::out|std::ios::binary); 
        if(!rst_file) 
        { 
            std::cerr<<"Can not create a new restart file"<<std::endl; 
            return -1; 
        } 
        int v1=2,v2=0;  // version 2.0 
        rst_file.write((char*)(&v1), sizeof(v1)); 
        rst_file.write((char*)(&v2), sizeof(v2)); 
         
        rst_file.write((char*)(&c_count), sizeof(c_count)); 
        rst_file.write((char*)(&c_time), sizeof(c_time)); 
        rst_file.write((char*)(&IDIM), sizeof(IDIM)); 
        rst_file.write((char*)(&JDIM), sizeof(JDIM)); 
        rst_file.write((char*)(&DIR), sizeof(DIR)); 
 
        
        rst_file.write((char*)Q, sizeof(Q[0][0][0])*ISIZE*JSIZE*DIR); 
         
        rst_file.close(); 
         
        #ifndef WIN32 
        std::ostringstream  tar_cmd, rm_cmd; 
        tar_cmd<<"tar -czf "<<rst_name.str()<<".tar.gz "<<rst_name.str(); 
        rm_cmd<<"rm -f "<<rst_name.str(); 
        filename = tar_cmd.str(); 
        system(filename.c_str()); 
        filename = rm_cmd.str(); 
        system(filename.c_str()); 
        #endif 
         
        return 1; 
    } 
    else if(write_or_read==1) 
    { 
        std::ifstream rst_file(prefix_or_filename.c_str(),std::ios::binary); 
        if(!rst_file) 
        { 
            std::cerr<<"Can not open a restart file"<<std::endl; 
            return -1; 
        } 
         
        int v1,v2; 
        rst_file.read((char*)(&v1), sizeof(v1)); 
        rst_file.read((char*)(&v2), sizeof(v2)); 
        std::cout<<"Restart file is V"<<v1<<"."<<v2<<" version"<<std::endl; 
 
        if(v1!=2 || v2!=0) 
        { 
            std::cerr<<"The version of this restart file dosen't match this function."<<std::endl; 
            return -2; 
        }         
        double t_time; 
        int t_count; 
        rst_file.read((char*)(&t_count), sizeof(t_count)); 
        rst_file.read((char*)(&t_time), sizeof(t_time)); 
        c_time  = t_time; 
        c_count = t_count; 
        std::cout<<"Restarting from Iteration: "<<c_count<<"       and Time: "<<c_time<<std::endl; 
         
        int idim=0,jdim=0,dir=0; 
        rst_file.read((char*)(&idim), sizeof(idim)); 
        rst_file.read((char*)(&jdim), sizeof(jdim)); 
        rst_file.read((char*)(&dir), sizeof(dir)); 
        std::cout<<"This restart file has a set of data whose size is idim="<<idim<<" jdim="<<jdim<<" dir="<<dir<<std::endl; 
     
        if(idim!=IDIM || jdim!=JDIM || dir!=DIR ) 
        { 
            std::cerr<<"This restart file dosen't match this case."<<std::endl; 
            return -1; 
        } 
     
        rst_file.read((char*)Q, sizeof(Q[0][0][0])*(idim+3)*(jdim+3)*dir); 
         
        rst_file.close(); 
 
        // compute rho/u/v/e from Q 
        for(int i=2;i<=IDIM;i++) 
        for(int j=2;j<=JDIM;j++) 
            macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
         
        std::cout<<myrank<<" CPU reading checkpoint file, OK!"<<std::endl; 
        return 1; 
    } 
     
} 
 
 
 
int compute_force(std::string& prefix,  
	                const int ibc,  
                  const double c_time,  
                  const int c_count, 
                  const int write_file, 
                  double* Cx, 
                  double* Cy) 
                 
{ 
    // only CPU 0 can write checkoint file when the whole Q is available on it 
	  if(myrank!=0) return 0; 
	  	 
	  if(get_GATHR()==0) 
	  { 
	     gather_array(); 
	      
       for(int i=1;i<=IDIM+1;i++) 
       for(int j=1;j<=JDIM+1;j++) 
          macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
    	 
    	 set_ALL_UPDATED(); 
    } 
    else if(get_ALL_UPDATED()==0) 
    { 
       for(int i=1;i<=IDIM+1;i++) 
       for(int j=1;j<=JDIM+1;j++) 
          macro(Q[i][j], M[i][j][0], M[i][j][1], M[i][j][2], M[i][j][3]); 
    	 
    	 set_ALL_UPDATED(); 
    } 
     
     
    std::ostringstream  frc_name; 
    frc_name<<prefix<<"_force_"<<std::setw(6)<<std::setfill('0')<<ibc<<"_"<<c_count 
                <<".frc"; 
    std::string filename=frc_name.str(); 
    std::ofstream frc_file; 
 
     if(write_file) 
     { 
        std::cout<<"open file to write force"<<std::endl; 
        frc_file.open(filename.c_str(),std::ios_base::out); 
        if(!frc_file) 
        { 
            std::cerr<<"Can not create a new force file"<<std::endl; 
            return -1; 
        } 
    } 
 
 
    double Fx=0.0,Fy=0.0,s; 
    double trr,trU,trL,tre,trF; 
    get_ref(trr, trU, tre, trL); 
 
     
    if(bcdat[ibc][1]==0)   // I direction 
    { 
        int jstrt,jend,inout; 
        jstrt = bcdat[ibc][5]; 
        jend  = bcdat[ibc][6]; 
        inout = bcdat[ibc][2]; 
         
        if(write_file)  
        	   frc_file<<"variables=\"x\",\"y\",\"nx\",\"ny\", \"s\", \"tx\", \"ty\",\"fx\",\"fy\",\"fn\",\"ft\"\n" 
                     <<" ZONE T=\""<<c_time<<"s\", I="<<1<<", J="<<jend-jstrt+1<<", F=POINT"<<std::endl; 
         
        for(int j=jstrt;j<=jend;j++) 
        { 
            double nx,ny,tx,ty; 
            if(inout==1)   // is inward 
            { 
                nx = -S[inout][j][0][0]; 
                ny = -S[inout][j][0][1]; 
            } 
            else if(inout==IDIM)// is outward 
            { 
         
                nx = S[inout][j][0][0]; 
                ny = S[inout][j][0][1]; 
            } 
            s = S[inout][j][0][2]; 
            tx =  ny; 
            ty = -nx; 
             
            double v_stress[3], fn, fnx, fny, ft, ftx,fty, bp,fx,fy, sign; 
            bp = get_pressure(M[inout][j][0], M[inout][j][3]); 
            get_vstress(M[inout][j], Q[inout][j], v_stress); 
            // normal stress 
            fnx = (v_stress[0]-bp)*nx + v_stress[1]*ny; 
            fny = v_stress[1]*nx + (v_stress[2]-bp)*ny; 
            fn = sqrt(fnx*fnx+fny*fny); 
            // shear stress 
            ftx = v_stress[0]*tx + v_stress[1]*ty; 
            fty = v_stress[1]*tx + v_stress[2]*ty; 
            sign = ftx*tx + fty*ty >0.0E0? 1.0E0: -1.0E0;
            ft = -sign * sqrt(ftx*ftx+fty*fty); 
            // x y components 
            fx = fnx+ftx; 
            fy = fny+fty; 
 
            trF = trr*tre;  
            fn *= trF; 
            ft *= trF; 
            fx *= trF; 
            fy *= trF;  
            // force integral 
            Fx += fx*s; 
            Fy += fy*s; 
 
            double bx=RC[inout][j][0]; 
            double by=RC[inout][j][1]; 
             
            if(write_file)  
            frc_file<<" "<<bx<<" "<<by 
                    <<" "<<nx<<" "<<ny<<" "<<S[inout][j][0][2] 
                    <<" "<<tx<<" "<<ty 
                    <<" "<<fx<<" "<<fy<<" "<<fn<<" "<<ft<<std::endl;  
        } 
    } 
    else if(bcdat[ibc][1]==1)   // J direction 
    { 
        int istrt,iend,inout; 
        istrt = bcdat[ibc][3]; 
        iend  = bcdat[ibc][4]; 
        inout = bcdat[ibc][2]; 
        //std::cout<<"Output boundary force from ["<<istrt<<"]["<<inout<<"] to ["<<iend<<"]["<<inout<<"]"<<std::endl; 
         
        if(write_file)  
        frc_file<<"variables=\"x\",\"y\",\"nx\",\"ny\",\"sx\", \"tx\", \"ty\",\"fx\",\"fy\",\"fn\",\"ft\"\n" 
        <<" ZONE T=\""<<c_time<<"s\", I="<<iend-istrt+1<<", J="<<1<<", F=POINT"<<std::endl; 
         
        for(int i=istrt;i<=iend;i++) 
        { 
            double nx,ny,tx,ty; 
            if(inout==1)   // is inward 
            { 
                nx = -S[i][inout][1][0]; 
                ny = -S[i][inout][1][1]; 
            } 
            else if(inout==JDIM)// is outward 
            { 
         
                nx = S[i][inout][1][0]; 
                ny = S[i][inout][1][1]; 
            } 
            s =  S[i][inout][1][2]; 
            tx =  ny; 
            ty = -nx; 
             
            double v_stress[3], fn, fnx, fny, ft, ftx,fty, bp,fx,fy; 
            bp = get_pressure(M[i][inout][0], M[i][inout][3]); 
            get_vstress(M[i][inout], Q[i][inout], v_stress); 
            // normal stress on the face 
            fnx = (v_stress[0]-bp)*nx + v_stress[1]*ny; 
            fny = v_stress[1]*nx + (v_stress[2]-bp)*ny; 
            fn = sqrt(fnx*fnx+fny*fny); 
            // shear  stress on the face 
            ftx = v_stress[0]*tx + v_stress[1]*ty; 
            fty = v_stress[1]*tx + v_stress[2]*ty; 
            ft = sqrt(ftx*ftx+fty*fty)*(tx>0.0E0?1.0E0:-1.0E0); 
            // x y components 
            fx = fnx + ftx; 
            fy = fny + fty;
 
            trF = trr*tre;  
            fn *= trF; 
            ft *= trF; 
            fx *= trF; 
            fy *= trF;  
            // force integral  
            Fx += fx*s; 
            Fy += fy*s; 
             
            double bx=RC[i][inout][0]; 
            double by=RC[i][inout][1]; 
             
            if(write_file)  
            frc_file<<" "<<bx<<" "<<by 
                    <<" "<<nx<<" "<<ny<<" "<<S[i][inout][1][2] 
                    <<" "<<tx<<" "<<ty 
                    <<" "<<fx<<" "<<fy<<" "<<fn<<" "<<ft<<std::endl;  
        } 
    } 
    double Vr, VL, VU, cx, cy; 
    get_ReV(Vr, VL, VU); 
    // force coeeficients 
    cx = Fx/(0.5E0 * Vr * VU * VU); 
    cy = Fy/(0.5E0 * Vr * VU * VU); 
 
    if(write_file)  
    { 
    frc_file<<"# on this boundary, C_x = "<<cx<<"  C_y = "<<cy<<std::endl; 
    frc_file.close(); 
    } 
 
    if(Cx!=NULL && Cy!=NULL)  
    {*Cx = cx;  *Cy = cy;} 
    return 1; 
} 
 
 
 
