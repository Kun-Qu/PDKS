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



#include "bc.h"
#include "parallel.h"
#include "IO/mutator.h"
#include "IO/control-device.h"

int nbcs;   // Number of bc


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
int bcdat[MAXBC][8];
double bcm[MAXBC][4];  // boundary value in rho, u, v, e

void read_bc(const int n_bcs, const std::string& bc_f)
{
    nbcs = n_bcs;
    ControlDevice Ctrl( GetFileControlDevice(bc_f.c_str(),"BC") );
    ostringstream* bc_entry = new ostringstream[n_bcs];
    ControlDevice* bc_ctrl = new ControlDevice[n_bcs];
    for(int i=0;i<n_bcs;i++)
    {
          bc_entry[i]<<"bc"<<i+1;
          bc_ctrl[i] = Ctrl.getSubDevice( bc_entry[i].str() );

          bc_ctrl[i].add(  "bctype", GetMutator(bcdat[i][0])  );
          bc_ctrl[i].add(     "dir", GetMutator(bcdat[i][1])  );
          bc_ctrl[i].add(   "index", GetMutator(bcdat[i][2])  );
          bc_ctrl[i].add( "start_i", GetMutator(bcdat[i][3])  );
          bc_ctrl[i].add(   "end_i", GetMutator(bcdat[i][4])  );
          bc_ctrl[i].add( "start_j", GetMutator(bcdat[i][5])  );
          bc_ctrl[i].add(   "end_j", GetMutator(bcdat[i][6])  );
          bc_ctrl[i].add(     "mis", GetMutator(bcdat[i][7])  );

          bc_ctrl[i].add(  "brho", GetMutator(bcm[i][0])  );
          bc_ctrl[i].add(    "bU", GetMutator(bcm[i][1])  );
          bc_ctrl[i].add(    "bV", GetMutator(bcm[i][2])  );
          bc_ctrl[i].add(    "be", GetMutator(bcm[i][3])  );
    }

    Ctrl.update();
    //std::cout<<"BC read OK!!"<<std::endl;

    double rr,rU,re,rp;    
    get_ref(rr,rU,re,rp); 
    for(int ibc=0;ibc<n_bcs;ibc++)
    {
      bcm[ibc][0] /= rr;
      bcm[ibc][1] /= rU;
      bcm[ibc][2] /= rU;
      if(bcdat[ibc][0]==204 || bcdat[ibc][0]==202)
      	   bcm[ibc][3] /= rr*re;   // pressure not energy
      else 
      	   bcm[ibc][3] /= re;
    }
      
    //Ctrl.print_values(cout);
    delete [] bc_ctrl;
    delete [] bc_entry;
    return;
}




void metbc(void)
{
    for(int ibc=0;ibc<nbcs;ibc++)
    {
        int bctype = bcdat[ibc][0];
        if(bctype == ISPOUT || bctype == ISPIN || bctype == ISYM ||
           bctype == IFAR)  // no need for S
           {
                int idir, in_out;
                idir = bcdat[ibc][1];
                in_out = bcdat[ibc][2];

                if(idir==0)   // I
                {
                    int jstrt,jend,gv,in1;                    
                    jstrt = bcdat[ibc][5];
                    jend =  bcdat[ibc][6];                    
                    gv  = in_out==1?1:IDIM+1;
                    in1 = in_out==1?2:IDIM;
                    for(int j=jstrt;j<=jend;j++) VOL[gv][j] = VOL[in1][j];
                }
                else   // J
                {
                    int istrt,iend,gv,in1;
                    istrt = bcdat[ibc][3];
                    iend =  bcdat[ibc][4];
                    gv  = in_out==1?1:JDIM+1;
                    in1 = in_out==1?2:JDIM;
                    for(int i=istrt;i<=iend;i++) VOL[i][gv] = VOL[i][in1];
                }

           }
    }// end for(ibc)
}





int set_bcdat(const int ibc, const int ibcpar[8], const double dbcpar[4])
{
    double rr,rU,re,rp;    
    get_ref(rr,rU,re,rp);  
      
    
    for(int i=0;i<8;i++)
        bcdat[ibc][i] = ibcpar[i];
    
//    if(ibcpar[0]==ISPIN)  //
//    {
      bcm[ibc][0] = dbcpar[0]/rr;
      bcm[ibc][1] = dbcpar[1]/rU;
      bcm[ibc][2] = dbcpar[2]/rU;
      bcm[ibc][3] = dbcpar[3]/re;
      
      if(bcm[ibc][3]<0.0E0)
      {
          std::cerr<<"Set_bcdat: negtive boundary energy for supersonic inflow"<<std::endl;
          return 0;
      }
//    }
    
    return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*
This function compute all variables on the boundaries of a grid and return them
in the four arrays, east/west/south/north
*/
void get_bc_macro(double east[JDIM+3][NVAR], 
                  double west[JDIM+3][NVAR],
                  double south[IDIM+3][NVAR],
                  double north[IDIM+3][NVAR])
{
    double c_east[JDIM+3][NVAR], c_west[JDIM+3][NVAR],c_south[IDIM+3][NVAR],c_north[IDIM+3][NVAR];

    for(int ibc=0;ibc<nbcs;ibc++)  // process every boundary
    {
         if(bcdat[ibc][1]==0)  // I dir
         {
             int inout, jstrt,jend,in0,in1;
             jstrt = bcdat[ibc][5];
             jend  = bcdat[ibc][6];
             inout = bcdat[ibc][2];

             if(inout==1)
             {
                in0 = 1;
                in1 = 2;
             }
             else if(inout==IDIM)
             {
                in0   = IDIM+1;
                in1   = IDIM;
             }

             // compute and store
             for(int j=jstrt;j<=jend;j++)
             {
                 double tM[NVAR];
                 // average the ghost and inner cells to get values on the boundary
                 for(int iv=0;iv<NVAR;iv++) 
                 	tM[iv] =  bcdat[ibc][0]>=200?
                                    0.5E0*(M[in0][j][iv] + M[in1][j][iv]):
                                    M[in0][j][iv];

                     if(inout==1)  // east
                         for(int iv=0;iv<NVAR;iv++) c_east[j][iv] = tM[iv];
                     else if(inout==IDIM)  // west
                         for(int iv=0;iv<NVAR;iv++)  c_west[j][iv] = tM[iv];
             }
             
             if(bcdat[ibc][0]==IPERD)
             {
             	   if(inout==1)  // east
             	   	 for(int j=jstrt;j<=jend;j++)
             	   	 for(int iv=0;iv<NVAR;iv++) 
             	              c_west[j][iv] = c_east[j][iv];
             	   else if(inout==IDIM)  // west
             	   	 for(int j=jstrt;j<=jend;j++)
                   for(int iv=0;iv<NVAR;iv++)
                            c_east[j][iv] = c_west[j][iv];
             }
             else if(bcdat[ibc][0]==ICGRID)
             {
             		int Ajstr,Ajend, Bjstr,Bjend;
             		Ajstr = bcdat[ibc][5];
             		Ajend = bcdat[ibc][6];
             		Bjstr = bcdat[ibc][3];
             		Bjend = bcdat[ibc][4];
             		
             		int jl = Ajend - Ajstr + 1;
             		
             		if(inout==1)  // east
             		  for(int j=0;j<=jl;j++)
             		  for(int iv=0;iv<NVAR;iv++) 
             	           c_east[Bjstr-j][iv] = c_east[Ajstr+j][iv];
             	  else if(inout==IDIM)
             	  	for(int j=0;j<=jl;j++)
             	  	for(int iv=0;iv<NVAR;iv++) 
             	  	       c_west[Bjstr-j][iv] = c_west[Ajstr+j][iv];
             }
         }//  endif I dir
         else if(bcdat[ibc][1]==1)  // J dir
         {
             int inout,istrt,iend,in0,in1;
             istrt = bcdat[ibc][3];
             iend  = bcdat[ibc][4];
             inout = bcdat[ibc][2];

             if(inout==1)
             {
                in0 = 1;
                in1 = 2;
             }
             else if(inout==JDIM)
             {
                in0   = JDIM+1;
                in1   = JDIM;
             }

             for(int i=istrt;i<=iend;i++)
             {
                 double tM[NVAR];
                 for(int iv=0;iv<NVAR;iv++) 
                         tM[iv] =  bcdat[ibc][0]>=200?
                                    0.5E0*(M[i][in0][iv] + M[i][in1][iv]):
                                    M[i][in0][iv];

                 if(inout==1)  // south 
                     for(int iv=0;iv<NVAR;iv++) c_south[i][iv] = tM[iv];
                 else if(inout==JDIM)  // north
                     for(int iv=0;iv<NVAR;iv++) c_north[i][iv] = tM[iv];
             }
             
             if(bcdat[ibc][0]==IPERD)
             {
             	   if(inout==1)  // south 
             	   	 for(int i=istrt;i<=iend;i++)
             	   	 for(int iv=0;iv<NVAR;iv++) 
             	              c_north[i][iv] = c_south[i][iv];
             	   else if(inout==JDIM)  // west
             	   	 for(int i=istrt;i<=iend;i++)
             	   	 for(int iv=0;iv<NVAR;iv++) 
                            c_south[i][iv] = c_north[i][iv];
             }
             else if(bcdat[ibc][0]==ICGRID)
             {
             	  //std::cout<<"C-grid"<<std::endl;
             		int Aistr,Aiend, Bistr,Biend;
             		Aistr = bcdat[ibc][3];
             		Aiend = bcdat[ibc][4];
             		Bistr = bcdat[ibc][5];
             		Biend = bcdat[ibc][6];
             		
             		int il = Aiend - Aistr + 1;
             		
             		if(inout==1)  // south 
             		  for(int i=0;i<=il;i++)
             		  for(int iv=0;iv<NVAR;iv++) 
             	           c_south[Bistr-i][iv] = c_south[Aistr+i][iv];
             	  else if(inout==IDIM)
             	  	for(int i=0;i<=il;i++)
             	  	for(int iv=0;iv<NVAR;iv++) 
             	  	       c_north[Bistr-i][iv] = c_north[Aistr+i][iv];
             }
         }//  J dir
    } // end loop ibc for cell type BC


    // computing values on corners
    for(int iv=0;iv<NVAR;iv++)
    {
        c_south[1][iv]= 1.5E0*c_south[2][iv] - 0.5E0*c_south[3][iv];
        c_north[1][iv]= 1.5E0*c_north[2][iv] - 0.5E0*c_north[3][iv];
        c_south[IDIM+1][iv]= 1.5E0*c_south[IDIM][iv] - 0.5E0*c_south[IDIM-1][iv];
        c_north[IDIM+1][iv]= 1.5E0*c_north[IDIM][iv] - 0.5E0*c_north[IDIM-1][iv];

        c_east[1][iv]= 1.5E0*c_east[2][iv] - 0.5E0*c_east[3][iv];
        c_west[1][iv]= 1.5E0*c_west[2][iv] - 0.5E0*c_west[3][iv];
        c_east[JDIM+1][iv]= 1.5E0*c_east[JDIM][iv] - 0.5E0*c_east[JDIM-1][iv];
        c_west[JDIM+1][iv]= 1.5E0*c_west[JDIM][iv] - 0.5E0*c_west[JDIM-1][iv];

        // store in the arrays for returning
        east[1][iv] = 0.5E0*(c_east[1][iv] + c_south[1][iv]);
        south[1][iv] = east[1][iv];
        east[JDIM][iv] = 0.5E0*(c_east[JDIM+1][iv] + c_north[1][iv]);
        north[1][iv] = east[JDIM+1][iv];
        west[1][iv] = 0.5E0*(c_west[1][iv] + c_south[IDIM+1][iv]);
        south[IDIM][iv] = west[1][iv];
        west[JDIM][iv] = 0.5E0*(c_west[JDIM+1][iv] + c_north[IDIM+1][iv]);
        north[IDIM][iv] = west[JDIM+1][iv];
    }

    // compute values on nodes(2..DIM-1) from values in cells
    for(int i=2;i<IDIM;i++)  
         for(int iv=0;iv<NVAR;iv++)
         {
               south[i][iv] = 0.5E0*(c_south[i][iv] + c_south[i+1][iv]); 
               north[i][iv] = 0.5E0*(c_north[i][iv] + c_north[i+1][iv]);
         }

    for(int j=2;j<JDIM;j++) 
         for(int iv=0;iv<NVAR;iv++)
         {     
               east[j][iv] = 0.5E0*(c_east[j][iv] + c_east[j+1][iv]);
               west[j][iv] = 0.5E0*(c_west[j][iv] + c_west[j+1][iv]);
         }

    return;
}

////////////////////////////////////////////////////////////////////////////////
void bc(double time)
{
    for(int ibc=0;ibc<nbcs;ibc++)
    {
        switch(bcdat[ibc][0])
        {
            case ITAN:
                        //std::cout<<myrank<<" CPU.....  ITAN"<<std::endl;
                        bc_tan(ibc);
                        break;
            case ISYM:
                        //std::cout<<myrank<<" CPU.....  ISYM"<<std::endl;
                        bc_tan(ibc);
                        break;
            case ITAN2:
                        //std::cout<<myrank<<" CPU.....  ITAN2"<<std::endl;
                        bc_tan2(ibc);
                        break;
            case IWALL:
                        //std::cout<<myrank<<" CPU.....  IWALL"<<std::endl;
                        bc_wall(ibc);
                        break;
            case ISUIN:
                        //std::cout<<myrank<<" CPU.....  ISUIN"<<std::endl;
                        bc_suin(ibc);
                        break;
            case ISUOUT:
                        //std::cout<<myrank<<" CPU.....  ISUOUT"<<std::endl;
                        bc_suout(ibc);
                        break;
            case ISPIN:
                        //std::cout<<myrank<<" CPU.....  ISPIN"<<std::endl;
                        bc_spin(ibc );
                        break;
            case ISPOUT:
                        //std::cout<<myrank<<" CPU.....  ISPOUT"<<std::endl;
                        bc_spout(ibc );
                        break;
            case ICGRID:
                        //std::cout<<myrank<<" CPU.....  ICGRID"<<std::endl;
                        bc_cgrid(ibc );
                        break;
            case IFAR:
                        //std::cout<<myrank<<" CPU.....  IFAR"<<std::endl;
                        bc_far(ibc );
                        break;
            case IPERD:
                        //std::cout<<myrank<<" CPU.....  IPERD"<<std::endl;
                        bc_perd(ibc);
                        break;
            case IUSRDEF:
            	          bc_usrdef(ibc,time);
            	          break;
            default:
                        std::cout<<"Unknow BC type: "<<bcdat[ibc][0]<<std::endl;
                        break;
        }
    }

    return;

}

/////////////////////////////////////////////////////////////////////
/*
subsonic inflow with free stream values
give frestream rho, u, v and e of free stream

compute total pressure and total temperature of free stream
extrapolate to get u and v on boundary
compute static temperature from total temperature and u,v
compute Mach number on boundary
compute static pressure from total pressure and Mach number
compute density from temperature and pressure

we get all rho, u, v and e on boundary
    
*/
int bc_suin(const int ibc)
{
#ifdef DEBUG
    std::cout<<"IFSUIN BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    double rfs,ufs,vfs,efs;
    rfs = bcm[ibc][0];
    ufs = bcm[ibc][1];
    vfs = bcm[ibc][2];
    efs = bcm[ibc][3];
    double t_un = ufs*ufs+vfs*vfs;
    double t_ms = t_un/(Gamma*Gamma_1*efs);
    double t_pt = 1+5.0E-1*Gamma_1*t_ms;
    double t_et = t_pt*efs;  // total intenal energy
    t_pt = pow(t_pt, Gamma/Gamma_1);
    t_pt = t_pt * (Gamma_1*efs*rfs);  // total pressure

    t_un = 1.0E0/sqrt(t_un);
    double t_sinV = vfs*t_un;
    double t_cosV = ufs*t_un;

    if(bcdat[ibc][1]==0)  // i dir
    {
        int jstrt,jend,inout,itan,gh,in1;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error subsonic inlet BC, jend<jtrt"<<std::endl;
            return 0;
        }
        
        inout = bcdat[ibc][2];
        for(int j=jstrt;j<=jend;j++)
        {
            // Get normal vector
            double nx,ny,vn,vt;
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
                nx = S[inout][j][0][0];
                ny = S[inout][j][0][1];
            }
            else  
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                itan = IDIM+1;  // bd "cell"
                
                nx = -S[inout][j][0][0];
                ny = -S[inout][j][0][1];
            }// endif inout
           
            double rin,uin,vin,ein; 
            rin = M[in1][j][0];
            uin = M[in1][j][1];
            vin = M[in1][j][2];
            ein = M[in1][j][3];
            
            double QEXT = uin*nx + vin*ny;
            double RHON = rin;
            double PN   = Gamma_1*rin*ein;

            // One Dimensional Riemann invariant
            double AEXT = sqrt(Gamma*PN/RHON);
            double REXT = QEXT - 2.0E0 * AEXT / Gamma_1 ;
            double HTEXT = (AEXT*AEXT/Gamma_1 + 0.5E0*QEXT*QEXT);
            double A = (- REXT 
                        + sqrt(REXT*REXT
                               +(Gamma+1.E0)*(HTEXT-0.5E0*REXT*REXT)) )
                        * Gamma_1 / (Gamma + 1.E0);
            double QDOTN = REXT + 2.0E0 * A / Gamma_1;
            double RMACH = QDOTN / A;
            
            // Isentropic relationships
            double ARG = 1.0E0 + 0.5E0 * Gamma_1 * RMACH * RMACH;
            double B_P = t_pt / pow(ARG, (Gamma/Gamma_1));
            double B_E = t_et / ARG;
            double B_RHO = B_P/(Gamma_1*B_E);
            double B_U = QDOTN * nx;
            double B_V = QDOTN * ny;
            //
            double cfeq[DIR];
            get_feq(B_RHO,B_U,B_V,B_E,cfeq);
            
            // set boundary macro-states just for outputing results
            M[itan][j][0] = M[gh][j][0] = B_RHO;
            M[itan][j][1] = M[gh][j][1] = B_U;
            M[itan][j][2] = M[gh][j][2] = B_V;
            M[itan][j][3] = M[gh][j][3] = B_E;
            //std::cout<<j<<"  @@@  "<<r[1]<<std::endl;
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[itan][j][dir] =  cfeq[dir];
                   Q[gh  ][j][dir] =  cfeq[dir];
            }
            
        }// end loop
        
    }
    else if(bcdat[ibc][1]==1)  // j dir
    {
        int istrt,iend,inout,itan,gh,in1;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error subsonic inlet BC, iend<itrt"<<std::endl;
            return 0;
        }
        
        inout = bcdat[ibc][2];
        for(int i=istrt;i<=iend;i++)
        {
            // Get normal vector
            double nx,ny,vn,vt;
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
                nx = S[i][inout][0][0];
                ny = S[i][inout][0][1];
            }
            else  
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                itan = JDIM+1;  // bd "cell"
                
                nx = -S[i][inout][0][0];
                ny = -S[i][inout][0][1];
            }// endif inout
           
            double rin,uin,vin,ein; 
            rin = M[i][in1][0];
            uin = M[i][in1][1];
            vin = M[i][in1][2];
            ein = M[i][in1][3];
            
            double QEXT = uin*nx + vin*ny;
            double RHON = rin;
            double PN   = Gamma_1*rin*ein;

            // One Dimensional Riemann invariant
            double AEXT = sqrt(Gamma*PN/RHON);
            double REXT = QEXT - 2.0E0 * AEXT / Gamma_1 ;
            double HTEXT = (AEXT*AEXT/Gamma_1 + 0.5E0*QEXT*QEXT);
            double A = (- REXT 
                        + sqrt(REXT*REXT
                               +(Gamma+1.E0)*(HTEXT-0.5E0*REXT*REXT)) )
                        * Gamma_1 / (Gamma + 1.E0);
            double QDOTN = REXT + 2.0E0 * A / Gamma_1;
            double RMACH = QDOTN / A;
            
            // Isentropic relationships
            double ARG = 1.0E0 + 0.5E0 * Gamma_1 * RMACH * RMACH;
            double B_P = t_pt / pow(ARG, (Gamma/Gamma_1) );
            double B_E = t_et / ARG;
            double B_RHO = B_P/(Gamma_1*B_E);
            double B_U = QDOTN * nx;
            double B_V = QDOTN * ny;
            //
            double cfeq[DIR];
            get_feq(B_RHO,B_U,B_V,B_E,cfeq);
            
            // set boundary macro-states just for outputing results
            M[i][itan][0] = M[i][gh][0] = B_RHO;
            M[i][itan][1] = M[i][gh][1] = B_U;
            M[i][itan][2] = M[i][gh][2] = B_V;
            M[i][itan][3] = M[i][gh][3] = B_E;
            //std::cout<<j<<"  @@@  "<<r[1]<<std::endl;
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][itan][dir] =  cfeq[dir];
                   Q[i][gh  ][dir] =  cfeq[dir];
            }
            
        }// end loop
        
    }
    
    return 1;
    
    
}

/////////////////////////////////////////////////////////////////////
/*
far field 
Here, I only implemented non-modified algorithm based on Reimann invariants.
This method is not very good for subsonic far-field boundary. later I will add
modification in order to use a smaller domain.
*/
/////////////////////////////////////////////////////////////////////
int bc_far(const int ibc )
{
    double r[3],u[3],v[3],e[3];
    r[0] = bcm[ibc][0];
    u[0] = bcm[ibc][1];
    v[0] = bcm[ibc][2];
    e[0] = bcm[ibc][3];
#ifdef DEBUG    
    std::cout<<"IFAR BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif    
    if(bcdat[ibc][1]==0)  // i dir
    {
        int jstrt,jend,inout,itan,gh,in1;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }
        
        inout = bcdat[ibc][2];
        for(int j=jstrt;j<=jend;j++)
        {
            // Get normal vector
            double nx,ny,vn,vt;
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
                nx = S[inout][j][0][0];
                ny = S[inout][j][0][1];
            }
            else  
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                itan = IDIM+1;  // bd "cell"
                
                nx = -S[inout][j][0][0];
                ny = -S[inout][j][0][1];
            }// endif inout
            
            r[2] = M[in1][j][0];
            u[2] = M[in1][j][1];
            v[2] = M[in1][j][2];
            e[2] = M[in1][j][3];
            
            // determin boundary value by 1-D Reimann invariants            
            reimann_farfield(r,u,v,e, nx, ny);
            
            double cfeq[DIR];
            get_feq(r[1],u[1],v[1],e[1],cfeq);
            
            // set boundary macro-states just for outputing results
            M[itan][j][0] = r[1];
            M[itan][j][1] = u[1];
            M[itan][j][2] = v[1];
            M[itan][j][3] = e[1];
            //std::cout<<j<<"  @@@  "<<r[1]<<std::endl;
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[itan][j][dir] =  cfeq[dir];
                   Q[gh  ][j][dir] =  cfeq[dir];
            }
            
        }// end loop
        
    }
    else if(bcdat[ibc][1]==1)  // j dir
    {
        int istrt,iend,inout,itan,gh,in1;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error tangent BC i range"<<std::endl;
            return 0;
        }
        
        inout = bcdat[ibc][2];
        for(int i=istrt;i<=iend;i++)
        {
            // get normal vector
            double nx,ny,vn,vt;
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
                nx = S[i][inout][1][0];
                ny = S[i][inout][1][1];
            }
            else  
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                itan = JDIM+1;  // bd "cell"
                
                nx = -S[i][inout][1][0];
                ny = -S[i][inout][1][1];
            }// endif inout
            r[2] = M[i][in1][0];
            u[2] = M[i][in1][1];
            v[2] = M[i][in1][2];
            e[2] = M[i][in1][3]; 
            // determin boundary values by means of 1-D Reimann invariants
            reimann_farfield(r,u,v,e, nx, ny);
            // set boundary macro-states just for outputing results
            M[i][itan][0] = r[1];
            M[i][itan][1] = u[1];
            M[i][itan][2] = v[1];
            M[i][itan][3] = e[1];
            //std::cout<<i<<"  @@@  "<<r[1]<<std::endl;
            
            double cfeq[DIR];
            get_feq(r[1],u[1],v[1],e[1],cfeq);
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][itan][dir] =  cfeq[dir];
                   Q[i][gh  ][dir] =  cfeq[dir];
            }
            
        }// end loop
        
    }
    
    return 1;
}

void reimann_farfield(double rho[3], double u[3], double v[3], double e[3], 
                double nx, double ny)
{
    static int i=2;
    double s[3],vn[3],a[3],R[3],mach;
    
    i++;
    // enthopy
    s[0] = Gamma_1*e[0]*pow(rho[0],-Gamma_1);
    s[2] = Gamma_1*e[2]*pow(rho[2],-Gamma_1);
    
    // normal velocity
    vn[0] = u[0]*nx  + v[0]*ny;
    vn[2] = u[2]*nx  + v[2]*ny;
    
    // sound speed
    a[0] = sqrt(Gamma*Gamma_1*e[0]);
    a[2] = sqrt(Gamma*Gamma_1*e[2]);
    
    // Reimann invariants
    R[0] = vn[0] - (a[0]+a[0])/Gamma_1;
    R[2] = vn[2] + (a[2]+a[2])/Gamma_1;
    
    // normal veclocity and sound speed on boundary
    vn[1] = 5.0E-1*   (R[0]+R[2]);
     a[1] = 0.25E0 * Gamma_1 * (R[2]-R[0]);
    //if(i%4==0)
    //    std::cout<<vn[1]<<"   "<<a[1]<<std::endl;
    
    mach = vn[0]/a[0];
    
    if(mach>1.0E0)  // supersonic
     {
         if(vn[1]>=0.0E0)  // inflow
         {
             rho[1] = rho[0]; u[1] = u[0]; 
             v[1] = v[0]; e[1] = e[0];
         }
         else // outflow
         {
             e[1] = e[2];
             rho[1] = rho[2];
             u[1] = u[2];
             v[1] = v[2];
         }
     }
     else  // subsonic
     {
         if(vn[1]>=0.0E0) // inflow
         { 
             u[1] = vn[1]*nx + (u[0] - vn[0]*nx);
             v[1] = vn[1]*ny + (v[0] - vn[0]*ny);
             e[1] = a[1]*a[1]/(Gamma*Gamma_1);
             s[1] = s[0];
             rho[1] = pow( s[1]/(e[1]*Gamma_1), -1.0E0/Gamma_1  );
         }
         else // outflow
         {
             s[1] = s[2];
             e[1] = a[1]*a[1]/(Gamma*Gamma_1);
             rho[1] = pow( s[1]/(e[1]*Gamma_1), -1.0E0/Gamma_1  );
             u[1] = vn[1]*nx + (u[2] - vn[2]*nx);
             v[1] = vn[1]*ny + (v[2] - vn[2]*ny);
         }
   }// end if super/sub sonic    
   
   return;
}


/////////////////////////////////////////////////////////////////////
/*
C-type grid
[0] bctype = ICGRID
[1] dir, 0 i / 1 j
[2] in_out 1 / DIM
[3],[4]  start index, end index of  B (I dir)  // start index, end index of  A (J dir)  
[5],[6]  start index, end index of  A (I dir)  // start index, end index of  B (J dir)
[8]  NA
*/
int bc_cgrid(const int ibc )
{
#ifdef DEBUG
    std::cout<<"ICGRID BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    if(bcdat[ibc][1]==0)  // i dir
    {
        int Ajstrt,Ajend, Bjstrt,Bjend,jl,   in, s1, s2, t1, t2;
        if(bcdat[ibc][2]==1) // i = 1
        {
            s1 = 2;
            s2 = 3;
            t1 = 1;
            t2 = 0;
        }
        else if(bcdat[ibc][2]==IDIM) // i = IDIM
        {
            s1 = IDIM;
            s2 = IDIM-1;
            t1 = IDIM+1;
            t2 = IDIM+2;
        }
        else
        {
            std::cerr<<"Error C-grid BC data, in_out != 1 or IDIM"<<std::endl;
            return 0;
        }
            
        Ajstrt = bcdat[ibc][3];
        Ajend =  bcdat[ibc][4];
        jl = Ajend - Ajstrt + 1;
        Bjstrt = bcdat[ibc][5];
        Bjend  = bcdat[ibc][6];

        if(Ajend<Ajstrt)
        {
            std::cerr<<"Error C-grid BC data, Ajend<Ajstrt"<<std::endl;
            return 0;
        }
        if(Bjend>Bjstrt)
        {
            std::cerr<<"Error C-grid BC data, Bjend>Bjstrt"<<std::endl;
            return 0;
        }

        if( abs(Ajstrt-Ajend)!=abs(Bjstrt-Bjend) 
            || Ajstrt<2 || Bjstrt>JDIM
            || Ajend>Bjend)
        {
            std::cerr<<"Error C-grid BC range"<<std::endl;
            return 0;
        }

        for(int j=0; j<jl;j++)     
        {
            
            
            for(int dir=0;dir<DIR;dir++)
            {
                Q[t1][Ajstrt+j][dir] = Q[s1][Bjstrt-j][dir];
                Q[t2][Ajstrt+j][dir] = Q[s2][Bjstrt-j][dir];
            
                Q[t1][Bjstrt-j][dir] = Q[s1][Ajstrt+j][dir];
                Q[t2][Bjstrt-j][dir] = Q[s2][Ajstrt+j][dir];
            }
            macro(Q[t1][Ajstrt+j],M[t1][Ajstrt+j][0],M[t1][Ajstrt+j][1],M[t1][Ajstrt+j][2],M[t1][Ajstrt+j][3]);
            macro(Q[t2][Ajstrt+j],M[t2][Ajstrt+j][0],M[t2][Ajstrt+j][1],M[t2][Ajstrt+j][2],M[t2][Ajstrt+j][3]);
            macro(Q[t1][Bjstrt-j],M[t1][Bjstrt-j][0],M[t1][Bjstrt-j][1],M[t1][Bjstrt-j][2],M[t1][Bjstrt-j][3]);
            macro(Q[t2][Bjstrt-j],M[t2][Bjstrt-j][0],M[t2][Bjstrt-j][1],M[t2][Bjstrt-j][2],M[t2][Bjstrt-j][3]);
        }
        
    }
    else if(bcdat[ibc][1]==1) // jdir
    {
        int Aistrt,Aiend,  Bistrt,Biend,il,  in, s1, s2, t1, t2;
        if(bcdat[ibc][2]==1) 
        {
            s1 = 2;
            s2 = 3;
            t1 = 1;
            t2 = 0;
        }
        else if(bcdat[ibc][2]==JDIM) 
        {
            s1 = JDIM;
            s2 = JDIM-1;
            t1 = JDIM+1;
            t2 = JDIM+2;
        }
        else
        {
            std::cerr<<"Error tangent BC data, in_out != 1 or IDIM"<<std::endl;
            return 0;
        }
        
        Aistrt = bcdat[ibc][3];
        Aiend =  bcdat[ibc][4];
        il = Aiend - Aistrt + 1;
        Bistrt = bcdat[ibc][5];
        Biend  = bcdat[ibc][6];

        if(Aiend<Aistrt)
        {
            std::cerr<<"Error C-grid BC data, Ajend<Ajstrt"<<std::endl;
            return 0;
        }
        if(Biend>Bistrt)
        {
            std::cerr<<"Error C-grid BC data, Bjend>Bjstrt"<<std::endl;
            return 0;
        }        
        if( abs(Aistrt-Aiend)!=abs(Bistrt-Biend) 
            || Aistrt<2 || Bistrt>IDIM
            || Aiend>Biend)
        {
            std::cerr<<"Error C-grid BC range"<<std::endl;
            return 0;
        }
        
        
        for(int i=0; i<il;i++)
        {
            // setting macro states on the boudanry is just for output
            // no compact on computation since we only need Q
            //std::cout<<Aistrt+i<<"  "<<Bistrt-i<<std::endl;
            
            for(int dir=0;dir<DIR;dir++)
            {            
               Q[Aistrt+i][t1][dir] = Q[Bistrt-i][s1][dir];
               Q[Aistrt+i][t2][dir] = Q[Bistrt-i][s2][dir];            
               Q[Bistrt-i][t1][dir] = Q[Aistrt+i][s1][dir];
               Q[Bistrt-i][t2][dir] = Q[Aistrt+i][s2][dir];
            }
            macro(Q[Aistrt+i][t1],M[Aistrt+i][t1][0],M[Aistrt+i][t1][1],M[Aistrt+i][t1][2],M[Aistrt+i][t1][3]);
            macro(Q[Aistrt+i][t2],M[Aistrt+i][t2][0],M[Aistrt+i][t2][1],M[Aistrt+i][t2][2],M[Aistrt+i][t2][3]);
            macro(Q[Bistrt-i][t1],M[Bistrt-i][t1][0],M[Bistrt-i][t1][1],M[Bistrt-i][t1][2],M[Bistrt-i][t1][3]);
            macro(Q[Bistrt-i][t2],M[Bistrt-i][t2][0],M[Bistrt-i][t2][1],M[Bistrt-i][t2][2],M[Bistrt-i][t2][3]);
        }
    }
    
    return 1;
}





/////////////////////////////////////////////////////////////////////
// Tangent boundary
int bc_tan(const int ibc)
{
#ifdef DEBUG    
    std::cout<<"ITAN BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    
    double rr,rU,re,rp;    
    get_ref(rr,rU,re,rp);    
    
    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1, in2;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        inout = bcdat[ibc][2];
        if(jend<jstrt || jend>JDIM || jstrt<2)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }        
        
        for(int j=jstrt;j<=jend;j++)
        {
            double nx,ny;
            // get inward normal vector
            if(inout==1)   // is inward
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                in2 = 3;
                itan = 1;  // wall "cell"
                nx = S[inout][j][0][0];
                ny = S[inout][j][0][1];
            }
            else if(inout==IDIM)// is outward
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                in2 = IDIM-1;
                itan = IDIM+1;  // wall "cell"
                nx = -S[inout][j][0][0];
                ny = -S[inout][j][0][1];
            }
            else
            {
                std::cerr<<"Error tangent BC data, in_out != 1 or IDIM"<<std::endl;
                return 0;
            }
            
            if(bcdat[ibc][0]==ITAN)
            {
                mirror( nx, ny, Q[in1][j], Q[itan][j]);
                //  on the wall
                for(int dir=0;dir<DIR;dir++)
                      Q[itan][j][dir] = 0.5E0*( Q[itan][j][dir] + Q[in1][j][dir] );
                macro(Q[itan][j], M[itan][j][0], M[itan][j][1], M[itan][j][2],
                                  M[itan][j][3]);
            }
            else if(bcdat[ibc][0]==ISYM)
            {
                mirror( nx, ny, Q[in1][j], Q[itan][j]);
                mirror( nx, ny, Q[in2][j], Q[  gh][j]);

                macro(Q[itan][j], M[itan][j][0], M[itan][j][1], M[itan][j][2],
                                  M[itan][j][3]);
                macro(Q[gh][j], M[gh][j][0], M[gh][j][1], M[gh][j][2],
                                  M[gh][j][3]);
            }
            else
            {
                std::cerr<<"Error! bc_tan can only be used for slip wall and "
                        <<"symetric boundaries."<<std::endl;
                return 0;
            }
        }// end for(j)
    }
    else if(bcdat[ibc][1]==1)   // J direction
        {
        int istrt,iend,inout,itan,gh,in1, in2;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        inout = bcdat[ibc][2];
        if(iend<istrt || iend>IDIM || istrt<2)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
        
        for(int i=istrt;i<=iend;i++)
        {
            double nx,ny;
            //std::cout<<i<<std::endl;
            // get inward normal vector
            if(inout==1)   // is inward
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                in2 = 3;
                itan = 1;  // wall "cell"
                nx = S[i][inout][1][0];
                ny = S[i][inout][1][1];
            }
            else if(inout==JDIM)  // is outward
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                in2 = JDIM-1;
                itan = JDIM+1;  // wall "cell"
                nx = -S[i][inout][1][0];
                ny = -S[i][inout][1][1];
            }
            else
            {
                std::cerr<<"Error tangent BC data, in_out != 1 or JDIM"<<std::endl;
                return 0;
            }

            if(bcdat[ibc][0]==ITAN)
            {
                mirror( nx, ny, Q[i][in1], Q[i][itan]);
                for(int dir=0;dir<DIR;dir++)
                     Q[i][itan][dir] = 0.5E0*( Q[i][itan][dir] + Q[i][in1][dir] );
                macro(Q[i][itan], M[i][itan][0], M[i][itan][1], M[i][itan][2], 
                                  M[i][itan][3]);
            }
            else if(bcdat[ibc][0]==ISYM)
            {
               mirror( nx, ny, Q[i][in1], Q[i][itan]);
               mirror( nx, ny, Q[i][in2], Q[i][gh]);

               macro(Q[i][itan], M[i][itan][0], M[i][itan][1], M[i][itan][2],
                                  M[i][itan][3]);
               macro(Q[i][gh], M[i][gh][0], M[i][gh][1], M[i][gh][2],
                                  M[i][gh][3]);
            }
            else
            {
               std::cerr<<"Error! bc_tan can only be used for slip wall and "
                       <<"symetric boundaries."<<std::endl;
                return 0;
            }
           
        }// end for(j)
    }
    else
    {
        std::cout<<"Invalid direction of ITAN boundary "<<std::endl;
        return -1;
    }
    
    return 1;
}




/////////////////////////////////////////////////////////////////////
// Tangent boundary 2
// compute macro vars with tranditional method, and then get feq
int bc_tan2(const int ibc)
{
#ifdef DEBUG    
    std::cout<<"ITAN BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    
    double rr,rU,re,rp;    
    get_ref(rr,rU,re,rp);    
    
    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1,in2;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        inout = bcdat[ibc][2];
        if(jend<jstrt || jend>JDIM || jstrt<2)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }        
        
        for(int j=jstrt;j<=jend;j++)
        {
            double nx,ny;
            // get inward normal vector
            if(inout==1)   // is inward
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                in2 = 3;
                itan = 1;  // wall "cell"
                nx = S[inout][j][0][0];
                ny = S[inout][j][0][1];
            }
            else if(inout==IDIM)// is outward
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                in2 = IDIM-1;
                itan = IDIM+1;  // wall "cell"
                nx = -S[inout][j][0][0];
                ny = -S[inout][j][0][1];
            }
            else
            {
                std::cerr<<"Error tangent BC data, in_out != 1 or IDIM"<<std::endl;
                return 0;
            }

//  --------  start of extrapolating wall            
            // compute pressure on the wall
            double p1 = get_pressure(M[in1][j][0],M[in1][j][3]);
            double p2 = get_pressure(M[in2][j][0],M[in2][j][3]);            
            double pw = (9.E0 * p1 - p2) / 8.0E0 ;  // 2nd order
            if (pw<0.0E0) pw = p1;
            double rho1,rhow, ew;
            rho1 = M[in1][j][0];
            //( PW / (PN/RHON**GAMMA) )**(1.E0/GAMMA)         
            rhow = pow( pw / ( p1/pow(rho1, Gamma) ) , 1.E0/Gamma );
            ew  = pw/(Gamma_1 * M[itan][j][0]);
            M[itan][j][0] = rhow;
            M[itan][j][3] = ew;
            
            double vn,vnx,vny;
            vn = nx*M[itan][j][1] + ny*M[itan][j][2];
            vnx = vn*nx;
            vny = vn*ny;
            M[itan][j][1] = M[in1][j][1] - vnx;
            M[itan][j][2] = M[in1][j][2] - vny;
            get_feq(M[itan][j][0], M[itan][j][1], M[itan][j][2], M[itan][j][3], Q[itan][j]);                                                    
           
        }// end for(j)
    }
    else if(bcdat[ibc][1]==1)   // J direction
        {
        int istrt,iend,inout,itan,gh,in1,in2;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        inout = bcdat[ibc][2];
        if(iend<istrt || iend>IDIM || istrt<2)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
        
        for(int i=istrt;i<=iend;i++)
        {
            double nx,ny;
            //std::cout<<i<<std::endl;
            // get inward normal vector
            if(inout==1)   // is inward
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                in2 = 3;
                itan = 1;  // wall "cell"
                nx = S[i][inout][1][0];
                ny = S[i][inout][1][1];
            }
            else if(inout==JDIM)  // is outward
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                in2 = JDIM-1;
                itan = JDIM+1;  // wall "cell"
                nx = -S[i][inout][1][0];
                ny = -S[i][inout][1][1];
            }
            else
            {
                std::cerr<<"Error tangent BC data, in_out != 1 or JDIM"<<std::endl;
                return 0;
            }
//  --------  start of extrapolating wall            
            // compute pressure on the wall
            double p1 = get_pressure(M[i][in1][0],M[i][in1][3]);
            double p2 = get_pressure(M[i][in2][0],M[i][in2][3]);            
            double pw = (9.E0 * p1 - p2) / 8.0E0 ;  // 2nd order
            if (pw<0.0E0) pw = p1;
            double rho1,rhow, ew;
            rho1 = M[i][in1][0];       
            rhow = pow( pw / ( p1/pow(rho1, Gamma) ) , 1.E0/Gamma );
            ew  = pw/(Gamma_1 * rhow);
            
            M[i][itan][0] = rhow;
            M[i][itan][3] = ew;
            
            double vn,vnx,vny;
            vn = nx*M[i][itan][1] + ny*M[i][itan][2];
            vnx = vn*nx;
            vny = vn*ny;
            M[i][itan][1] = M[i][in1][1] - vnx;
            M[i][itan][2] = M[i][in1][2] - vny;
            get_feq(M[i][itan][0], M[i][itan][1], M[i][itan][2], M[i][itan][3], Q[i][itan]);
           
        }// end for(j)
    }
    else
    {
        std::cout<<"Invalid direction of ITAN boundary "<<std::endl;
        return -1;
    }
    
    return 1;
}

////////////////////////////////////////////////////////////////////////////////
// Supersonic outflow,
int bc_spout(const int ibc )
{
#ifdef DEBUG
    std::cout<<"ISPOUT BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif

    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }
        inout = bcdat[ibc][2];
        
        
        for(int j=jstrt;j<=jend;j++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                itan = IDIM+1;  // bd "cell"
            }
            
            M[itan][j][0] = M[in1][j][0];
            M[itan][j][1] = M[in1][j][1];
            M[itan][j][2] = M[in1][j][2];
            M[itan][j][3] = M[in1][j][3];
            
            M[gh][j][0] = M[in1][j][0];  
            M[gh][j][1] = M[in1][j][1];  
            M[gh][j][2] = M[in1][j][2];  
            M[gh][j][3] = M[in1][j][3];  
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[itan][j][dir] =  Q[in1][j][dir];
                   Q[gh  ][j][dir] =  Q[in1][j][dir];
            }
            
        }//end for(j)
        
    }
    else if(bcdat[ibc][1]==1)   // J direction
    {
        int istrt,iend,inout,itan,gh,in1;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
        inout = bcdat[ibc][2];
        
        
        for(int i=istrt;i<=iend;i++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                itan = JDIM+1;  // bd "cell"
            }
            
            M[i][itan][0] = M[i][in1][0];
            M[i][itan][1] = M[i][in1][1];
            M[i][itan][2] = M[i][in1][2];
            M[i][itan][3] = M[i][in1][3];
                                         
            M[i][gh][0] = M[i][in1][0];  
            M[i][gh][1] = M[i][in1][1];  
            M[i][gh][2] = M[i][in1][2];  
            M[i][gh][3] = M[i][in1][3];              
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][itan][dir] =  Q[i][in1][dir];
                   Q[i][gh  ][dir] =  Q[i][in1][dir];
            }
            
        }//end for(i)
        
    }
    else
    {
        std::cout<<"Invalid in/out of ISPOUT boundary "<<std::endl;
        return -1;
    }    
    
    return 1;
}


// Supersonic outflow,
int bc_spin(const int ibc )
{
#ifdef DEBUG
    std::cout<<"IFSPIN BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif

    double cr,cu,cv,ce;
    double cfeq[DIR];
    cr = bcm[ibc][0];
    cu = bcm[ibc][1];
    cv = bcm[ibc][2];
    ce = bcm[ibc][3];
    
    get_feq(cr,cu,cv,ce,cfeq);

    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }    
        inout = bcdat[ibc][2];
        
        for(int j=jstrt;j<=jend;j++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                itan = IDIM+1;  // bd "cell"
            }
            // bd
            M[itan][j][0] = cr;
            M[itan][j][1] = cu;
            M[itan][j][2] = cv;
            M[itan][j][3] = ce;
            // ghost
            M[gh][j][0] = cr;
            M[gh][j][1] = cu;
            M[gh][j][2] = cv;
            M[gh][j][3] = ce;
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[itan][j][dir] =  cfeq[dir];
                   Q[gh  ][j][dir] =  cfeq[dir];
            }
            
        }//end for(j)
        
    }
    else if(bcdat[ibc][1]==1)   // J direction
    {
        int istrt,iend,inout,itan,gh,in1;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
        inout = bcdat[ibc][2];
                
        for(int i=istrt;i<=iend;i++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                itan = JDIM+1;  // bd "cell"
            }
            
            M[i][itan][0] = cr;
            M[i][itan][1] = cu;
            M[i][itan][2] = cv;
            M[i][itan][3] = ce;            
            
            M[i][gh][0] = cr;
            M[i][gh][1] = cu;
            M[i][gh][2] = cv;
            M[i][gh][3] = ce;
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][itan][dir] =  cfeq[dir];
                   Q[i][gh  ][dir] =  cfeq[dir];
            }
            
        }//end for(i)
        
    }
    else
    {
        std::cout<<"Invalid in/out of ISPIN boundary "<<std::endl;
        return -1;
    }    
    
    return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Subsonic outflow,
int bc_suout(const int ibc )
{
#ifdef DEBUG
    std::cout<<"ISUOUT BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    double cp;
    cp = bcm[ibc][3];   // Pressure
    double t1,t2,t3,t4;
    get_ref(t1,t2,t3,t4);
    //std::cout<<myrank<<" CPU, subsonic outflow with pressure = "<<cp*t1*t3
    //	               <<" internal energy ="<<std::endl;

    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error subsonic outflow BC, jend<jtrt"<<std::endl;
            return 0;
        }
        inout = bcdat[ibc][2];
        
        
        for(int j=jstrt;j<=jend;j++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = IDIM+2;    // ghost cell
                in1 = IDIM;     // 1st internal cell
                itan = IDIM+1;  // bd "cell"
            }
            
            M[itan][j][0] = M[in1][j][0];
            M[itan][j][1] = M[in1][j][1];
            M[itan][j][2] = M[in1][j][2];
            M[itan][j][3] = cp/(Gamma_1*M[itan][j][0]);
            
            M[gh][j][0] = M[itan][j][0];  
            M[gh][j][1] = M[itan][j][1];  
            M[gh][j][2] = M[itan][j][2];  
            M[gh][j][3] = M[itan][j][3];  
            
            double cfeq[DIR];
            get_feq(M[itan][j][0], M[itan][j][1],M[itan][j][2],M[itan][j][3],cfeq);
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[itan][j][dir] =  cfeq[dir];
                   Q[gh  ][j][dir] =  cfeq[dir];
            }
            
        }//end for(j)
        
    }
    else if(bcdat[ibc][1]==1)   // J direction
    {
        int istrt,iend,inout,itan,gh,in1;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
        inout = bcdat[ibc][2];
        
        
        for(int i=istrt;i<=iend;i++)
        {
            if(inout==1)
            {
                gh = 0;    // ghost cell
                in1 = 2;   // 1st internal cell
                itan = 1;  // bd "cell"
            }
            else  
            {
                gh = JDIM+2;    // ghost cell
                in1 = JDIM;     // 1st internal cell
                itan = JDIM+1;  // bd "cell"
            }
            
            M[i][itan][0] = M[i][in1][0];
            M[i][itan][1] = M[i][in1][1];
            M[i][itan][2] = M[i][in1][2];
            M[i][itan][3] = cp/(Gamma_1*M[i][itan][0]);
                                         
            M[i][gh][0] = M[i][itan][0];  
            M[i][gh][1] = M[i][itan][1];  
            M[i][gh][2] = M[i][itan][2];  
            M[i][gh][3] = M[i][itan][3];              
            
            double cfeq[DIR];
            get_feq(M[i][itan][0], M[i][itan][1],M[i][itan][2],M[i][itan][3],cfeq);
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][itan][dir] =  cfeq[dir];
                   Q[i][gh  ][dir] =  cfeq[dir];
            }
            
        }//end for(i)
        
    }
    else
    {
        std::cout<<"Invalid in/out of ISUOUT boundary "<<std::endl;
        return -1;
    }    
    
    return 1;
}



/////////////////////////////////////////////////////////////////
/*
bcdat[ibc][0]: BC type
          [1]: BC direction, 0=i, 1=j
          [2]: NA (1 and DIM)
          [3]: istrt
          [4]: iend
          [5]: jstrt
          [6]: jend
          [7]: profile index
*/
int bc_perd(const int ibc)
{
#ifdef DEBUG
    std::cout<<"Periodic BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif

    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        if(jend<jstrt || jend>JDIM+1 || jstrt<1)
        {
            std::cerr<<"Error tangent BC, jend<jtrt"<<std::endl;
            return 0;
        }

        for(int j=jstrt;j<=jend;j++)
        {
            // just for output
            for(int iv=0;iv<NVAR;iv++)
            {
                M[1     ][j][iv] = M[IDIM][j][iv];
                M[IDIM+1][j][iv] = M[1   ][j][iv];
            }
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[0][j][dir] =  Q[IDIM-1][j][dir];
                   Q[1][j][dir] =  Q[IDIM  ][j][dir];
                   Q[IDIM+1][j][dir] = Q[2][j][dir];
                   Q[IDIM+2][j][dir] = Q[3][j][dir];
            }
            
        }//end for(j)
        
    }
    else if(bcdat[ibc][1]==1)   // J direction
    {
        int istrt,iend;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        if(iend<istrt || iend>IDIM+1 || istrt<1)
        {
            std::cerr<<"Error tangent BC, iend<itrt"<<std::endl;
            return 0;
        }
                
        
        for(int i=istrt;i<=iend;i++)
        {
            for(int iv=0;iv<NVAR;iv++)
            {
                M[i][1     ][iv] = M[i][JDIM][iv];
                M[i][JDIM+1][iv] = M[i][1     ][iv];
            }
            
            for(int dir=0;dir<DIR;dir++)
            {
                   Q[i][0][dir] =  Q[i][JDIM-1][dir];
                   Q[i][1][dir] =  Q[i][JDIM  ][dir];
                   Q[i][JDIM+1][dir] = Q[i][2][dir];
                   Q[i][JDIM+2][dir] = Q[i][3][dir];
            }
            
        }//end for(i)
        
    }
    else
    {
        std::cout<<"Invalid in/out of ISPIN boundary "<<std::endl;
        return -1;
    }    
    
    return 1;
}


/////////////////////////////////////////////////////////////////////
// nonslip wall boundary
/*
bcdat[ibc][0]: BC type
          [1]: BC direction, 0=i, 1=j
          [2]: 1 or DIM
          [3]: istrt
          [4]: iend
          [5]: jstrt
          [6]: jend
          [7]: adiabadic=-1, const temperature=-2, const heat flux=-3
               prescribed temperture in a array = aux-array index
*/
int bc_wall(const int ibc)
{
#ifdef DEBUG    
    std::cout<<"IWALL BC ";
    for(int dbi=1;dbi<8;dbi++) std::cout<<" "<<bcdat[ibc][dbi];
    std::cout<<std::endl;
#endif
    double rr,rU,re,rp;    
    get_ref(rr,rU,re,rp); 
    rp = rr*re;
    
    double cr,cu,cv,ce;
    cr = bcm[ibc][0];   // NA
    cu = bcm[ibc][1];  
    cv = bcm[ibc][2];
    ce = bcm[ibc][3];
    //std::cout<<"state on wall "<<cr*rr<<" "<<cu*rU<<" "<<cv*rU<<" "<<ce*re<<std::endl;
    
    if(bcdat[ibc][1]==0)   // I direction
    {
        int jstrt,jend,inout,itan,gh,in1, in2, in3,wthermal;
        jstrt = bcdat[ibc][5];
        jend  = bcdat[ibc][6];
        inout = bcdat[ibc][2];
        wthermal = bcdat[ibc][7];
        if(jend<jstrt || jend>JDIM || jstrt<2)
        {
            std::cerr<<"Error nonslip BC, jend<jtrt"<<std::endl;
            return 0;
        }  
              
        // get inward normal vector
        if(inout==1)   // is inward
        {
            gh = 0;    // ghost cell
            in1 = 2;   // 1st internal cell
            in2 = 3;
            in3 = 4;
            itan = 1;  // wall "cell"
        }
        else if(inout==IDIM)// is outward
        {
            gh = IDIM+2;    // ghost cell
            in1 = IDIM;     // 1st internal cell
            in2 = IDIM-1;
            in3 = IDIM-2;
            itan = IDIM+1;  // wall "cell"
        }
        else
        {
            std::cerr<<"Error nonslip BC data, in_out != 1 or IDIM"<<std::endl;
            return 0;
        }    
                
        for(int j=jstrt;j<=jend;j++)
        {

//  --------  start of extrapolating wall            
            // compute pressure on the wall
            double p1 = get_pressure(M[in1][j][0],M[in1][j][3]);
            double p2 = get_pressure(M[in2][j][0],M[in2][j][3]);

//            double pw = p1;
            double pw = (9.E0 * p1 - p2) / 8.0E0 ;  // 2nd order
            if (pw<0.0E0) pw = p1;  // 1st order
            
            double ew;
            if(wthermal==-1)  // adiabadic wall
            {
                ew = (9.E0 * M[in1][j][3] - M[in2][j][3]) / 8.0E0 ; // 2nd order
                if(ew<0.0E0) ew = M[in1][j][3]; // 1st order 
            }
            else if(wthermal==-2)  // specified temp
                ew = ce;
            else if(wthermal==-3)  // const flux
                ew = (9.E0 * M[in1][j][3] - M[in2][j][3]) / 8.0E0  - 0.5E0*ce;
            else
            {
                std::cerr<<"Sorry, only adiabadic or constant temperature wall"
                         <<std::endl;
                return 0;
            }
            
            double uw,vw;
            uw=cu;
            vw=cv;
            double rhow = pw/(ew*Gamma_1); 
            M[itan][j][0] = rhow;       
            M[itan][j][1] = uw;
            M[itan][j][2] = vw;
            M[itan][j][3] = ew;           
            get_feq(rhow,uw,vw,ew, Q[itan][j]);
            
            double neq1[DIR],neq2[DIR];
            get_feq(M[in1][j][0],M[in1][j][1],M[in1][j][2],M[in1][j][3],
                    neq1);
            //get_feq(M[in2][j][0],M[in2][j][1],M[in2][j][2],M[in2][j][3],
            //        neq2);
            // extrapolating non-equlibrium part
            for(int idir=0;idir<DIR;idir++)
            {
                neq1[idir] = Q[in1][j][idir] - neq1[idir];
                //neq2[idir] = Q[in2][j][idir] - neq2[idir];
                Q[itan][j][idir] += neq1[idir]; //1.5E0*neq1[idir] -0.5E0*neq2[idir];
                Q[gh][j][idir]  = Q[itan][j][idir];
            }
//  --------  end of extrapolating wall
           
        }// end for(j)
    }
    
    
    else if(bcdat[ibc][1]==1)   // J direction
    {
        int istrt,iend,inout,itan,gh,in1, in2, in3,wthermal;
        istrt = bcdat[ibc][3];
        iend  = bcdat[ibc][4];
        inout = bcdat[ibc][2];
        wthermal = bcdat[ibc][7];
        if(iend<istrt || iend>IDIM || istrt<2)
        {
            std::cerr<<"Error nonslip BC, iend<istrt"<<std::endl;
            return 0;
        }        

            // get inward normal vector
        if(inout==1)   // is inward
        {
            gh = 0;    // ghost cell
            in1 = 2;   // 1st internal cell
            in2 = 3;
            in3 = 4;
            itan = 1;  // wall "cell"
        }
        else if(inout==JDIM)// is outward
        {
            gh = JDIM+2;    // ghost cell
            in1 = JDIM;     // 1st internal cell
            in2 = JDIM-1;
            in3 = JDIM-2;
            itan = JDIM+1;  // wall "cell"
        }
        else
        {
            std::cerr<<"Error nonslip BC data, in_out != 1 or JDIM"<<std::endl;
            return 0;
        }
                    
        for(int i=istrt;i<=iend;i++)
        {

//  --------  start of extrapolating wall            
            // compute pressure on the wall
            double p1 = get_pressure(M[i][in1][0],M[i][in1][3]);
            double p2 = get_pressure(M[i][in2][0],M[i][in2][3]);
            
//            double pw = p1;
          double pw = (9.E0 * p1 - p2) / 8.0E0 ;  // 2nd order
          if (pw<0.0E0) pw = p1;
            
            double ew;
            if(wthermal==-1)  // adiabadic wall
            {
                ew = (9.E0 * M[i][in1][3] - M[i][in2][3]) / 8.0E0 ;  // 2nd extrapolate
                if (ew<0.0E0) ew = M[i][in1][3]; 
            }
            else if(wthermal==-2)  // specified temp
                ew = ce;
            else if(wthermal==-3)  // const flux, here ce is the normal gradient of the internal energy, @e/@n
                ew = (9.E0 * M[i][in1][3] - M[i][in2][3]) / 8.0E0  - 0.5E0*ce;
            else
            {
                std::cerr<<"Sorry, only adiabadic or constant temperature wall"
                         <<std::endl;
                return 0;
            }
            
            double uw,vw;
            uw=cu;
            vw=cv;
            double rhow = pw/(ew*Gamma_1);   
            //std::cout<<"rho= "<<rhow*rr<<"   e= "<<ew*re<<"  p= "<<pw*rp<<std::endl;
            
            M[i][itan][0] = rhow;       
            M[i][itan][1] = uw;
            M[i][itan][2] = vw;
            M[i][itan][3] = ew;
            
            
            //std::cout<<"Wall "<<rhow<<" "<<pw<<" "<<ew<<std::endl;
            get_feq(rhow,uw,vw,ew, Q[i][itan]);
            
            double neq1[DIR],neq2[DIR];
            get_feq(M[i][in1][0],M[i][in1][1],M[i][in1][2],M[i][in1][3],
                    neq1);
            //get_feq(M[i][in2][0],M[i][in2][1],M[i][in2][2],M[i][in2][3],
            //        neq2);
            
            // extrapolating non-equlibrium part
            
            for(int idir=0;idir<DIR;idir++)
            {
                neq1[idir] = Q[i][in1][idir] - neq1[idir];
                //neq2[idir] = Q[i][in2][idir] - neq2[idir];
                Q[i][itan][idir] += neq1[idir]; //(1.5E0*neq1[idir] -0.5E0*neq2[idir]);
                Q[i][gh][idir]  = Q[i][itan][idir];
            }

//  --------  end of extrapolating wall
                                          
           
        }// end for(j)
    }
    else
    {
        std::cout<<"Invalid direction of IWALL boundary "<<std::endl;
        return -1;
    }
    
    return 1;
}


/////////////////////////////////////////////////////////////////////
/*
moving shock on the top BC
*/
int bc_usrdef(const int ibc, double time)
{

double rr,rU,re,rp;
get_ref(rr,rU,re,rp);

int bdir,bside,istrt, iend, jstrt, jend;

bdir  = bcdat[ibc][1];
bside = bcdat[ibc][2];
istrt = bcdat[ibc][3];
iend  = bcdat[ibc][4];
jstrt = bcdat[ibc][5];
jend  = bcdat[ibc][6];

int g1, g2;
if(bdir==0)  // I dir
  if(bside==1)
  {
    g1 = 0;
    g2 = 1;
  }
  else
  {
    g1 = IDIM+1;
    g2 = IDIM+2;
  }
else   // J dir
  if(bside==1)
  {
    g1 = 0;
    g2 = 1;
  }
  else
  {
    g1 = JDIM+1;
    g2 = JDIM+2;
  }

double cx,cy;
if(bdir==1)  // J dir
  for(int i=istrt;i<=iend;i++)
  for(int j=g1;j<=g2;j++)
  {
     cx = RC[i][2][0];
     if(bside==1) 
     {
         if(j==g2)
            cy = RC[i][2][1] - (RC[i][3][1] - RC[i][2][1]);
         else // j==g1
            cy = RC[i][2][1] - 2.0E0*(RC[i][3][1] - RC[i][2][1]);
     }
     else // bside==JDIM
     {
         if(j==g1)
            cy = RC[i][JDIM][1] + (RC[i][JDIM][1]-RC[i][JDIM-1][1]);
         else
            cy = RC[i][JDIM][1] + 2.0E0*(RC[i][JDIM][1]-RC[i][JDIM-1][1]);
     }

     get_anm(time, cx, cy, M[i][j]);
     M[i][j][0] /= rr;
     M[i][j][1] /= rU;
     M[i][j][2] /= rU;
     M[i][j][3] /= re;
     get_feq(M[i][j][0],M[i][j][1],M[i][j][2],M[i][j][3], Q[i][j]);
  }
else // I dir
  for(int i=g1;i<=g2;i++)
  for(int j=2;j<=JDIM;j++)
  {
     cy = RC[2][j][1];
     if(bside==1)
     {
        if(i==g2) // i=1
           cx = RC[2][j][0] - (RC[3][j][0] - RC[2][j][0]);
        else // i=0
           cx = RC[2][j][0] - 2.0E0*(RC[3][j][0] - RC[2][j][0]);
     }
     else // 
     {
        if(i==g1) // i =IDIM+1
          cx = RC[IDIM][j][0] + (RC[IDIM][j][0] - RC[IDIM-1][j][0]);
        else
          cx = RC[IDIM][j][0] + 2.0E0*(RC[IDIM][j][0] - RC[IDIM-1][j][0]);
     }

     get_anm(time, cx, cy, M[i][j]);
     M[i][j][0] /= rr;
     M[i][j][1] /= rU;
     M[i][j][2] /= rU;
     M[i][j][3] /= re;
     get_feq(M[i][j][0],M[i][j][1],M[i][j][2],M[i][j][3], Q[i][j]);
  }
    return 1;
}



void get_anm(const double ct, const double cx, const double cy, double* cruve)
{
double cr,cd,crho;
cr = sqrt(cx*cx+cy*cy);
cruve[0] = 1.0E0 + ct/cr;
cruve[1] = -cx/cr;
cruve[2] = -cy/cr;
cruve[3] = 1.0E-6/(Gamma_1*cruve[0]);

return;
}

