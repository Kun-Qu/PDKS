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

#include "grid_2d.h"

double R[IDIM+2][JDIM+2][2], VOL[IDIM+2][JDIM+2],RC[IDIM+2][JDIM+2][2];

double S[IDIM+2][JDIM+2][2][3];  // {XI_x,XI_y,S_XI}, {ETA_x,ETA_y,S_ETA}


int metric(void)
{
/*
C 1.  Calculate xsi(x), xsi(y). 
C     Calculate S[I][J][0][2] = |del(xsi)|.
C
*/
      for(int J = 1;J<= JDIM+1;J++)
      for(int I = 0;I<= IDIM+1;I++)
      {
       S[I][J][0][0] =   R[I][J][1] - R[I][J-1][1];                                
       S[I][J][0][1] = -(R[I][J][0] - R[I][J-1][0]);
                                    
       S[I][J][0][2] = sqrt( S[I][J][0][0]*S[I][J][0][0] +
                             S[I][J][0][1]*S[I][J][0][1] );
      }

      for(int J = 1;J<= JDIM+1;J++)
      for(int I = 0;I<= IDIM+1;I++)
      {
               if (S[I][J][0][2] > 0.0E0) 
               {
                  S[I][J][0][0] = S[I][J][0][0] / S[I][J][0][2];
                  S[I][J][0][1] = S[I][J][0][1] / S[I][J][0][2];
               }else
               {
                  S[I][J][0][0] = -1.0E0;
                  S[I][J][0][1] = -1.0E0;
                  S[I][J][0][2] = -1.0E0;
               }
      }
/*
C 2.  Calculate eta(x), eta(y) 
C     Calculate S[I][J][K][1][2] = |del(eta)|.
*/
      for(int I = 1;I<=IDIM+1;I++)
      for(int J = 0;J<=JDIM+1;J++)
      {
        S[I][J][1][0] =   R[I-1][J][1] - R[I][J][1];
        S[I][J][1][1] = -(R[I-1][J][0] - R[I][J][0]);
                                     
        S[I][J][1][2] = sqrt( S[I][J][1][0]*S[I][J][1][0] + 
                              S[I][J][1][1]*S[I][J][1][1] );
      }

      for(int I = 1;I<=IDIM+1;I++)
      for(int J = 0;J<=JDIM+1;J++)
      {
               if (S[I][J][1][2] > 0.0E0)
               {
                  S[I][J][1][0] = S[I][J][1][0] / S[I][J][1][2];
                  S[I][J][1][1] = S[I][J][1][1] / S[I][J][1][2];
               }else
               {
                  S[I][J][1][0] = -1.0E0;
                  S[I][J][1][1] = -1.0E0;
                  S[I][J][1][2] = -1.0E0;
               }
               
      }

      
/*
C 4.  Calculate the cell volume
*/

      for(int I = 1;I<= IDIM+1;I++)
      for(int J = 1;J<= JDIM+1;J++)
      {
          
      VOL[I][J] = 0.5E0 * ( (R[I  ][J][0] - R[I-1][J-1][0])*
                            (R[I-1][J][1] - R[I  ][J-1][1]) -
                            (R[I  ][J][1] - R[I-1][J-1][1])*
                            (R[I-1][J][0] - R[I  ][J-1][0]) );
      //if(I==2) std::cout<<"  "<<VOL[I][J]
      //                  <<"   "<<S[I][J][0][2]
      //                  <<"   "<<S[I][J][1][2]<<std::endl;
      }


/*
C 5.  Check for negative or zero areas
C
C I faces
*/
      int ISNEG  = 0;
      for(int J = 2;J<= JDIM;J++)
      for(int I = 0;I<= IDIM+1;I++)
               if (S[I][J][0][2] < 0.0E0)
               {
                  ISNEG++;
                  std::cout<<I<<" "<<J<<" neg I face"<<std::endl;
               }
                  
      if (ISNEG > 0)
         std::cout<<ISNEG<<" negtive I-faces"<<std::endl;
/*
C J faces
*/
      ISNEG  = 0;
      for(int I = 2;I<= IDIM;I++)
      for(int J = 0;J<= JDIM+1;J++)
               if (S[I][J][1][2] < 0.0E0)
                  ISNEG++;
                  
      if (ISNEG > 0)
         std::cout<<ISNEG<<" negtive J-faces"<<std::endl;
/*
C 6. Check for negative or zero volumes - left-handed coordinate
C    systems or collapsed cells
*/
      int IERROR = 0;
      for(int I = 2;I<= IDIM;I++)
      for(int J = 2;J<= JDIM;J++)
               if (VOL[I][J] <= 0.0E0)
                  IERROR ++;
                  
      if (IERROR > 0)
         std::cout<<IERROR<<" negtive vol"<<std::endl;
         
}


void extend_grid(void)
{
    // extend on the faces of j=1,JDIM
    for(int i=1;i<=IDIM;i++)
    for(int cd=0;cd<2;cd++)
    {
        R[i][0][cd] = R[i][1][cd];
        R[i][JDIM+1][cd] = R[i][JDIM][cd];
    }
    
    // extend on on the faces of i=1,IDIM
    for(int j=1;j<=JDIM;j++)
    for(int cd=0;cd<2;cd++)
    {
        R[0][j][cd] = R[1][j][cd];
        R[IDIM+1][j][cd] = R[IDIM][j][cd];
    }
    
    R[0][0][0] = R[1][1][0];
    R[0][0][1] = R[1][1][1];
    
    R[0][JDIM+1][0] = R[1][JDIM][0];
    R[0][JDIM+1][1] = R[1][JDIM][1];
    
    R[IDIM+1][JDIM+1][0] = R[IDIM][JDIM][0];
    R[IDIM+1][JDIM+1][1] = R[IDIM][JDIM][1];
    
    R[IDIM+1][0][0] = R[IDIM][1][0];
    R[IDIM+1][0][1] = R[IDIM][1][1];
    
    return;
}


void user_create_grid(void)
{
/*   
    const double RX=3.0E0, RY=6.0E0, PI =atan(1.0E0)*4.0E0,  tht = PI*5.0E0/12.0E0;
    double x,y,xi,et,dxi,det;
    dxi = 1.0E0/(IDIM-1);
    det = 1.0E0/(JDIM-1);
    
    for(int i=1;i<=IDIM;i++)
    for(int j=1;j<=JDIM;j++)
    {
        xi = (i-1.0E0)*dxi;
        et = (j-1.0E0)*det;
        
        x =  (RX - (RX-1.0E0)*xi)*cos(tht*(1.0E0 - et - et) + PI);
        y =  (RY - (RY-1.0E0)*xi)*sin(tht*(1.0E0 - et - et) + PI);
        
        R[i][j][0] = x;
        R[i][j][1] = y;
    }
*/


    const double dth=atan(1.0E0)*3.0E0/(IDIM-1), dr = 3.0E0/(JDIM-1),PI=atan(1.0E0)*4.0E0;
    double th,r,x,y;
    for(int i=1;i<=IDIM;i++)
    for(int j=1;j<=JDIM;j++)
    {
        th = -(i-1)*dth + PI*11.0E0/8.0E0;
        r =  (j-1)*dr + 1.0E0;
        
        x = r*cos(th);
        y = r*sin(th);
        
        R[i][j][0] = x;
        R[i][j][1] = y;
    }


/*
    double dxi,det,x,y,PI;
    PI = atan(1.0E0)*4.0E0;
    dxi = 1.0E0/(IDIM-1);
    det = 1.0E0/(JDIM-1);
    for(int i=1;i<=IDIM;i++)
    for(int j=1;j<=JDIM;j++)
    {
        x =(tan(((i-1.0E0)*dxi-0.5E0)*PI*0.5E0)+1.0E0)*2.0E0;
        y = (j-1.0E0)*det;
        R[i][j][0] = x;
        R[i][j][1] = y;
    }
*/
    return;
}

void cell_center(void)
{
// corners
    // i=1,j=1
    RC[1][1][0] = R[1][1][0];
    RC[1][1][1] = R[1][1][1];
    
    // i=1,j=JDIM+1
    RC[1][JDIM+1][0] = R[1][JDIM][0];
    RC[1][JDIM+1][1] = R[1][JDIM][1];
    
    // i=IDIM+1,j=1
    RC[IDIM+1][1][0] = R[IDIM][1][0];
    RC[IDIM+1][1][1] = R[IDIM][1][1];
    
    // i=IDIM+1,j=JDIM+1 
    RC[IDIM+1][JDIM+1][0] = R[IDIM][JDIM][0];
    RC[IDIM+1][JDIM+1][1] = R[IDIM][JDIM][1];
    
    
// The interior of the domain
    for(int i=2;i<=IDIM;i++)
    for(int j=2;j<=JDIM;j++)
    {
        RC[i][j][0] = 0.25E0* (R[i  ][j][0] + R[i-1][j-1][0] + 
                               R[i-1][j][0] + R[i  ][j-1][0]);
        RC[i][j][1] = 0.25E0* (R[i  ][j][1] + R[i-1][j-1][1] + 
                               R[i-1][j][1] + R[i  ][j-1][1]);
    }
    
// edges
    // i=1 & IDIM+1,  j=2..JDIM
    for(int j=2;j<=JDIM;j++)
    {   // i=1
        RC[1][j][0] = 0.5E0* (R[1][j-1][0] + R[1][j][0]); 
        RC[1][j][1] = 0.5E0* (R[1][j-1][1] + R[1][j][1]);
        // i=IDIM+1
        RC[IDIM+1][j][0] = 0.5E0* (R[IDIM+1][j-1][0] + R[IDIM+1][j][0]);
        RC[IDIM+1][j][1] = 0.5E0* (R[IDIM+1][j-1][1] + R[IDIM+1][j][1]);
    }

    // i=1..IDIM,  j=1 & JDIM+1
    for(int i=2;i<=IDIM;i++)
    {   // j=1
        RC[i][1][0] = 0.5E0* (R[i-1][1][0] + R[i][1][0]);
        RC[i][1][1] = 0.5E0* (R[i-1][1][1] + R[i][1][1]);
        // j=JDIM+1
        RC[i][JDIM+1][0] = 0.5E0* (R[i-1][JDIM][0] + R[i][JDIM][0]);
        RC[i][JDIM+1][1] = 0.5E0* (R[i-1][JDIM][1] + R[i][JDIM][1]);
    }



    return;
}

/*
A single-block 2D grid file of PLOT3D format reading function
Please make sure the grid file is right, for this function doesn't
do all format checking.
The format is
------------------------------
NX NY
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
------------------------------
*/
int read_grid(const char* plot3d_name)
{
    std::ifstream plot3d(plot3d_name);
    if(!plot3d)
    {
        std::cerr<<"Read_grid: can not open the plot3d file "
                 <<plot3d_name<<std::endl;
        return 0;
    }
    
    int  nx, ny;
    plot3d>>nx>>ny;
    if(nx!=IDIM)
    {
        std::cerr<<"Read_grid: NX of the grid file isn't equal to IDIM. Please reset "
         <<"IDIM in dks.h and rebuild the code again."<<std::endl;
        return 0;
    }
    else if(ny!=JDIM)
    {
        std::cerr<<"Read_grid: NY of the grid file isn't equal to JDIM. Please reset "
        <<"JDIM in dks.h and rebuild the code again."<<std::endl;
        return 0;
    }
    
    // read X
    for(int j=1;j<=JDIM;j++)
    for(int i=1;i<=IDIM;i++)
        plot3d>>R[i][j][0];
    
    // read Y
    for(int j=1;j<=JDIM;j++)
    for(int i=1;i<=IDIM;i++)
        plot3d>>R[i][j][1];
    
    return 1;
}
