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
#include "parallel.h"
#include "dks.h"
#include "mpi.h"


int myrank, cpusize, range[MAX_CPU][4];

int SYNCH=1, GATHR=1, ALL_UPDATED=1;

double ghost_left[(MAXDIM+MAXDIM)*DIR], 
       ghost_right[(MAXDIM+MAXDIM)*DIR],
       donor_left[(MAXDIM+MAXDIM)*DIR], 
       donor_right[(MAXDIM+MAXDIM)*DIR];
/*
A simple partion function, just partion along I index
*/ 
int partition(void)
{
	if( cpusize>MAX_CPU ) return 0;
		
  int istride = (IDIM-1)/cpusize;
  for(int rank=0;rank<cpusize;rank++)  // partition along I
  {
	range[rank][0] = rank==0 ? 2 : 1+istride*rank;
	range[rank][1] = 2;
	range[rank][2] = rank==cpusize-1 ? IDIM : (rank+1)*istride;
	range[rank][3] = JDIM;
	//std::cout<<myrank<<" ---- CPU "<<rank<<" I range is from "<<range[rank][0]<<" to "<<range[rank][2]<<std::endl;
  }
  
  return 1;
}


/* ////////////////////////////////////////////////////////////////////////////
Gather all elements of Q array on different partitions to CPU 0
*/

int gather_array(void)
{
   if(get_GATHR()) return 0;   
   double buff[(MAXDIM+3)*DIR];
   
   for(int rank=1;rank<cpusize;rank++)
   {
   	 int iis, iit, jjs, jjt;   
     my_range(rank, iis, iit, jjs, jjt);
     int data_len = (jjt-jjs+1)*DIR;
     
     for(int i=iis;i<=iit;i++)
     {     	   
       	 if(myrank == rank)
       	 {
       	   int index=0;
       	   for(int j=jjs;j<=jjt;j++)
     	     for(int dir=0;dir<DIR;dir++)
     	     {
     	   	    buff[index] = Q[i][j][dir];
     	   	    index ++;
     	     }     	   
     	   
     	     MPI_Send(buff, data_len, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
     	   }
     	   	 
     	   
         if(myrank == 0)
         {
           	MPI_Status mpi_err;
            MPI_Recv(buff, data_len, MPI_DOUBLE, rank, i, MPI_COMM_WORLD, &mpi_err);
            //std::cout<<"Recieved data["<<i<<"] from CPU "<<rank<<std::endl;
            int index = 0;
            for(int j=jjs;j<=jjt;j++)
            for(int dir=0;dir<DIR;dir++)
       	    {
       	     	  Q[i][j][dir] = buff[index];
       	   	    index ++;
       	    }
   //    	    std::cout<<"Recieved Q["<<i<<"] from CPU["<<rank<<"]"<<std::endl;
         }
         
         MPI_Barrier(MPI_COMM_WORLD);
     	   
     }// loop i
   	
   }// loop rank
   
   set_GATHR();
   return 1;   
}

/* //////////////////////////////////////////////////////
synchronize all ghost zones of all CPU
Here, the width of a ghost zone is assumed as 2
*/
int synchronize(int i_prd)
{
	  if(get_SYNCH()) return 1;
	  	
	   MPI_Barrier(MPI_COMM_WORLD);
	  
	   
        int box_left[4], box_right[4];
        box_left[0] = range[myrank][0];
        box_left[1] = 2;
        box_left[2] = range[myrank][0]+1;
        box_left[3] = JDIM;

        box_right[0] = range[myrank][2]-1;
        box_right[1] = 2;
        box_right[2] = range[myrank][2];
        box_right[3] = JDIM;
        
        int data_len = (box_left[2]-box_left[0]+1)*(box_left[3]-box_left[1]+1)*DIR;
        
        // set up donor_left
        int index = 0;
        for(int i=box_left[0];i<=box_left[2];i++)
        for(int j=box_left[1];j<=box_left[3];j++)
        for(int dir=0;dir<DIR;dir++)
        {
           donor_left[index] = Q[i][j][dir];
           index++;
        }

        // set up donor_right
        index = 0;
        for(int i=box_right[0];i<=box_right[2];i++)
        for(int j=box_right[1];j<=box_right[3];j++)
        for(int dir=0;dir<DIR;dir++)
        {
           donor_right[index] = Q[i][j][dir];
           index++;
        }
        
        // transfer data from my donor to neighbor's ghost
        MPI_Status mpi_err;
        int left=myrank-1,right=myrank+1;
        if(myrank==0)
        	  if(!i_prd)        	  	
               left = MPI_PROC_NULL;
            else
            	 left = cpusize-1;
        else if(myrank==cpusize-1)
        	  if(!i_prd) 
        	     right = MPI_PROC_NULL;
        	  else
        	  	 right = 0;
/*
 int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Status *status) 
*/            
        // send donnor to left and receive ghost from right
        MPI_Sendrecv(donor_left,  data_len, MPI_DOUBLE, left, 1, 
                     ghost_right, data_len, MPI_DOUBLE, right, 1, 
                     MPI_COMM_WORLD, &mpi_err);
        // send donnor to right and recieve ghost from left
        MPI_Sendrecv(donor_right,  data_len, MPI_DOUBLE, right, 2, 
                     ghost_left,   data_len, MPI_DOUBLE, left, 2, 
                     MPI_COMM_WORLD, &mpi_err);
                     
        // copy from ghost into Q
        // left ghost
        //if(myran>0)
        //{
           index = 0;
           for(int i=box_left[0]-2;i<=box_left[0]-1;i++)
           for(int j=box_left[1]  ;j<=box_left[3]  ;j++)
           for(int dir=0;dir<DIR;dir++)
           {
               Q[i][j][dir] = ghost_left[index];
               index ++;
           }
        //}
        
        
        // right ghost
        //if(myrank<cpusize-1)
        //{
           index = 0;
           for(int i=box_right[2]+1;i<=box_right[2]+2;i++)
           for(int j=box_right[1]  ;j<=box_right[3]  ;j++)
           for(int dir=0;dir<DIR;dir++)
           {
               Q[i][j][dir] = ghost_right[index];
               index ++;
           }
        //}
        
    update_BC_node();    
    
    set_SYNCH();
    MPI_Barrier(MPI_COMM_WORLD);
    return 1;
}



void update_BC_node(void)
{
	const int row=2;
	double buff[MAXDIM*row*DIR];
	int data_len;
	
	for(int rank=0;rank<cpusize;rank++)
	{
		int iis, iit, jjs, jjt;   
        my_range(rank, iis, iit, jjs, jjt);
/*        
        if(iis==2)
     	{
     		data_len = (jjt-jjs+1)*DIR*row;
     		if(myrank==rank)  // setup buff for send
     		{
     		   int index = 0;
     		   for(int i=iis;i<=iis+row-1;i++)
     		   for(int j=jjs;j<=jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
			   {
     		      buff[index] = Q[i][j][dir];
     		      index ++;
     		   }
     		}
     		MPI_Bcast( buff, data_len, MPI_DOUBLE, rank, MPI_COMM_WORLD);
     		//std::cout<<myrank<<" Complete broadcast of LEFT from "<<rank<<std::endl;
     		   int index = 0;
     		   for(int i=iis;i<=iis+row-1;i++)
     		   for(int j=jjs;j<=jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
     		   {
     		       Q[i][j][dir] = buff[index];
     		       index ++;
     		   }
     		
     	}
     	//MPI_Barrier(MPI_COMM_WORLD);
     	
     	
        if(iit==IDIM)
     	{
     		data_len = (jjt-jjs+1)*DIR*row;
     		if(myrank==rank)
     		{
     		   int index = 0;
     		   for(int i=iit-(row-1);i<=iit;i++)
     		   for(int j=jjs;j<=jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
			   {
     		      buff[index] = Q[i][j][dir];
     		      index ++;
     		   }
     		}
     		MPI_Bcast( buff, data_len, MPI_DOUBLE, rank, MPI_COMM_WORLD);
     		//std::cout<<myrank<<" Complete broadcast of RIGHT from "<<rank<<std::endl;
     		   int index = 0;
     		   for(int i=iit-(row-1);i<=iit;i++)
     		   for(int j=jjs;j<=jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
     		   {
     		       Q[i][j][dir] = buff[index];
     		       index ++;
     		   }
     		
     	}     	
     	//MPI_Barrier(MPI_COMM_WORLD);
*/     	
        if(jjs==2)
     	{
     		data_len = (iit-iis+1)*DIR*row;
     		if(myrank==rank)
     		{
     		   int index = 0;
     		   for(int i=iis;i<=iit;i++)
     		   for(int j=jjs;j<=jjs+(row-1);j++)
     		   for(int dir=0;dir<DIR;dir++)
			   {
     		      buff[index] = Q[i][j][dir];
     		      index ++;
     		   }
     		}
     		MPI_Bcast( buff, data_len, MPI_DOUBLE, rank, MPI_COMM_WORLD);
     		//std::cout<<myrank<<" Complete broadcast of BOTTOM from "<<rank<<std::endl;
     		   int index = 0;
     		   for(int i=iis;i<=iit;i++)
     		   for(int j=jjs;j<=jjs+(row-1);j++)
     		   for(int dir=0;dir<DIR;dir++)
     		   {
     		       Q[i][j][dir] = buff[index];
     		       index ++;
     		   }
     		
     	}
     	//MPI_Barrier(MPI_COMM_WORLD);

/*     	
        if(jjt==JDIM)
     	{
     		data_len = (iit-iis+1)*DIR*row;
     		if(myrank==rank)
     		{
     		   int index = 0;
     		   for(int i=iis;i<=iit;i++)
     		   for(int j=jjt-(row-1);j<jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
			   {
     		      buff[index] = Q[i][j][dir];
     		      index ++;
     		   }
     		}
     		MPI_Bcast( buff, data_len, MPI_DOUBLE, rank, MPI_COMM_WORLD);
     		//std::cout<<myrank<<" Complete broadcast of TOP from "<<rank<<std::endl;
     		   int index = 0;
     		   for(int i=iis;i<=iit;i++)
     		   for(int j=jjt-(row-1);j<jjt;j++)
     		   for(int dir=0;dir<DIR;dir++)
     		   {
     		       Q[i][j][dir] = buff[index];
     		       index ++;
     		   }
     		
     	}	
     	//MPI_Barrier(MPI_COMM_WORLD);	 
*/
	}
	
	
	
}
