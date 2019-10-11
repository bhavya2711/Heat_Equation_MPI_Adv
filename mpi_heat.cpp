#include <mpi.h>
#include <math.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
 #endif

 double generate2DHeat(long n, long global_i, long global_j);
 int check2DHeat(long n, long global_i, long global_j, double v, long k); //this function return 1 on correct. 0 on incorrect. Note that it may return 1 on incorrect. But a return of 0 means it is definitely incorrect

 #ifdef __cplusplus
}
#endif

void showHeat(double **H, int m, int n) {
  for (int i=0; i<m; ++i) {
    for (int j=0; j<n; ++j) {
      std::cout<<H[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
}




int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1;
  }

  // command line parameters
  MPI_Init(&argc,&argv);
  long n, K;
  n = atol(argv[1]);
  K = atol(argv[2]);

  int rank,size;
//rank and size
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

    //sqrt of size
   int root_size = sqrt(size);
    //chunksize
   long chunksize = n/root_size; //n/sqrt(p)

   int row_div = rank/root_size,col_div = rank%root_size;

   double**  H = new double*[chunksize];
   double**  G = new double*[chunksize];
   for(long i=0;i<chunksize;i++)
     {
     H[i] = new double[chunksize];
     G[i] = new double[itchunksize];

   }
//row start and end
   long rowstart = (row_div*chunksize);
   long colstart = (col_div*chunksize);
   long rowend = rowstart+chunksize;
   long colend = colstart+chunksize;

  for (long grow = rowstart,rset=0; grow<rowend; grow++,rset++) {
    for (long gcol= colstart,cset=0; gcol<colend; gcol++,cset++) {
       H[rset][cset] = generate2DHeat(n, grow,gcol);
    }
  }
  // write code here

    //all send and receive buffer
  double *rec_buffer1 = new double[chunksize];
  double *rec_buffer3 = new double[chunksize];
  double *rec_buffer2 = new double[chunksize];
  double *rec_buffer4 = new double[chunksize];
  double *send_buffer1 = new double[chunksize];
  double *send_buffer3 = new double[chunksize];
  double *send_buffer2 = new double[chunksize];
  double *send_buffer4 = new double[chunksize];

  int top,bottom,left,right;

  left =col_div?rank-1:-1;

  right = (col_div == (root_size-1))?-1:rank+1;

  top = rank-root_size;
  bottom = rank+root_size;
    
    //mpi status array
  MPI_Status status_array[4];
    
    //mpi request array
  MPI_Request request_array[8];
    //counter
  long count = 0;

    //start time using MPI_Wtime
  double start = MPI_Wtime();

  for (long it = 0; it<K; it++)
  {
     for(long ind =0,ite = 0;ind < chunksize;ind++)
       {
   rec_buffer4[ite] = H[ind][0];
   rec_buffer2[ite] = H[ind][chunksize-1];
   rec_buffer1[ite]  = H[0][ind];
   rec_buffer3[ite] = H[chunksize-1][ind];
   send_buffer4[ite] = H[ind][0];

   send_buffer2[ite] = H[ind][chunksize-1];

   send_buffer1[ite]  = H[0][ind];
   //send_buffer2
   send_buffer3[ite] = H[chunksize-1][ind];
   //send_buffer4
   ite++;
       }
     count = 0;
     if(top>=0){
         //top, then right, bottom, left
       MPI_Isend(send_buffer1,chunksize,MPI_DOUBLE,top,0,MPI_COMM_WORLD,&request_array[count]);
       count++;
     }
     if(right != -1){
       MPI_Isend(send_buffer2,chunksize,MPI_DOUBLE,right,1,MPI_COMM_WORLD,&request_array[count]);
       count++;
     }
     if(bottom < size){
       MPI_Isend(send_buffer3,chunksize,MPI_DOUBLE,bottom,2,MPI_COMM_WORLD,&request_array[count]);
       count++;
     }
     if(left != -1){
       MPI_Isend(send_buffer4,chunksize,MPI_DOUBLE,left,3,MPI_COMM_WORLD,&request_array[count]);
       count++;
     }
     for(long i=1;(i+1)<chunksize;i++)
       {
   for(long j=1;(j+1)<chunksize;j++)
     {
       G[i][j] = (H[i][j] + H[i-1][j] + H[i+1][j] + H[i][j-1] + H[i][j+1])/5;
     }
       }
     if(top >= 0){
         //top, right bottom,left
       MPI_Recv(rec_buffer1,chunksize,MPI_DOUBLE,top,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

     }
     if(right != -1){
       MPI_Recv(rec_buffer2,chunksize,MPI_DOUBLE,right,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

     }
     if(bottom < size){
       MPI_Recv(rec_buffer3,chunksize,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

     }
     if(left != -1){
       MPI_Recv(rec_buffer4,chunksize,MPI_DOUBLE,left,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

     }
     G[0][0] = (H[0][0] + rec_buffer1[0] + rec_buffer4[0] + H[0][1] + H[1][0])/5;
     G[0][chunksize-1] = (H[0][chunksize-1] + rec_buffer1[chunksize-1] + rec_buffer2[0]+H[0][chunksize-2]+H[1][chunksize-1])/5;
     G[chunksize-1][0] = (H[chunksize-1][0] + rec_buffer3[0] + rec_buffer4[chunksize-1]+H[chunksize-1][1]+H[chunksize-2][0])/5;
     G[chunksize-1][chunksize-1] = (H[chunksize-1][chunksize-1] + rec_buffer3[chunksize-1] +
          rec_buffer2[chunksize-1]+H[chunksize-1][chunksize-2]+H[chunksize-2][chunksize-1])/5;
     for(long i=1,j=1;i<chunksize-1;i++,j++)
   {
      G[0][i] = (H[0][i]+H[1][i]+rec_buffer1[i]+H[0][i-1]+H[0][i+1])/5;
      G[chunksize-1][i] = (H[chunksize-1][i]+H[chunksize-2][i]+rec_buffer3[i]+H[chunksize-1][i-1]+H[chunksize-1][i+1])/5;
      G[j][0] = (H[j][0]+rec_buffer4[j]+H[j][1]+H[j-1][0]+H[j+1][0])/5;
      G[j][chunksize-1] = (H[j][chunksize-1]+H[j-1][chunksize-1]+H[j+1][chunksize]+rec_buffer2[j]+H[j][chunksize-2])/5;
   }
    //showHeat(G[0], i, j);

      H = G;
      //showHeat(H[0], i, j);
      //Wait for all
      MPI_Waitall(count,request_array,status_array);
      //showHeat(H[0], i, j);


   }
  if(rank == 0)
    {
       double end = MPI_Wtime();
       cerr<<end-start<<endl;
    }
    //delete and free all
 for(long i=0;i<chunksize;i++)
 {
   delete[] H[i];
 }

  delete[] H;
  delete[] rec_buffer1;
  delete[] rec_buffer3;
  delete[] rec_buffer2;
  delete[] rec_buffer4;

//close the environment
  MPI_Finalize();

  return 0;
}
