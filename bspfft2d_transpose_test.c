#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "bspfft2d_transpose.c"

#define NITERS 1
#define MEGA 1000000.0

int P,n0,n1;

void bspfft2d_transpose_test(n0,n1){
  int i,j,p,s,nlr,nlc,it;
  double **a, *Error;
  double time0,time1;
  
  bsp_begin(P);
  p= bsp_nprocs();
  s= bsp_pid();
  
  bsp_push_reg(&n0,SZINT);
  bsp_push_reg(&n1,SZINT);
  
  bsp_sync();
  
  if (s==0){
      for (j=0; j<p; j++){ 
        bsp_put(j,&n0,&n0,0,SZINT);
        bsp_put(j,&n1,&n1,0,SZINT);
      }
   
  }
    
  bsp_sync();
  
  // number of local rows (cyclic distribution w/ powers of 2) and columns (the length of each fft)
  nlr = n0/p;
  nlc = n1; 
   
  // memory allocation and variable registration
  a = matallocd(nlr,2*nlc); 
  bsp_push_reg(a,2*nlr*nlc*SZDBL);

  
  if (s==0) printf("2D FFT of a matrix %d-by-%d using %d processors\n",n0,n1,p);
  int k=0;
  // matrix creation
  for(i=0;i<nlr;i++)
  for(j=0;j<nlc;j++){
    a[i][2*j]= 10*s+k;//2.0;
    a[i][2*j+1]= 10*s+k;//1.0;
    k++;
  }
  
  // pointer to the beginning of the matrix, used as starting point for communication
  
  bsp_sync();
  time0 = bsp_time();
  
  for(it=0;it<NITERS;it++){
   printm(a,nlr,nlc,s);
    sleep(1);
    printf("-----\n");
    sleep(1);
    bspfft2d_transpose(a,nlr,nlc,1);    // forward 2D fft
   // printm(a,nlr,nlc,s);
//  bspfft2d_transpose(a,nlr,nlc,-1);   //backward 2d fft
  }

  bsp_sync();
  time1 = bsp_time();
 
//printm(a,nlr,nlc,s);  

  printf("%d: It took exactly %f seconds for each FFT\n",s,(time1-time0)/NITERS);

  
  //free memory and de-register variables
  
  matfreed(a);
  bsp_pop_reg(a);
  bsp_pop_reg(&n0);
bsp_pop_reg(&n1);

  bsp_end();
}

int main(int argc, char **argv){
  int i,n0,n1;
  bsp_init(bspfft2d_transpose_test, argc, argv);
    if (argc>0){
      P = atoi(argv[1]);
      if(P>bsp_nprocs()) bsp_abort("**Sorry, not enough processors available.**\n");
      n0 = atoi(argv[2]);
      n1 = atoi(argv[3]);
    }

  bspfft2d_transpose_test(n0,n1);
  
  exit(0);
  
}
