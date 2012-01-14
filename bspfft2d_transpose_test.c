#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "bspfft2d_transpose.c"

int P;

void bspfft2d_transpose_test(){
  int n0,n1,i,j,p,s,nlr,nlc;
  double **a;
  
  bsp_begin(P);
  p= bsp_nprocs();
  s= bsp_pid();
  
  bsp_push_reg(&n0,SZINT);
  bsp_push_reg(&n1,SZINT);
  
  bsp_sync();
  
  if (s==0){
    //printf("Please enter matrix row size n0:\n");
    // scanf("%d",&n0);
    //  if(n0<2*M) bsp_abort("Error in input: n0 < 2M");
    //  printf("Please enter matrix column size n1:\n");
    //  scanf("%d",&n1);
    //  if(n1<2*N) bsp_abort("Error in input: n1 < 2N");
        n0=4;
        n1=4;
      for (j=0; j<p; j++){ 
        bsp_put(j,&n0,&n0,0,SZINT);
        bsp_put(j,&n1,&n1,0,SZINT);
      }
      //if(n0<P) bsp_abort("Too small for this many procs!");
  }
    
  bsp_sync();
  
  // number of local rows (cyclic distribution w/ powers of 2)
  nlr = n0/p;
  nlc = n1; 
   
  a = matallocd(nlr,2*nlc);
  bsp_push_reg(a,2*nlr*nlc*SZDBL);
  
  //matrix creation
  
  for(i=0;i<nlr;i++)
  for(j=0;j<nlc;j++){
    a[i][2*j]= 10*s+nlc*i+j;
    a[i][2*j+1]= 10*s+nlc*i+j;
  }
  
  double *pa = a[0];
  bsp_push_reg(pa,2*nlr*nlc*SZDBL);
  
  bsp_sync();
  
 /* 
  printf("before:\n");
  printm(a,nlr,nlc,s);
  */
  bspfft2d_transpose(a,nlr,nlc,1,pa);
 // bspfft2d_transpose(a,nlr,nlc,-1,pa);

  
  //sleep(1);
// printf("after:\n");
  
  printm(a,nlr,nlc,s);
  
  matfreed(a);
  bsp_pop_reg(pa);
  bsp_pop_reg(a);
  bsp_pop_reg(&n0);
  bsp_pop_reg(&n1);

  bsp_end();
}

int main(int argc, char **argv){
  
  bsp_init(bspfft2d_transpose_test, argc, argv);
  /*
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%d",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, not enough processors available.\n"); fflush(stdout);
        exit(1);
    }
    */
    P = 2;

  bspfft2d_transpose_test();
  
  exit(0);
  
}
