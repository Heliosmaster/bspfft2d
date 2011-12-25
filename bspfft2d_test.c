#include "bspfft2d.c"

/**
  * test on 2-Dimensional Fast Fourier Transform
  * author: Davide Taviani
  **/

int M,N;

void bspfft2d_test(){
  bsp_begin(M*N);
  p=bsp_nprocs(); /* p=M*N */
  pid=bsp_pid();
  bsp_push_reg(&M,SZINT);
  bsp_push_reg(&N,SZINT);
  bsp_push_reg(&n,SZINT);
  bsp_sync();
  
   if (pid==0){
        printf("Please enter matrix size n:\n");
        scanf("%d",&n);
        for (q=0; q<p; q++){
            bsp_put(q,&M,&M,0,SZINT);
            bsp_put(q,&N,&N,0,SZINT);
            bsp_put(q,&n,&n,0,SZINT);
        }
    }
    bsp_sync();
    bsp_pop_reg(&n); /* not needed anymore */
    bsp_pop_reg(&N);
    bsp_pop_reg(&M);

    /* Compute 2D processor numbering from 1D numbering */
    s= pid%M;  /* 0 <= s < M */
    t= pid/M;  /* 0 <= t < N */

    /* Allocate and initialize matrix */
    nlr=  nloc(M,s,n); /* number of local rows */
    nlc=  nloc(N,t,n); /* number of local columns */
    
    a= matallocd(nlr,nlc);
  
}

int main(int argc, char **argv){
  
  bsp_init(bspfft2d_test, argc, argv);
  
  printf("Please enter number of processor rows M:\n");
    scanf("%d",&M);
    printf("Please enter number of processor columns N:\n");
    scanf("%d",&N);
    if (M*N > bsp_nprocs()){
        printf("Sorry, not enough processors available.\n"); 
        fflush(stdout);
        exit(1);
    }

    bspfft2d_test();
  
  exit(0);
  
}