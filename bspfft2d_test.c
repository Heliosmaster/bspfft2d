#include "bspfft2d.c"

/**
  * test on 2-Dimensional Fast Fourier Transform
  * author: Davide Taviani
  **/

int M,N,n0,n1;

void bspfft2d_test(){
  int p, pid, q, s, t, n0, n1, nlr, nlc, i, j;//,k1, *rho_np, *rho_p;
  double **a;
  
  bsp_begin(M*N);
  p=bsp_nprocs(); /* p=M*N */
  pid=bsp_pid();
  bsp_push_reg(&M,SZINT);
  bsp_push_reg(&N,SZINT);
  bsp_push_reg(&n0,SZINT);
  bsp_push_reg(&n1,SZINT);  
  bsp_sync();
  
   if (pid==0){
       /* printf("Please enter matrix row size n0:\n");
        scanf("%d",&n0);
        if(n0<2*M) bsp_abort("Error in input: n0 < 2M");
        printf("Please enter matrix column size n1:\n");
        scanf("%d",&n1);
        if(n1<2*N) bsp_abort("Error in input: n1 < 2N");*/
        n0=8;
        n1=8;
        for (q=0; q<p; q++){
            bsp_put(q,&M,&M,0,SZINT);
            bsp_put(q,&N,&N,0,SZINT);
            bsp_put(q,&n0,&n0,0,SZINT);
            bsp_put(q,&n1,&n1,0,SZINT);
        }
    }
    bsp_pop_reg(&n1);
    bsp_pop_reg(&n0); //not needed anymore 
    bsp_pop_reg(&N);
    bsp_pop_reg(&M);
    bsp_sync();
   

    /* Compute 2D processor numbering from 1D numbering */
    s= pid%M;  /* 0 <= s < M */
    t= pid/M;  /* 0 <= t < N */
    
   // printf("%d, s=%d,t=%d\n",pid,s,t);

    /* Allocate and initialize matrix */
    nlr=nloc(M,s,n0); /* number of local rows */
    nlc=nloc(N,t,n1); /* number of local columns */
    
    a= matallocd(nlr,2*nlc);
    //bsp_push_reg(a,2*nlr*nlc*SZDBL);
    if (s==0 && t==0){
        printf("2Dimensional FFT of a matrix with %d rows and %d columns\n",n0,n1);
        printf("using the %d by %d cyclic distribution\n",M,N);
    }
    
   
  // old: np = n/p , now every local fft is of the length nlc
    
    /*
    Matrix creation
    */
    double k=1.0; 
    
    for(i=0;i<nlr;i++)
      for(j=0; j<nlc; j++){
        a[i][2*j]= 2.0;
        a[i][2*j+1]= 1.0;
      }

a = bspfft2d(a,n0,n1,M,N,s,t,1);
   // bspfft2d(a,n0,n1,M,N,s,t,-1,w0,w,tw,rho_np,rho_p);
    

printm(a,nlr,nlc,s,t);

    //bsp_sync();
  
    matfreed(a);
    
    bsp_end();
  
}

int main(int argc, char **argv){
  
  bsp_init(bspfft2d_test, argc, argv);
  /*
  printf("Please enter number of processor rows M:\n");
    scanf("%d",&M);
    printf("Please enter number of processor columns N:\n");
    scanf("%d",&N);
    if (M*N > bsp_nprocs()){
        printf("Sorry, not enough processors available.\n"); 
        fflush(stdout);
        exit(1);
        }*/
        M=4;
        N=2;

    bspfft2d_test();
  
  exit(0);
  
}