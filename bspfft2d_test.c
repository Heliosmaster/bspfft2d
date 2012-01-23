#include "bspfft2d.c"

/**
  * test on 2-Dimensional Fast Fourier Transform
  * author: Davide Taviani
  **/

int M,N,n0,n1;

void bspfft2d_test(n0,n1){
  int p, pid, q, s, t, nlr, nlc, i, j;
  double **a,*times, tim,time0,time1;
  
  bsp_begin(M*N);
  p=bsp_nprocs(); /* p=M*N */
  pid=bsp_pid();
  times = vecallocd(p);
  bsp_push_reg(times,p*SZDBL);
  bsp_push_reg(&M,SZINT);
  bsp_push_reg(&N,SZINT);
  bsp_push_reg(&n0,SZINT);
  bsp_push_reg(&n1,SZINT);  
  bsp_sync();
  
   if (pid==0){
        for (q=0; q<p; q++){
            bsp_put(q,&M,&M,0,SZINT);
            bsp_put(q,&N,&N,0,SZINT);
            bsp_put(q,&n0,&n0,0,SZINT);
            bsp_put(q,&n1,&n1,0,SZINT);
        }
    }
    bsp_pop_reg(&n1);
    bsp_pop_reg(&n0);
    bsp_pop_reg(&N);
    bsp_pop_reg(&M);
    bsp_sync();
   

    /* Compute 2D processor numbering from 1D numbering */
    s= pid%M;  /* 0 <= s < M */
    t= pid/M;  /* 0 <= t < N */
    

    /* Allocate and initialize matrix */
    nlr=nloc(M,s,n0); /* number of local rows */
    nlc=nloc(N,t,n1); /* number of local columns */
    
    a= matallocd(nlr,2*nlc);
    if (s==0 && t==0){
        printf("2Dimensional FFT of a matrix with %d rows and %d columns\n",n0,n1);
        printf("using the %d by %d cyclic distribution\n",M,N);
    }
    
   
  // old: np = n/p , now every local fft is of the length nlc
    
    /*
    * Matrix creation
    */
    
    for(i=0;i<nlr;i++)
      for(j=0; j<nlc; j++){
        a[i][2*j]= 2.0;
        a[i][2*j+1]= 1.0;
      }

  bsp_sync();
  time0 = bsp_time();
  


  a=bspfft2d(a,n0,n1,M,N,s,t,1);
  a=bspfft2d(a,n0,n1,M,N,s,t,-1);
    
  bsp_sync();
  time1=bsp_time();
  
  tim = time1-time0;
  bsp_put(0,&tim,times,pid*SZDBL,SZDBL);
  bsp_sync();
  
  double tia = 0;
  
  if(pid==0){
    tim = 0;
    tia = 0;
    for(i=0;i<p;i++){
      tia += times[i];
      if (times[i]>tim) tim=times[i];
    }
    tia = tia/p;
    printf("Average time: %f\n",tia);
    printf("Maximum time: %f\n",tim);
  }

  
    matfreed(a);
    vecfreed(times);
    bsp_pop_reg(times);
    bsp_end();
  
}

int main(int argc, char **argv){
  
  bsp_init(bspfft2d_test, argc, argv);
  if (argc>0){
      M = atoi(argv[1]);
      N = atoi(argv[2]);
      if (M*N > bsp_nprocs()) bsp_abort("Sorry, not enough processors.\n");
      n0 = atoi(argv[3]);
      n1 = atoi(argv[4]);
    }
    
    bspfft2d_test(n0,n1);
  
  exit(0);
  
}