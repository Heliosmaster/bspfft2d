#include "BSPedupack/bspedupack.c"
//#include "BSPedupack/bspfft.c"

/**
  * 2-Dimensional Fast Fourier Transform
  * author: Davide Taviani
  * x
  * M, N # of procs row and columns
  * P(s,t) current proc
  * sign,....
  **/

//double **a; 
int nloc(int p, int s, int n){
  /* Compute number of local components of processor s for vector
  of length n distributed cyclically over p processors. */
  return  (n+p-s-1)/p ;
} /* end nloc */
  
  
  /****************** Sequential functions ********************************/
void ufft(double *x, int n, int sign, double *w){
    
    /* This sequential function computes the unordered discrete Fourier
      transform of a complex vector x of length n, stored in a real array
      of length 2n as pairs (Re x[j], Im x[j]), 0 <= j < n.
      n=2^m, m >= 0.
    If sign = 1, then the forward unordered dft FRx is computed;
    if sign =-1, the backward unordered dft conjg(F)Rx is computed,
      where F is the n by n Fourier matrix and R the n by n bit-reversal
      matrix. The output overwrites x.
      w is a table of n/2 complex weights, stored as pairs of reals,
      exp(-2*pi*i*j/n), 0 <= j < n/2,
      which must have been initialized before calling this function.
      */
    
    int k, nk, r, rk, j, j0, j1, j2, j3;
    double wr, wi, taur, taui;
    
    for(k=2; k<=n; k *=2){
      nk= n/k;
      for(r=0; r<nk; r++){
        rk= 2*r*k;
        for(j=0; j<k; j +=2){
          wr= w[j*nk];
          if (sign==1) {
            wi= w[j*nk+1];
          } else {
            wi= -w[j*nk+1];
          }
          j0= rk+j;
          j1= j0+1;
          j2= j0+k;
          j3= j2+1;
          taur= wr*x[j2] - wi*x[j3];
          taui= wi*x[j2] + wr*x[j3];
          x[j2]= x[j0]-taur;
          x[j3]= x[j1]-taui;
          x[j0] += taur;
          x[j1] += taui;
        }
      }
    }
    
  } /* end ufft */
  
  void ufft_init(int n, double *w){
    
    /* This function initializes the n/2 weights to be used
      in a sequential radix-2 FFT of length n.
      n=2^m, m >= 0.
      w is a table of n/2 complex weights, stored as pairs of reals,
      exp(-2*pi*i*j/n), 0 <= j < n/2.
      */
    
    int j, n4j, n2j;
    double theta;
    
    if (n==1)
    return;
    theta= -2.0 * M_PI / (double)n;
    w[0]= 1.0;
    w[1]= 0.0;
    if (n==4){
      w[2]=  0.0;
      w[3]= -1.0;
    } else if (n>=8) {
      /* weights 1 .. n/8 */
      for(j=1; j<=n/8; j++){
        w[2*j]=   cos(j*theta);
        w[2*j+1]= sin(j*theta);
      }
      /* weights n/8+1 .. n/4 */
      for(j=0; j<n/8; j++){
        n4j= n/4-j;
        w[2*n4j]=   -w[2*j+1];
        w[2*n4j+1]= -w[2*j];
      }
      /* weights n/4+1 .. n/2-1 */
      for(j=1; j<n/4; j++){
        n2j= n/2-j;
        w[2*n2j]=   -w[2*j];
        w[2*n2j+1]=  w[2*j+1];
      }
    }
    
  } /* end ufft_init */
  
  void twiddle(double *x, int n, int sign, double *w){
    
    /* This sequential function multiplies a complex vector x
      of length n, stored as pairs of reals, componentwise
      by a complex vector w of length n, if sign=1, and
      by conjg(w), if sign=-1. The result overwrites x.
      */
    
    int j, j1;
    double wr, wi, xr, xi;
    
    for(j=0; j<2*n; j +=2){
      j1= j+1;
      wr= w[j];
      if (sign==1) {
        wi= w[j1];
      } else {
        wi= -w[j1];
      }
      xr= x[j];
      xi= x[j1];
      x[j]=  wr*xr - wi*xi;
      x[j1]= wi*xr + wr*xi;
    }
    
  } /* end twiddle */
  
  void twiddle_init(int n, double alpha, int *rho, double  *w){
    
    /* This sequential function initializes the weight table w
      to be used in twiddling with a complex vector of length n,
      stored as pairs of reals.
      n=2^m, m >= 0.
      alpha is a real shift parameter.
      rho is the bit-reversal permutation of length n,
      which must have been initialized before calling this function.
      The output w is a table of n complex values, stored as pairs of reals,
      exp(-2*pi*i*rho(j)*alpha/n), 0 <= j < n.
      */
    
    int j;
    double theta;
    
    theta= -2.0 * M_PI * alpha / (double)n;
  //  printf("%d: theta=%f\n",bsp_pid(),theta);
    
    for(j=0; j<n; j++){
      w[2*j]=   cos(rho[j]*theta);
      w[2*j+1]= sin(rho[j]*theta);
    }
    
  } /* end twiddle_init */
  
  void permute(double *x, int n, int *sigma){
    
    /* This in-place sequential function permutes a complex vector x
      of length n >= 1, stored as pairs of reals, by the permutation sigma,
      y[j] = x[sigma[j]], 0 <= j < n.
      The output overwrites the vector x.
      sigma is a permutation of length n that must be decomposable
      into disjoint swaps.
      */
    
    int j, j0, j1, j2, j3;
    double tmpr, tmpi;
    
    for(j=0; j<n; j++){
      if (j<sigma[j]){
        /* swap components j and sigma[j] */
        j0= 2*j;
        j1= j0+1;
        j2= 2*sigma[j];
        j3= j2+1;
        tmpr= x[j0];
        tmpi= x[j1];
        x[j0]= x[j2];
        x[j1]= x[j3];
        x[j2]= tmpr;
        x[j3]= tmpi;
      }
    }
    
  } /* end permute */
  
  void bitrev_init(int n1, int *rho){
    
    /* This function initializes the bit-reversal permutation rho
      of length n, with n=2^m, m >= 0.
      */
    
    int j;
    unsigned int n, rem, val, k, lastbit, one=1;
    
    if (n1==1){
      rho[0]= 0;
      return;
    }
    n= n1;
    for(j=0; j<n1; j++){
      rem= j; /* j= (b(m-1), ... ,b1,b0) in binary */
      val= 0;
      for (k=1; k<n; k <<= 1){
        lastbit= rem & one; /* lastbit = b(i) with i= log2(k) */
        rem >>= 1;          /* rem = (b(m-1), ... , b(i+1)) */
        val <<= 1;
        val |= lastbit;     /* val = (b0, ... , b(i)) */
      }
      rho[j]= (int)val;
    }
    
  } /* end bitrev_init */
  
  /****************** Parallel functions ********************************/
  int k1_init(int n1, int N,int nlc){
    
    /* This function computes the largest butterfly size k1 of the first
      superstep in a parallel FFT of length n on p processors with p < n.
      */
    
    int c, k1;
    
    for(c=1; c<N; c *=nlc)
        ;
    k1= n1/c;
    
    return k1;
    
  } /* end k1_init */
  
  void bspredistr(double *x, int i, int n0, int n1,
  int M, int N, int pid, int s, int t,
  int c0, int c1,char rev, int *rho_p, double *pa){
    
    /* This function redistributes the complex vector x of length n,
       stored as pairs of reals, from group-cyclic distribution
       over p processors with cycle c0 to cycle c1, where
       c0, c1, p, n are powers of two with 1 <= c0 <= c1 <= p <= n.
       s is the processor row, t is the processor column and pid = P(s,t)
      
       If rev=true, the function assumes the processor numbering
       is bit reversed on input.
      
       rho_p is the bit-reversal permutation of length p.
    */
      
    
    double *tmp;
    int j0, j2, j, jglob, ratio, size,nlr,nlc,nlc_dest;
    int npackets, destproc, destindex, r;
    
    nlc= nloc(N,t,n1);
    ratio= c1/c0;
    size= MAX(nlc/ratio,1);
    npackets= nlc/size;
    tmp= vecallocd(2*size);
    nlr=nloc(M,s,n0);    
    
    if (rev) {
      j0= rho_p[t]%c0;
      j2= rho_p[t]/c0;
    } else {
      j0= t%c0;
      j2= t/c0;
    }    
    for(j=0; j<npackets; j++){
      jglob= j2*c0*nlc + j*c0 + j0;
      
      destproc = (jglob/(c1*nlc))*c1 + jglob%c1;
      //if(s==1 && t == 0) printf("P(%d,%d), destproc = %d,%d\n",s,t,destproc,s+M*destproc);
      // we now have P(s,destproc), we have to convert from 2D to 1D numbering
      destproc= s+M*destproc;
    
      // compute the number of local columns for the destproc
      nlc_dest = nloc(N,(jglob/(c1*nlc))*c1 + jglob%c1,n1);
      
      /*
      * the second term of the sum is because we don't really know
      * the address of a[i] in the destproc, so we start from the
      * beginning of a and jump
      */
      destindex= (jglob%(c1*nlc))/c1+i*nlc_dest;
      for(r=0; r<size; r++){
        tmp[2*r]=x[2*(j+r*ratio)];
        tmp[2*r+1]= x[2*(j+r*ratio)+1];
      }
      //printf("I am trying to put stuff in proc %d, at address %d\n",destproc,destindex*2*SZDBL);
      bsp_put(destproc,tmp,pa,destindex*2*SZDBL,size*2*SZDBL);
    }
    //bsp_sync(); // deleted to avoid unnecessary syncs
    vecfreed(tmp);
    
  } /* end bspredistr nlr=  nloc(M,s,n0); np=  nloc(N,t,n1); */
  
  void bspfft1d(double **a, int n0, int n1, int M, int N, int pid, int s,
                int t, int sign, double *w0, double *w,
                double *tw, int *rho_np, int *rho_p, double *pa){
    
    /* This parallel function computes the discrete Fourier transform on the rows
      of a complex matrix a:
      
      the n0 rows of a are vectors of length n1=2^m, m >= 1, stored in a real array
      of length 2n1 as pairs (Re x[j], Im x[j]), 0 <= j < n.
      a must have been registered before calling this function.
      p is the number of processors, p=2^q, 0 <= q < m.
      s is the processor number, 0 <= s < p.
    The function uses three weight tables:
     w0 for the unordered fft of length k1,
      w  for the unordered fft of length n/p,
      tw for a number of twiddles, each of length n/p.
    The function uses two bit-reversal permutations:
    rho_np of length n/p,
      rho_p of length p.
      The weight tables and bit-reversal permutations must have been
      initialized before calling this function.
      If sign = 1, then the dft is computed,
      y[k] = sum j=0 to n-1 exp(-2*pi*i*k*j/n)*x[j], for 0 <= k < n.
      If sign =-1, then the inverse dft is computed,
      y[k] = (1/n) sum j=0 to n-1 exp(+2*pi*i*k*j/n)*x[j], for 0 <= k < n.
      Here, i=sqrt(-1). The output vector y overwrites x. (in-place fft)
    */
    char rev;
    int nlc, nlr, k1, r, c0, c, ntw, j,i;
    double ninv;
   
    nlr=  nloc(M,s,n0); /* number of local rows */
    nlc=  nloc(N,t,n1); /* number of local columns = length of 1d fft that every proc has */
    printf("(%d,%d): %d,%d\n",s,t,nlr,nlc);
    k1= k1_init(n1,N,nlc);
   
    for(i=0;i<nlr;i++){// do this stuff for every row
      permute(a[i],nlc,rho_np);
      rev= TRUE;
      for(r=0; r<nlc/k1; r++) ufft(&a[i][2*r*k1],k1,sign,w0);
    }
  
    c0= 1;
    ntw= 0;
    for (c=k1; c<=N; c *=nlc){
      for(i=0;i<nlr;i++){
        bspredistr(a[i],i,n0,n1,M,N,pid,s,t,c0,c,rev,rho_p,pa);        
      }
      bsp_sync();  //sync is done only after every row has been redistributed
      rev= FALSE;
      for(i=0;i<nlr;i++){ //do this stuff for every row
        twiddle(a[i],nlc,sign,&tw[2*ntw*nlc]); 
        ufft(a[i],nlc,sign,w);
      }      
      c0= c;
      ntw++;
    }
    
    if (sign==-1){
      ninv= 1 / (double)n1;
      for(i=0;i<nlr;i++) { //do this for every row
        for(j=0; j<2*nlc; j++)  a[i][j] *= ninv;
      } 
    }
  } /* end bspfft */
  
  void bspfft1d_init(int n1, int N, int s, int t, double *w0, double *w, double *tw,
  int *rho_np, int *rho_p){
    
    /* This parallel function initializes all the tables used in the FFT. */
    
    int nlc, k1, ntw, c;
    double alpha;
    
    nlc= nloc(N,t,n1);
    bitrev_init(nlc,rho_np);
    bitrev_init(N,rho_p);
    
    k1= k1_init(n1,N,nlc);
    ufft_init(k1,w0);
    ufft_init(nlc,w);
    
    ntw= 0;
    for (c=k1; c<=N; c *=nlc){
      alpha= (t%c) / (double)(c);
      twiddle_init(nlc,alpha,rho_np,&tw[2*ntw*nlc]);
      ntw++;
    }
    
  } /* end bspfft_init */
  
//**a matrix to be 2d-ffted
void bspfft2d(double **a, int n0, int n1, int M, int N, int s,
              int t,int sign, double *w0, double *w, double *tw, int *rho_np, int *rho_p){
  
    double *pa;
  int nlr, nlc, pid,i,j;
  
  nlr=  nloc(M,s,n0); /* number of local rows */
  nlc=  nloc(N,t,n1); /* number of local columns */
  pid = bsp_pid();
  
  pa= (nlr>0 ? a[0] : NULL);
  bsp_push_reg(pa,nlr*nlc*SZDBL);
  bsp_sync();
  bspfft1d(a,n0,n1,M,N,pid,s,t,sign,w0,w,tw,rho_np,rho_p,pa);
  
  //todo fft1d on the cols
  
  pa = a[0]; // this assignment is a workaround (otherwise it complains about bsp_pop_reg without a bsp_push_reg )
  bsp_pop_reg(pa);
  bsp_sync();  
}