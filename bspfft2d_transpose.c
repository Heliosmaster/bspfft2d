#include "BSPedupack/bspedupack.c"
#include "BSPedupack/bspfft.c"

  /**
  * Prints the given m-by-n matrix a (mainly for debugging and code reduction)
  * s = current processor
  */
void printm(double **a, int m, int n,int s){
  int i,j;
  for(i=0;i<m;i++){
        printf("%d-%d: ",s,i);
      for(j=0;j<n;j++){
        printf("%f %f | ",a[i][2*j],a[i][2*j+1]);
      }
      printf("\n");
    }
}


  /**
  * Performs 1D FFT on every row of the nlr-by-nlc matrix a
  * sign = 1/-1 for forward/backward fft
  */
void bspfft_rows(double **a,int nlr,int nlc, int sign){
  double *w;
  int *rho;
  
  //memory allocation for the required weight & bit reversal tables
  rho = vecalloci(nlc);
  w=  vecallocd(nlc);
  
  // initialization of the tables
  ufft_init(nlc,w);
  bitrev_init(nlc,rho);
  
  int i,j;
  
  // permutes every row and performs 1D fft according to sign
  for(i=0;i<nlr;i++){
    permute(a[i],nlc,rho);
    ufft(a[i],nlc,sign,w);
  }
  
  // freeing memory
  vecfreei(rho);
  vecfreed(w);
  
  // if inverse FFT 
  if (sign==-1){
    double ninv = 1/(double)nlc;
    for(i=0;i<nlr;i++) for(j=0;j<nlc;j++){
      a[i][2*j] *= ninv;
      a[i][2*j+1] *= ninv;
    }
  }
}


void transpose(double **a,int nlr,int nlc, double *pa){
  int p = bsp_nprocs();
  int s = bsp_pid();
  double *tmp;
  int i,j,destproc,destrow,destindex,iglob;
  
  //allocates a temporary vector
  tmp = vecallocd(2);
  
  
  for(i=0;i<nlr;i++)
  for(j=0;j<nlc;j++){
    
    // computes the destination processor and rows starting from the local index j
    destproc = j%p;
    destrow = j/p;
    
    // computes global index
    iglob = s+i*p;
    
    // computes the actual index on the remote destination
    destindex = destrow*nlc+iglob;
    
    // stores in tmp the stuff to be sent
    tmp[0] = a[i][2*j];
    tmp[1] = a[i][2*j+1];
    
    // performs the actual comunication
    bsp_put(destproc,tmp,pa,destindex*2*SZDBL,2*SZDBL);
  }
  bsp_sync();
}


void bspfft2d_transpose(double **a,int nlr,int nlc, int sign, double *pa){
  //performs fft on the rows
  bspfft_rows(a,nlr,nlc,sign);
  
  //transpose and performs again fft on the rows
  transpose(a,nlr,nlc,pa);
  bspfft_rows(a,nlr,nlc,sign);
  
  //transpose it back
  transpose(a,nlr,nlc,pa);
}