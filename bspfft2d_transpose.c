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


void transpose(double **a,int nlr,int nlc, int nlc_dest, double *pm){
  int p = bsp_nprocs();
  int s = bsp_pid();
  double *tmp;
  int i,j,destproc,destrow,destindex,iglob,count=0;
  
  //allocates a temporary vector
  tmp = vecallocd(2);
  
  
  for(i=0;i<nlr;i++)
  for(j=0;j<nlc;j++){
    
    // computes the destination processor and rows starting from the local index j
    destproc = j%p;
    if (s == 0 && destproc != s) count++;
    destrow = j/p;
    
    // computes global index
    iglob = s+i*p;
    
    // computes the actual index on the remote destination
    destindex = destrow*nlc_dest+iglob;
    
    // stores in tmp the stuff to be sent
    tmp[0] = a[i][2*j];
    tmp[1] = a[i][2*j+1];
    
    //printf("%d: (%d,%d) i put a%d%d=%d into %d, row %d, index %d, total index %d\n",s,nlr,nlc,i,2*j,(int)a[i][2*j],destproc,destrow,iglob,destindex);
    // performs the actual comunication
    bsp_put(destproc,tmp,pm,destindex*2*SZDBL,2*SZDBL);
  }
  if (s == 0) printf("%d: Moving %d elements \n",s,count);
  bsp_sync();
}


void bspfft2d_transpose(double **a,int nlr,int nlc, int sign){
  //performs fft on the rows
  bspfft_rows(a,nlr,nlc,sign);
  
  //initialisation of the transpose
  double **t;
  int p = bsp_nprocs();
  
  /**
  * the transpose has a different size (generally): every proc has n0 columns (nlr*p=n0/p*p), and n1(=nlc)/p rows
  */
  int nlc_t = nlr*p;
  int nlr_t = nlc/p;

  //memory allocation
  t = matallocd(nlr_t,2*nlc_t);
  
  //set the pointer to the beginning of the transpose for communication purposes
  double *pt = t[0];
  double *pa = a[0];
  
  //i let everybody know about my new stuff
  bsp_push_reg(pt,2*nlr_t*nlc_t*SZDBL);
  bsp_push_reg(pa,2*nlr*nlc*SZDBL);
  bsp_push_reg(t,2*nlr_t*nlc_t*SZDBL);
  bsp_sync();
  
  //transpose a
  transpose(a,nlr,nlc,nlc_t,pt);
  
  // performs fft on the rows of the transpose of a
  bspfft_rows(t,nlr_t,nlc_t,sign);
  

  //transpose it back
  transpose(t,nlr_t,nlc_t,nlc,pa);
  

  // free memory
  matfreed(t);
  bsp_pop_reg(t);
  bsp_pop_reg(pa);
  bsp_pop_reg(pt);

}