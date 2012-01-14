#include "BSPedupack/bspedupack.c"
#include "BSPedupack/bspfft.c"


void printm(double **a, int m, int n,int s){
  int i,j;
  for(i=0;i<m;i++){
        printf("%d-%d: ",s,i);
      for(j=0;j<n;j++){
        printf("%d %d | ",(int)a[i][2*j],(int)a[i][2*j+1]);
      }
      printf("\n");
    }
}

void bspfft2d_rows(double **a,int nlr,int nlc, int sign){
  double *w;
  int *rho;
  rho = vecalloci(nlc);
  w=  vecallocd(nlc);
  ufft_init(nlc,w);
  bitrev_init(nlc,rho);
  int i;
  for(i=0;i<nlr;i++){
    permute(a[i],nlc,rho);
    ufft(a[i],nlc,sign,w);
  }
  vecfreed(w); 
}

void bspfft2d_move(double **a,int nlr,int nlc, double *pa){
  int p = bsp_nprocs();
  int s = bsp_pid();
  double *tmp;
  int i,j,destproc,destrow,destindex,iglob;
  
  
  tmp = vecallocd(2);
  
  for(i=0;i<nlr;i++)
  for(j=0;j<nlc;j++){
    iglob = s+i*p;
    destproc = j%p;
    destrow = j/p;
    destindex = destrow*nlc+iglob;
    tmp[0] = a[i][2*j];
    tmp[1] = a[i][2*j+1];
//    printf(" %d: putting a%d%d=%d in %d, row %d col %d, index %d\n",s,i,j,(int)a[i][2*j],destproc,destrow,iglob,destindex);
bsp_put(destproc,tmp,pa,destindex*2*SZDBL,2*SZDBL);
  }
  bsp_sync();
}


void bspfft2d_transpose(double **a,int nlr,int nlc, int sign, double *pa){
  bspfft2d_rows(a,nlr,nlc,sign);
bspfft2d_move(a,nlr,nlc,pa);
bspfft2d_rows(a,nlr,nlc,sign);
bspfft2d_move(a,nlr,nlc,pa);
}