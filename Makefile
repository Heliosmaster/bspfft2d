all: transpose mxn	

transpose:
	mpicc -o transpose bspfft2d_transpose_test.c -lbsponmpi -lm
mxn:
	mpicc -o fft2d bspfft2d_test.c -lbsponmpi -lm

	
