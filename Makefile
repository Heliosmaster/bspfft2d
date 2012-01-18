all: transpose mxn	

transpose:
	mpicc -o transpose_bin bspfft2d_transpose_test.c -lbsponmpi -lm
mxn:
	mpicc -o fft2d bspfft2d_test.c -lbsponmpi -lm
