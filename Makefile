all: transpose mxn	

transpose:
	mpcc -o transpose_bin bspfft2d_transpose_test.c -lbsponmpi -lm -O3
mxn:
	mpcc -o fft2d bspfft2d_test.c -lbsponmpi -lm -O3
