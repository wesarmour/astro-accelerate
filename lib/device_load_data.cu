//#include "AstroAccelerate/device_dedispersion_kernel.h"

//extern "C" void load_data(int i, float *device_pointer, float *host_pointer, size_t size, int nsamp, int maxshift, int nchans, int t_processed_s, int t_processed_c, float *dmshifts);

//{{{ load_data_from_host_to_device

void load_data(int i, int *inBin, float *device_pointer, float *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts) {

	//{{{ Copy data and set up the GPU constants/variables.
	if(i==-1) {
		int length=(t_processed+maxshift);
		size_t size=nchans*length*sizeof(float);
		cudaMemcpyToSymbol(dm_shifts, dmshifts, nchans * sizeof(float));
		cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(i_nchans, &nchans, sizeof(int));
		cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
		cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
	} else if (i > 0) {
		int length=(t_processed+maxshift);
		cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
		cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
	}
	//}}}
}

//}}}