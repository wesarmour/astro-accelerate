//Added by Karel Adamek

#include <vector>

#include "headers/params.h"
#include "headers/device_BC_plan.h"
#include "device_SPS_long_kernel.cu"

void Assign_parameters(int f, std::vector<PulseDetection_plan> *PD_plan, int *decimated_timesamples, int *dtm, int *iteration, int *nBoxcars, int *nBlocks, int *output_shift, int *shift, int *startTaps, int *unprocessed_samples, int *total_ut){
	*decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
	*dtm                   = PD_plan->operator[](f).dtm;
	*iteration             = PD_plan->operator[](f).iteration;
	*nBoxcars              = PD_plan->operator[](f).nBoxcars;
	*nBlocks               = PD_plan->operator[](f).nBlocks;
	*output_shift          = PD_plan->operator[](f).output_shift;
	*shift                 = PD_plan->operator[](f).shift;           
	*startTaps             = PD_plan->operator[](f).startTaps; 
	*unprocessed_samples   = PD_plan->operator[](f).unprocessed_samples;
	*total_ut              = PD_plan->operator[](f).total_ut;	
}

void PD_SEARCH_LONG_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}


int PD_SEARCH_LONG_BLN_IF(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	if(nBlocks>0) PD_GPU_1st_float1_BLN_IF<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		if( (f%2) == 0 ) {
			if(nBlocks>0) PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			if(nBlocks>0) PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
	}

	unprocessed_samples = unprocessed_samples*(1<<(max_iteration-1));
	return(unprocessed_samples);
}




int PD_SEARCH_LONG_BLN_IF_LINAPPROX(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	if(nBlocks>0) PD_GPU_1st_float1_BLN_IF_LA<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		if( (f%2) == 0 ) {
			if(nBlocks>0) PD_GPU_Nth_float1_BLN_IF_LA<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			if(nBlocks>0) PD_GPU_Nth_float1_BLN_IF_LA<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
	}

	unprocessed_samples = unprocessed_samples*(1<<(max_iteration-1));
	return(unprocessed_samples);
}




