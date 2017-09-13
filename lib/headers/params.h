//P100 8,14,12,40
#define ACCMAX 350
#define ACCSTEP 11
#define CARD 0
#define NOPSSHIFT 5
#define NOPSLOOP 3
#define NDATAPERLOOP 1
#define BINDIVINT 8
#define BINDIVINF 32
#define CT 32
#define CF 8
#define NOPS 4.0
#define STATST 128
#define STATSLOOP 8

//Added by Karel Adamek
#define WARP 32
#define HALF_WARP 16
#define MSD_ELEM_PER_THREAD 8
#define MSD_WARPS_PER_BLOCK 16
#define THR_ELEM_PER_THREAD 4
#define THR_WARPS_PER_BLOCK 4
#define PD_NTHREADS 512
#define PD_NWINDOWS 2
#define PD_MAXTAPS 32
#define PD_SMEM_SIZE 1280
#define PD_FIR_ACTIVE_WARPS 2
#define PD_FIR_NWINDOWS 2
#define MIN_DMS_PER_SPS_RUN 64
#define MSD_PW_NTHREADS 512

/**** FDAS parameters ******/
/*Params for benchmarks */
#define SLIGHT 299792458.0
#define RADIX 1
#define NEXP 10
#define POTWO (1 << NEXP)
#define KERNLEN RADIX*POTWO
#define ACCEL_STEP (float)(2.0) //1 //default acceleration step
#define ACCEL_STEP_R (float)(1.0f/ACCEL_STEP)
#define ZMAX 3
#define NKERN (ZMAX + 1)
//#define ZLO  -(int)((ZMAX/ACCEL_STEP) )
#define TBSIZEX 32
#define TBSIZEY 1
#define PTBSIZEX 64
#define PTBSIZEY 1
//custom fft params (K. Adamek)
#define TAPS 8

#define UNROLLS 8
#define SNUMREG 8
#define SDIVINT 14
#define SDIVINDM 40
#define SFDIVINDM 40.0f

// for corner turn in shared memory corner_turn_SM(...)
#define CT_NTHREADS 512
#define CT_ROWS_PER_WARP 2
#define CT_CORNER_BLOCKS 1

// for periodicity harmonic summing
#define PHS_NTHREADS 64

// for power and interbin calculation
#define PAI_NTHREADS 512