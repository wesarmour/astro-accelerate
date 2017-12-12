#ifndef ASTROACCELERATE_MSD_LEGACY_H_
#define ASTROACCELERATE_MSD_LEGACY_H_

extern int MSD_linear_approximation(float *d_input, float *d_MSD_T, int nTaps, int nDMs, int nTimesamples, int offset);
extern int MSD_LA_Nth(float *d_input, float *d_bv_in, float *d_MSD_T, float *d_MSD_DIT, int nTaps, int nDMs, int nTimesamples, int offset, int DIT_value);

#endif
