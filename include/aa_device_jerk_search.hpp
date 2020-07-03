#ifndef ASTRO_ACCELERATE_AA_DEVICE_JERK_SEARCH_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_JERK_SEARCH_HPP

namespace astroaccelerate {

  /** \brief Function that performs analysis component on the GPU. */  
  int jerk_search(float ***dedispersed_data, JERK_Plan user_plan, float *dm_low, float *dm_step, int *list_of_ndms, float sampling_time, int *inBin, int nRanges);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP






