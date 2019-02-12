#ifndef ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP
#define ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP

#include "aa_pipeline_api.hpp"
#include "aa_py_filterbank_metadata.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_pipeline_api<unsigned short>* aa_py_pipeline_api(const aa_py_filterbank_metadata_struct metadata, unsigned short const*const input_data, const int card_number);
      void aa_py_pipeline_api_delete(aa_pipeline_api<unsigned short> const*const obj);
      bool aa_py_pipeline_api_bind_ddtr_plan(aa_pipeline_api<unsigned short> *const obj, aa_ddtr_plan const*const plan);

      aa_ddtr_strategy* aa_py_pipeline_api_ddtr_strategy(aa_pipeline_api<unsigned short> *const obj);
      
      bool aa_py_pipeline_api_bind_analysis_plan(aa_pipeline_api<unsigned short> *const obj, aa_analysis_plan const*const plan);
      aa_analysis_strategy* aa_py_pipeline_api_analysis_strategy(aa_pipeline_api<unsigned short> *const obj);

      bool aa_py_pipeline_api_bind_periodicity_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const int nHarmonics, const int export_powers, const bool candidate_algorithm, const bool enable_msd_baseline_noise);

      bool aa_py_pipeline_api_bind_fdas_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const int num_boots, const int num_trial_bins, const int navdms, const float narrow, const float wide, const int nsearch, const float aggression, const bool enable_msd_baseline_noise);
      
      bool aa_py_pipeline_api_run(aa_pipeline_api<unsigned short> *const obj);
    }
  }
}

#endif
