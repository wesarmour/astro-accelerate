#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"
#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_host_fake_signal_generator.hpp"

#include "aa_permitted_pipelines_3.hpp"

#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"

#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"

using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 120, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).

  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double total_bandwidth = 300.0f;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  const unsigned int nsamples = 1.5/tsamp; // 2s of data in current tsamp at minimum; must be more than signal_start + maxshift
  const double fch1 = 1550;
  const int nchans = 128;
  const double foff = -total_bandwidth/nchans;

  // setting the signal metadata
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);

  // Init the GPU card
  aa_device_info& device_info = aa_device_info::instance();
  if(device_info.check_for_devices()) {
    LOG(log_level::notice, "Checked for devices.");
  }
  else {
    LOG(log_level::error, "Could not find any devices.");
  }

  aa_device_info::CARD_ID selected_card = 0;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info.init_card(selected_card, selected_card_info)) {
    LOG(log_level::notice, "init_card complete. Selected card " + std::to_string(selected_card) + ".");
  }
  else {
    LOG(log_level::error, "init_card incomplete.")
  }

  aa_device_info::print_card_info(selected_card_info);
  
  const size_t free_memory = selected_card_info.free_memory; // Free memory on the GPU in bytes
  bool enable_analysis = true; 

  aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis);

  if(!(strategy.ready())) {
    std::cout << "There was an error" << std::endl;
    return 0;
  }

  const float sigma_cutoff = 1.0;
  const float sigma_constant = 4.0;
  const float max_boxcar_width_in_sec = 0.05;
  const aa_analysis_plan::selectable_candidate_algorithm algo = aa_analysis_plan::selectable_candidate_algorithm::off;

  aa_analysis_plan analysis_plan(strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, algo, false);
  aa_analysis_strategy analysis_strategy(analysis_plan);

  if(!(analysis_strategy.ready())) {
    std::cout << "ERROR: analysis_strategy not ready." << std::endl;
    return 0;
  }

  const float periodicity_sigma_cutoff = 0.0;
  const float periodicity_sigma_constant = sigma_constant;
  const int   nHarmonics = 3;
  const int   export_powers = 0;
  const bool  candidate_algorithm = true;
  const bool  enable_outlier_rejection = false;

  aa_periodicity_plan periodicity_plan(periodicity_sigma_cutoff, periodicity_sigma_constant, nHarmonics, export_powers, candidate_algorithm, enable_outlier_rejection);
  aa_periodicity_strategy periodicity_strategy(periodicity_plan);

  if(!periodicity_strategy.ready()) {
    std::cout << "ERROR: periodicity_strategy not ready." << std::endl;
  }

  // params needed by the fake signal function
  double dm_position = 90; // at what dm put the signal
  const int func_width = 1/(tsamp*10); // width of the signal in terms of # of samples; now at 1% of samling rate
  const int period = 200; // pulsar period with ms
  bool dump_to_disk = true;
  const float sigma = 0.5;

  // setting the metadata for running fake generator
  aa_fake_signal_metadata f_meta(dm_position, func_width, sigma, period);

  aa_fake_signal_generator signal;
  signal.create_signal(strategy, f_meta, metadata, dump_to_disk);

  if(!(signal.ready())) {
	std::cout << "Error in creating fake signal" << std::endl;
	return 0;
  }

  std::vector<unsigned short> input_data;
  input_data = signal.signal_data();

  aa_permitted_pipelines_3<aa_pipeline::component_option::empty, false> runner(strategy, analysis_strategy, periodicity_strategy, input_data.data());
  
  if(runner.setup()) {
    while(runner.next()) {
      LOG(log_level::notice, "Pipeline running over next chunk.");
    }
  }

  signal.print_info(f_meta);
  strategy.print_info(strategy);

  LOG(log_level::notice, "Finished.");

  return 0;
}