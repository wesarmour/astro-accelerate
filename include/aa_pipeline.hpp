#ifndef ASTRO_ACCELERATE_AA_PIPELINE_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_HPP

#include <iostream>
#include <stdio.h>
#include <memory>
#include <map>

#include "aa_compute.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_permitted_pipelines_2.hpp"
#include "aa_permitted_pipelines_3.hpp"
#include "aa_permitted_pipelines_4.hpp"
#include "aa_permitted_pipelines_5.hpp"

namespace astroaccelerate {

  /**
   * \class aa_pipeline aa_pipeline.hpp "include/aa_pipeline.hpp"
   * \brief Class to manage pipelines and their constituent modules, plans, and strategies.
   * \brief This class also delegates the movement of host memory into and out of the pipeline, 
   * \brief This class also manages the movement of device memory into and out of the device.
   * \details The class is templated over the input and output data type.
   * \details The class receives plan objects and calculates strategies.
   * \details The user may obtain the strategies.
   * \details The pipeline will not run unless all plans and strategies are successfully calculates for the pipeline that the user provided at construction.
   * \todo If any plan is offered to a bind method after a bind has already happened, then all flags indicating readiness must be invalidated.
   * \todo Nice to have but not needed: Add a way to transfer ownership of the data between aa_pipeline objects.
   * \todo Multiple bind calls which result in valid strategy objects should result in removal of the existing strategy so that only one strategy of each type exists.
   * \author Cees Carels.
   * \date: 23 October 2018.
   */

  template<typename T, typename U>
  class aa_pipeline {
  public:

    /** \brief Constructor for aa_pipeline that takes key parameters required on construction. */
    aa_pipeline(const aa_compute::pipeline &requested_pipeline,
		const aa_compute::pipeline_detail &pipeline_details,
		const aa_filterbank_metadata &filterbank_metadata,
		T const*const input_data,
		const aa_device_info::aa_card_info &card_info) : m_pipeline_details(pipeline_details),
								 m_card_info(card_info),
								 m_filterbank_metadata(filterbank_metadata),
								 bound_with_raw_ptr(true),
								 input_data_bound(true),
								 data_on_device(false),
								 pipeline_ready(false),
								 ptr_data_in(input_data) {
      
      //Add requested pipeline modules
      for(auto i : requested_pipeline) {
	required_plans.insert(std::pair<aa_compute::modules, bool>(i, true));
	supplied_plans.insert(std::pair<aa_compute::modules, bool>(i, false));
      }
      m_all_strategy.reserve(requested_pipeline.size());
      m_requested_pipeline = requested_pipeline;
    }
    
    /** \brief Destructor for aa_pipeline, performs cleanup as needed. */
    ~aa_pipeline() {
      pipeline_ready = false;
      if(data_on_device) {
	//There is still data on the device and the user has forgotten about it.
	//The device memory should be freed.
	std::cout << "ERROR:  Data still on device." << std::endl;
      }
        
      if(!bound_with_raw_ptr && input_data_bound) {
	//There is still (host) input data bound to the pipeline and the user has forgotten about it.
	//This (host) memory should be freed.
	std::cout << "ERROR:  Host data still bound to pipeline." << std::endl;
      }
    }
    
    /**
     * \brief Bind an aa_ddtr_plan to the aa_pipeline instance.
     * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
     */
    bool bind(aa_ddtr_plan plan) {
      pipeline_ready = false;
        
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::dedispersion) != required_plans.end()) {
	m_ddtr_plan = plan;
            
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::dedispersion) = true;
            
	//ddtr_strategy needs to know if analysis will be required
	if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
	  aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, true);
	  if(ddtr_strategy.ready()) {
	    m_ddtr_strategy = std::move(ddtr_strategy);
	    m_all_strategy.push_back(&m_ddtr_strategy);
	  }
	  else {
	    return false;
	  }
	}
	else {
	  aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, false);
	  if(ddtr_strategy.ready()) {
	    m_ddtr_strategy = std::move(ddtr_strategy);
	    m_all_strategy.push_back(&m_ddtr_strategy);
	  }
	  else {
	    std::cout << "ddtr_strategy not ready" << std::endl;
	    return false;
	  }
	}
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }
    

    /**
     * \brief Bind an aa_analysis_plan to the aa_pipeline instance.
     * \details If your pipeline includes analysis, then you must bind a valid aa_analysis_plan.
     * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
     */
    bool bind(aa_analysis_plan plan) {
      pipeline_ready = false;
      
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
	//Is the ddtr_strategy provided by this analysis_plan ready?
	if(!plan.ddtr_strategy().ready()) {
	  //This ddtr_strategy is not ready, so ignore this analysis_plan.
	  return false;
	}
	
	m_analysis_plan = plan;

	aa_analysis_strategy analysis_strategy(m_analysis_plan);
	if(analysis_strategy.ready()) {
	  m_analysis_strategy = std::move(analysis_strategy);
	  m_all_strategy.push_back(&m_analysis_strategy);
	}
	else {
	  return false;
	}
	
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::analysis) = true;
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }

    /**
     * \brief Bind an aa_periodicity_plan to the aa_pipeline instance.
     * \details If your pipeline includes periodicity, then you must bind a valid aa_periodicity_plan.
     * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
     */
    bool bind(aa_periodicity_plan plan) {
      pipeline_ready = false;
        
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::periodicity) != required_plans.end()) {
	m_periodicity_plan = plan;
	aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
	if(periodicity_strategy.ready()) {
	  m_periodicity_strategy = std::move(periodicity_strategy);
	  m_all_strategy.push_back(&m_periodicity_strategy);
	}
	else {
	  return false;
	}
	  
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::periodicity) = true;
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }
    
    /**
     * \brief Bind an aa_fdas_plan to the aa_pipeline instance.
     * \details If your pipeline includes fdas, then you must bind a valid aa_fdas_plan.
     * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
     */
    bool bind(aa_fdas_plan plan) {
      pipeline_ready = false;
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::fdas) != required_plans.end()) {
        m_fdas_plan = plan;
        aa_fdas_strategy fdas_strategy(m_fdas_plan);
        if(fdas_strategy.ready()) {
          m_fdas_strategy = std::move(fdas_strategy);
          m_all_strategy.push_back(&m_fdas_strategy);
        }
        else {
          return false;
        }
	
        //If the plan is valid then the supplied_plan becomes true 
        supplied_plans.at(aa_compute::modules::fdas) = true;
      }
      else {
        //The plan is not required, ignore.
        return false;
      }
      
      return true;
    }
    
    /** \returns The aa_ddtr_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_ddtr_strategy does not yet exist. */
    aa_ddtr_strategy ddtr_strategy() {
      //Does the pipeline actually need this strategy?
      if(required_plans.find(aa_compute::modules::dedispersion) != required_plans.end()) {
	//It does need this strategy.
	//Is it already computed?
	if(m_ddtr_strategy.ready()) {
	  //Return since it was already computed.
	  return m_ddtr_strategy;
	}
	else {
	  //ddtr_strategy was not yet computed, do it now.
	  //ddtr_strategy needs to know if analysis will be required
	  if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) { //analysis will be required
	    aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, true);
	    if(ddtr_strategy.ready()) {
	      m_ddtr_strategy = std::move(ddtr_strategy);
	      m_all_strategy.push_back(&m_ddtr_strategy);
	    }
	    else { //Tried to calculate ddtr_strategy with analysis enabled, but failed.
	      aa_ddtr_strategy empty_strategy;
	      return empty_strategy;
	    }
	  }
	  else { //analysis will not be required
	    aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, false);
	    if(ddtr_strategy.ready()) {
	      m_ddtr_strategy = std::move(ddtr_strategy);
	      m_all_strategy.push_back(&m_ddtr_strategy);
	    }
	    else { //Tried to calculate ddtr_strategy with analysis disabled, but failed.
	      std::cout << "ddtr_strategy not ready" << std::endl;
	      aa_ddtr_strategy empty_strategy;
	      return empty_strategy;
	    }
	  }
	}
	
	return m_ddtr_strategy;
      }
      else {
	//The pipeline does not need this strategy
	aa_ddtr_strategy empty_strategy;
	return empty_strategy;
      }
    }

    /** \returns The aa_analysis_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_analysis_strategy does not yet exist. */
    aa_analysis_strategy analysis_strategy() {
      //Does the pipeline actually need this strategy? 
      if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
        //It does need this strategy.                                                                                                                                                                                        
        //Is it already computed?
        if(m_analysis_strategy.ready()) { //Return since it was already computed.
          return m_analysis_strategy;
        }
	else {
	  //analysis_strategy was not yet computed, do it now.
	  aa_analysis_strategy analysis_strategy(m_analysis_plan);
	  if(analysis_strategy.ready()) {
	    m_analysis_strategy = std::move(analysis_strategy);
	    m_all_strategy.push_back(&m_analysis_strategy);
	  }
	  else { //Tried to calculate analysis_strategy, but failed.
	    aa_analysis_strategy empty_strategy;
	    return empty_strategy;
	  }
	}
      }
      else {
	//The pipeline does not need this strategy
	aa_analysis_strategy empty_strategy;
	return empty_strategy;
      }
    }

    /** \returns The aa_periodicity_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_periodicity_strategy does not yet exist. */
    aa_periodicity_strategy periodicity_strategy() {
      //Does the pipeline actually need this strategy?
      if(required_plans.find(aa_compute::modules::periodicity) != required_plans.end()) {
	//It does need this strategy.
	//Is it already computed?
	if(m_periodicity_strategy.ready()) { //Return since it was already computed.
	  return m_periodicity_strategy;
	}
	else {
	  //periodicity_strategy was not yet computed, do it now.
	  aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
	  if(periodicity_strategy.ready()) {
	    m_periodicity_strategy = std::move(periodicity_strategy);
	    m_all_strategy.push_back(&m_periodicity_strategy);
	  }
	  else { //Tried to calculate periodicity strategy, but failed.
	    aa_periodicity_strategy empty_strategy;
	    return empty_strategy;
	  }
	}
      }
      else {
	//The pipeline does not need this strategy
	aa_periodicity_strategy empty_strategy;
	return empty_strategy;
	
      }
    }

    /** \returns The aa_fdas_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_fdas_strategy does not yet exist. */
    aa_fdas_strategy fdas_strategy() {
      //Does the pipeline actually need this strategy?
      if(required_plans.find(aa_compute::modules::fdas) != required_plans.end()) {
        //It does need this strategy.
        //Is it already computed?
        if(m_fdas_strategy.ready()) { //Return since it was already computed.                                                                                                                                                          
          return m_fdas_strategy;
        }
        else {
          //fdas_strategy was not yet computed, do it now.
          aa_fdas_strategy fdas_strategy(m_fdas_plan);
          if(periodicity_strategy.ready()) {
            m_fdas_strategy = std::move(fdas_strategy);
            m_all_strategy.push_back(&m_fdas_strategy);
          }
          else { //Tried to calculate fdas strategy, but failed.
            aa_fdas_strategy empty_strategy;
            return empty_strategy;
          }
        }
      }
      else {
        //The pipeline does not need this strategy
        aa_fdas_strategy empty_strategy;
        return empty_strategy;
      }
    }
    
    /** \deprecated This functionality is delegated to the permitted_pipelines implementations. */
    bool transfer_data_to_device() {
      pipeline_ready = false;
      if(input_data_bound) {
	data_on_device = true;
      }
        
      return true;
    }
    
    /** \deprecated This functionality is delegated to the permitted_pipelines implementations. */
    bool transfer_data_to_host(std::vector<U> &data) {
      if(data_on_device) {
	data = std::move(data_out);
	data_on_device = false;
	pipeline_ready = false;
	return true;
      }
        
      return false;
    }
    
    /** \deprecated This functionality is delegated to the permitted_pipelines implementations. */
    bool transfer_data_to_host(U *&data) {
      if(data_on_device) {
	data = ptr_data_out;
	data_on_device = false;
	pipeline_ready = false;
	return true;
      }
        
      return false;
    }
    
    /** \deprecated This functionality is delegated to the permitted_pipelines implementations. */
    bool unbind_data() {
      /**
       * If the data is managed, then either it was either moved out via the transfer,
       * or it will be freed at de-allocation if the user forgot to transfer.
       *
       * If the data is unmanaged, then the data unbind call does not apply.
       * If the data came via a raw pointer, then the unbind call does not apply.
       *
       */
      pipeline_ready = false;
      return true;
    }
    
    /**
     * \brief Check whether all modules and strategies are ready for the pipeline to be able to execute.
     * \returns A boolean to indicate whether the pipeline is ready (true) or not (false).
     * \details If false, the user should check whether all plans have been bound and whether all strategies resulting from those plans are valid.
     */
    bool ready() {
      pipeline_ready = false;
        
      if(!data_on_device) {
	std::cout << "ERROR:  Data is not on device or properly bound to API." << std::endl;
	return false;
      }

      // Check if all plans are supplied.
      // If not, then one or more strategies will not be ready.
      // If so, return false.
      for(auto const& i : supplied_plans) {
	if(i.second == false) {
	  std::cout << "ERROR:  " << aa_compute::module_name(i.first) << " plan is not ok." << std::endl;
	  return false;
	}
      }

      // Check if all strategies are ready.
      // If not, then return false.
      bool all_strategies_ready = true;
      for(auto strategy : m_all_strategy) {
	if(strategy->setup()) {
	  //Memory allocations and setup happened successfully
	}
	else {
	  //Setup failed.
	  std::cout << "ERROR:  " << strategy->name() << " setup failed." << std::endl;
	  all_strategies_ready = false;
	}
      }

      if(!all_strategies_ready) {
	return false;
      }

      // Start configuring the pipeline that will be run.
      constexpr aa_compute::module_option zero_dm               = aa_compute::module_option::zero_dm;
      constexpr aa_compute::module_option zero_dm_with_outliers = aa_compute::module_option::zero_dm_with_outliers;

      if(m_pipeline_details.find(zero_dm) == m_pipeline_details.end() &&
	 m_pipeline_details.find(zero_dm_with_outliers) == m_pipeline_details.end()) {
	std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were selected. Selection OFF." << std::endl;
      }
      
      if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end() &&
	 m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	std::cout << "ERROR:  Both zero_dm and zero_dm_with_outliers cannot simultaneously be selected." << std::endl;
	return false;
      }

      //If fdas acceleration is used, these flags must be set
      bool fdas_enable_custom_fft        = false;
      bool fdas_enable_inbin             = false;
      bool fdas_enable_norm              = false;
      bool fdas_enable_output_ffdot_plan = false;
      bool fdas_enable_output_list       = false;

      if(m_pipeline_details.find(aa_compute::module_option::fdas_custom_fft) != m_pipeline_details.end()) {
	fdas_enable_custom_fft = true;
      }
	
      if(m_pipeline_details.find(aa_compute::module_option::fdas_inbin) != m_pipeline_details.end()) {
	fdas_enable_inbin = true;
      }

      if(m_pipeline_details.find(aa_compute::module_option::fdas_norm) != m_pipeline_details.end()) {
	fdas_enable_norm = true;
      }
	
      if(m_pipeline_details.find(aa_compute::module_option::output_ffdot_plan) != m_pipeline_details.end()) {
	fdas_enable_output_ffdot_plan = true;
      }
	
      if(m_pipeline_details.find(aa_compute::module_option::output_fdas_list) != m_pipeline_details.end()) {
	fdas_enable_output_list = true;
      }
	
      /**
       * The following code requires some explanation, especially for developers unfamiliar with modern C++.
       * The purpose of the following block of code is to be a runtime dispatch. What this does is look up
       * which version of the pipeline function matches what the user requested.
       *
       * All pipelines are contained inside classes of the form aa_permitted_pipelines_x, where x is a number.
       * These pipeline classes are templated over two parameters, a module_option relating to which zero_dm
       * version will be used, and another to indicate which candidate algorithm will be used.
       * The zero_dm option is contained inside aa_compute::module_option::zero_dm
       * and aa_compute::module_option::zero_dm_with_outliers.
       * To make the code less verbose, constexpr are used to make create "abbreviations" for this syntax.
       *
       * All aa_permitted_pipeline_x classes are derived from a base class aa_pipeline_runner.
       * A std::unique_ptr is used because this provides automatic storage (when the std::unique_ptr goes out of
       * scope, the destructor is called for the object to which it points).
       * A base class pointer of type aa_pipeline_runner is used so that it can point to any of the aa_permitted_pipelines_x.
       * Runtime polymorphism ensures the methods for the selected aa_permitted_pipeline_x class are called.
       *
       * The code checks the m_requested_pipeline against the possible pipelines, and selects the one matching the user
       * settings. The syntax of std::unique_ptr implies that the templated class and the template parameters must be specified
       * both when specifying the type of the std::unique_ptr, and when the type is "new"'d.
       *
       * Lastly, because the type is being constructed, the parameters for the constructor must also be provided.
       *
       */

      constexpr aa_compute::module_option off                   = aa_compute::module_option::empty;
      constexpr aa_compute::module_option old_rfi               = aa_compute::module_option::old_rfi;
      constexpr bool use_old_rfi = true;
      constexpr bool use_new_rfi = false;
      //Check which pipeline the user has requested (given by m_requested_pipeline) against the possible permitted pipelines.
      //Then, assign a new object of that type to the base class pointer.
      bool is_pipeline_set_to_runner = false;
      if(m_requested_pipeline == aa_permitted_pipelines::pipeline1) {
	if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end()) {
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm,               use_old_rfi>>(new aa_permitted_pipelines_1<zero_dm,               use_old_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 1 with zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm and do not contain old_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm,               use_new_rfi>>(new aa_permitted_pipelines_1<zero_dm,               use_new_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 1 with zero_dm, and new_rfi" << std::endl;
	  }
	}
	else if(m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	  //details contain zero_dm_with_outliers
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm_with_outliers and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_1<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 1 with zero_dm_with_outliers, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm_with_outliers, use_new_rfi>>(new aa_permitted_pipelines_1<zero_dm_with_outliers, use_new_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 1 with zero_dm_with_outliers, and new_rfi" << std::endl;
	  }
	}
	else {
	  std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF." << std::endl;
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<off, use_old_rfi>>(new aa_permitted_pipelines_1<off, use_old_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 1 without zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_1<off, use_new_rfi>>(new aa_permitted_pipelines_1<off, use_new_rfi>(m_ddtr_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 1 without zero_dm, and new_rfi" << std::endl;
	  }
	}	
      }
      else if(m_requested_pipeline == aa_permitted_pipelines::pipeline2) {
	if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end()) {
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm,               use_old_rfi>>(new aa_permitted_pipelines_2<zero_dm,               use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 with zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm and do not contain old_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm,               use_new_rfi>>(new aa_permitted_pipelines_2<zero_dm,               use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 with zero_dm, and new_rfi" << std::endl;
	  }
	}
	else if(m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	  //details contain zero_dm_with_outliers
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm_with_outliers and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_2<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 with zero_dm_with_outliers, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm_with_outliers, use_new_rfi>>(new aa_permitted_pipelines_2<zero_dm_with_outliers, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 with zero_dm_with_outliers, and new_rfi" << std::endl;
	  }
	}
	else {
	  std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF." << std::endl;
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<off, use_old_rfi>>(new aa_permitted_pipelines_2<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 without zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_2<off, use_new_rfi>>(new aa_permitted_pipelines_2<off, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 2 without zero_dm, and new_rfi" << std::endl;
	  }
	}
      }
      else if(m_requested_pipeline == aa_permitted_pipelines::pipeline3) {
	if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end()) {
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm,               use_old_rfi>>(new aa_permitted_pipelines_3<zero_dm,               use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 with zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm and do not contain old_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm,               use_new_rfi>>(new aa_permitted_pipelines_3<zero_dm,               use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 with zero_dm, and new_rfi" << std::endl;
	  }
	}
	else if(m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	  //details contain zero_dm_with_outliers
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm_with_outliers and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_3<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 with zero_dm_with_outliers, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm_with_outliers, use_new_rfi>>(new aa_permitted_pipelines_3<zero_dm_with_outliers, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 with zero_dm_with_outliers, and new_rfi" << std::endl;
	  }
	}
	else {
	  std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF." << std::endl;
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<off, use_new_rfi>>(new aa_permitted_pipelines_3<off, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 without zero_dm, and new_rfi" << std::endl;
	  }
	  else {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_3<off, use_old_rfi>>(new aa_permitted_pipelines_3<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
	    is_pipeline_set_to_runner	= true;
	    std::cout << "NOTICE: Selected Pipeline 3 without zero_dm, and old_rfi" << std::endl;
	  }
	}
      }
      else if(m_requested_pipeline == aa_permitted_pipelines::pipeline4) {
	if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end()) {
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm,               use_old_rfi>>(new aa_permitted_pipelines_4<zero_dm,               use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 with zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm and do not contain old_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm,               use_new_rfi>>(new aa_permitted_pipelines_4<zero_dm,               use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 with zero_dm, and new_rfi" << std::endl;
	  }
	}
	else if(m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	  //details contain zero_dm_with_outliers
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm_with_outliers and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_4<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 with zero_dm_with_outliers, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm_with_outliers, use_new_rfi>>(new aa_permitted_pipelines_4<zero_dm_with_outliers, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 with zero_dm_with_outliers, and new_rfi" << std::endl;
	  }
	}
	else {
	  std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF." << std::endl;
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<off, use_old_rfi>>(new aa_permitted_pipelines_4<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 without zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_4<off, use_new_rfi>>(new aa_permitted_pipelines_4<off, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 4 without zero_dm, and new_rfi" << std::endl;
	  }
	}
      }
      else if(m_requested_pipeline == aa_permitted_pipelines::pipeline5) {
	if(m_pipeline_details.find(zero_dm) != m_pipeline_details.end()) {
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm,               use_old_rfi>>(new aa_permitted_pipelines_5<zero_dm,               use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 with zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm and do not contain old_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm,               use_new_rfi>>(new aa_permitted_pipelines_5<zero_dm,               use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 with zero_dm, and new_rfi" << std::endl;
	  }
	}
	else if(m_pipeline_details.find(zero_dm_with_outliers) != m_pipeline_details.end()) {
	  //details contain zero_dm_with_outliers
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    //details contain zero_dm_with_outliers and old_rfi
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_5<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 with zero_dm_with_outliers, and old_rfi" << std::endl;
	  }
	  else {
	    //details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm_with_outliers, use_new_rfi>>(new aa_permitted_pipelines_5<zero_dm_with_outliers, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 with zero_dm_with_outliers, and new_rfi" << std::endl;
	  }
	}
	else {
	  std::cout << "NOTICE: Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF." << std::endl;
	  if(m_pipeline_details.find(old_rfi) != m_pipeline_details.end()) {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<off, use_old_rfi>>(new aa_permitted_pipelines_5<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 without zero_dm, and old_rfi" << std::endl;
	  }
	  else {
	    m_runner = std::unique_ptr<aa_permitted_pipelines_5<off, use_new_rfi>>(new aa_permitted_pipelines_5<off, use_new_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
	    is_pipeline_set_to_runner = true;
	    std::cout << "NOTICE: Selected Pipeline 5 without zero_dm, and new_rfi" << std::endl;
	  }
	}
      }
      else {
	//Pipeline 0
      }
	
      //Do any last checks on the plans as a whole
      if(is_pipeline_set_to_runner) {
	pipeline_ready = true;
      }

      std::cout << "---PIPELINE DIAGNOSTIC INFORMATION---" << std::endl;
      aa_device_info::print_card_info(m_card_info);

      if(required_plans.find(aa_compute::modules::dedispersion) != required_plans.end()) {
	aa_ddtr_strategy::print_info(m_ddtr_strategy);
      }

      if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
	aa_analysis_strategy::print_info(m_analysis_strategy);
      }

      if(required_plans.find(aa_compute::modules::periodicity) != required_plans.end()) {
	aa_periodicity_strategy::print_info(m_periodicity_strategy);
      }

      if(required_plans.find(aa_compute::modules::fdas) != required_plans.end()) {
	aa_fdas_strategy::print_info(m_fdas_strategy);
      }
      
      pipeline_ready = true;
      return true;
    }
    

    /** \brief Runs the pipeline end-to-end. */
    bool run() {
      /**
       * This method to be overloaded with all possible combinations of
       * data that the user may wish to extract from any pipeline.
       * Any methods that are not supported are compile-time errors because
       * the base class must provide a method for it.
       */
      if(pipeline_ready && m_runner->setup()) {
	while(m_runner->next()) {
	  std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
	}
	return true;
      }
      else {
	return false;
      }
    }

    /**
     * \brief Function pass input/output data from one aa_pipeline instance to another.
     * \warning Not yet implemented.
     */
    bool handoff(aa_pipeline &next_pipeline) {
      /**
       * Handoff control over the data to the next pipeline.
       *
       * The user must always retrieve the data from the last pipeline in the chain.
       */
      return true;
    }
    
  private:
    std::map<aa_compute::modules, bool> required_plans; /** Plans required to configure the pipeline. */
    std::map<aa_compute::modules, bool> supplied_plans; /** Plans supplied by the user. */
    std::vector<aa_strategy*>   m_all_strategy; /** Base class pointers to all strategies bound to the pipeline. */
    aa_compute::pipeline        m_requested_pipeline; /** The user requested pipeline that was bound to the aa_pipeline instance on construction. */
    const aa_compute::pipeline_detail &m_pipeline_details; /** The user requested pipeline details containing module options for the aa_pipeline instance. */
    aa_device_info::aa_card_info m_card_info; /** The user provided GPU card information for the aa_pipeline instance. */
    std::unique_ptr<aa_pipeline_runner> m_runner; /** A std::unique_ptr that will point to the correct class instantation of the selected aa_permitted_pipelines_ when the pipeline must be made ready to run. */
    
    aa_filterbank_metadata      m_filterbank_metadata; /** The filterbank file metadata that the user provided for the aa_pipeline instance on construction. */
    
    aa_ddtr_plan                m_ddtr_plan; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    aa_ddtr_strategy            m_ddtr_strategy; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    
    aa_analysis_plan            m_analysis_plan; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    aa_analysis_strategy        m_analysis_strategy; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    
    aa_periodicity_plan         m_periodicity_plan; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    aa_periodicity_strategy     m_periodicity_strategy; /** The instance of this type that is currently bound to the aa_pipeline instance. */

    aa_fdas_plan                m_fdas_plan; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    aa_fdas_strategy            m_fdas_strategy; /** The instance of this type that is currently bound to the aa_pipeline instance. */
    
    bool bound_with_raw_ptr; /** Flag to indicate whether the input data is bound via a raw pointer (true) or not (false). */
    bool input_data_bound; /** Flag to indicate whether the input data is bound already (true) or not (false). */
    bool data_on_device; /** \deprecated Flag to indicate whether the input data is on the device (true) or not (false). */
    bool pipeline_ready; /** Flag to indicate whether the pipeline is ready to execute (true) or not (false).  */
    
    std::vector<T>              data_in; /** Input data buffer. */
    std::vector<U>              data_out; /** Output data buffer. */
    T const*const               ptr_data_in; /** Input data pointer if bound_with_raw_ptr is true. */
    U*                          ptr_data_out; /** Input data pointer if bound_with_raw_ptr is true. */
  };

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_PIPELINE_HPP
