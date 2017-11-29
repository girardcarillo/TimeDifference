// - Implementation of TimeDifference
// Ourselves
#include "TimeDifference.h"
// Standard Library
// Third Party
// - A
// This Project
#include <snemo/datamodels/data_model.h>
#include <snemo/datamodels/topology_data.h>
#include <mctools/simulated_data.h>

#include <snemo/cuts/energy_measurement_cut.h>

//using namespace std

// Macro which automatically implements the interface needed
// to enable the module to be loaded at runtime
// The first argument is the typename
// The second is the string key used to access the module in pipeline
// scripts. This must be globally unique.
DPP_MODULE_REGISTRATION_IMPLEMENT(TimeDifference,"TimeDifference")

// Construct
TimeDifference::TimeDifference() : dpp::base_module()
{
  this->_set_defaults();
}

// Destruct
TimeDifference::~TimeDifference() {
  // MUST reset module at destruction
  this->reset();
}

// Initialize
void TimeDifference::initialize(const datatools::properties& setup_,
                          datatools::service_manager& /*flServices*/,
                          dpp::module_handle_dict_type& /*moduleDict*/) {
  dpp::base_module::_common_initialize(setup_);
  this->_set_initialized(true);
}

// Process
dpp::base_module::process_status
TimeDifference::process(datatools::things& data_record_) {
  DT_THROW_IF(! is_initialized(), std::logic_error,
              "Module '" << get_name () << "' is not initialized !");
  data_record_.tree_dump();

  ////Defining labels
  //td_label
  //const std::string & td_label = snemo::datamodel::data_info::default_topology_data_label();

  //sd_label
  const std::string & sd_label = snemo::datamodel::data_info::default_simulated_data_label();

  ////Storing data bases
  //topology data base
  if (data_record_.has("TD")) {
    const snemo::datamodel::topology_data & a_td
      =  data_record_.get<snemo::datamodel::topology_data>("TD");
    a_td.tree_dump();
  }
  else {
    DT_LOG_DEBUG(get_logging_priority(),"Data has no TD");
  }

  //simulated data base
  if (data_record_.has(sd_label)) {
    const mctools::simulated_data & a_sd
      = data_record_.get<mctools::simulated_data>(sd_label);
    a_sd.tree_dump();
  }
  else {
    DT_LOG_DEBUG(get_logging_priority(),"Data has no SD");
  }

  ////Applying cuts on data bases
  ///Cut on TD base

  //Energy cut



  //Topology cut

  ///Cut on SD base

  //Topology cut

  // MUST return a status, see ref dpp::base_module::process_status
  return PROCESS_OK;
}

// Reset
void TimeDifference::reset() {
  this->_set_initialized(false);
}

void TimeDifference::_set_defaults() {
}
