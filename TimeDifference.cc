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
#include <snemo/datamodels/topology_2e_pattern.h>

// #include <sstream>

// #include <snemo/cuts/energy_measurement_cut.h>

// using namespace std

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

  std::cout << "Enter output file name (ex.: flRec)" << std::endl;
  std::cin >> _output_file_name_;
}


// Process
dpp::base_module::process_status
TimeDifference::process(datatools::things& data_record_) {
  DT_THROW_IF(! is_initialized(), std::logic_error,
              "Module '" << get_name () << "' is not initialized !");

  ////Root export
  // std::stringstream ss;
  // ss.str("");
  // ss << "sd_" << _output_file_name_;
  // ss.str();
  _sd_output_file_ = new TFile ("", "recreate", "Output file of Simulation data");
  _sd_output_file_->cd();
  _sd_tree_= new TTree("calorimeter hit tree", "calorimeter hit tree");
  _sd_tree_->Branch("time", &_time_,"time/D");

  ////Defining data base labels
  //td_label
  //const std::string & td_label = snemo::datamodel::data_info::default_topology_data_label();

  //sd_label
  const std::string & sd_label = snemo::datamodel::data_info::default_simulated_data_label();

  ////Storing data bases

  //Topology data base
  DT_THROW_IF(! data_record_.has("TD"), std::logic_error,
              "Data has no TD !");
  const snemo::datamodel::topology_data & a_td
    =  data_record_.get<snemo::datamodel::topology_data>("TD");

  //Simulated data base
  DT_THROW_IF(! data_record_.has(sd_label), std::logic_error,
              "Data has no SD !");
  const mctools::simulated_data & a_sd
    = data_record_.get<mctools::simulated_data>(sd_label);
  const mctools::simulated_data::primary_event_type & a_primary_event
    = a_sd.get_primary_event();
  const genbb::primary_event::particles_col_type & the_primary_particles
    = a_primary_event.get_particles();

  ////Applying cuts on data bases
  //Cut on TD base
    if (a_td.has_pattern() && a_td.has_pattern_as<snemo::datamodel::topology_2e_pattern>()) {
      const snemo::datamodel::topology_2e_pattern & a_2e_topology
        = a_td.get_pattern_as<snemo::datamodel::topology_2e_pattern>();
      const double & a_energy_sum
        = a_2e_topology.get_electrons_energy_sum();
      if (a_energy_sum/CLHEP::MeV >= 2.7 && a_energy_sum/CLHEP::MeV <= 3.2) {
        std::cout << "Cut on TD base OK" << std::endl;
      }
    }

    //Cut on SD base
    int count = 0;
    for (const auto & iparticle : the_primary_particles) {
      if (iparticle.is_electron()) {
        count++;
        // std::cout << "count = " << count << std::endl;
      }
      if (count == 2) {
        std::cout << "Cut on SD base OK" << std::endl;
      }
    }

    ////Storing data

  // MUST return a status, see ref dpp::base_module::process_status
  return PROCESS_OK;
}

// Reset
void TimeDifference::reset() {
  this->_set_initialized(false);
}

void TimeDifference::_set_defaults() {
}
