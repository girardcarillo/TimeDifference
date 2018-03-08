// - Implementation of TimeDifference
// Ourselves
#include "TimeDifference.h"
// Standard Library
// Third Party
// - A
// This Project
// #include <datatools/bit_mask.h>
#include <snemo/datamodels/data_model.h>
#include <snemo/datamodels/topology_data.h>
#include <snemo/datamodels/particle_track_data.h>
// #include <snemo/reconstruction/charged_particle_tracking_module.h>
#include <snemo/reconstruction/tof_driver.h>
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

  ////Root export
  // std::stringstream ss;
  // ss.str("");
  // ss << "sd_" << _output_file_name_;
  // ss.str();
  _sd_output_file_ = new TFile ("sd_tree.root", "recreate", "Output file of Simulation data");
  _sd_output_file_->cd();
  _sd_tree_= new TTree("calorimeter_hit", "calorimeter_hit");
  _sd_tree_->Branch("time", &_time_,"time/D");
  _sd_tree_->Branch("time_difference_E", &_time_difference_E_,"time_difference_E/D");
  _sd_tree_->Branch("probability", &_internal_probability_,"probability/D");
  _sd_tree_->Branch("length_Emin", &_length_Emin_,"length_Emin/D");
  _sd_tree_->Branch("length_Emax", &_length_Emax_,"length_Emax/D");
  // _sd_tree_->Branch("Energy", &_energy_,"energy/D");
  _sd_tree_->Branch("minimal_energy", &_minimal_energy_,"minimal_energy/D");
  _sd_tree_->Branch("maximal_energy", &_maximal_energy_,"maximal_energy/D");
  _sd_tree_->Branch("time_Emin", &_time_Emin_,"time_Emin/D");
  _sd_tree_->Branch("time_Emax", &_time_Emax_,"time_Emax/D");
  _sd_tree_->Branch("sigma_time_Emin", &_sigma_time_Emin_,"sigma_time_Emin/D");
  _sd_tree_->Branch("sigma_time_Emax", &_sigma_time_Emax_,"sigma_time_Emax/D");
  this->_set_initialized(true);

  // std::cout << "Enter output file name (ex.: flRec)" << std::endl;
  // std::cin >> _output_file_name_;
}

std::string const other_events("other_events.txt");
std::ofstream other_events_flux(other_events.c_str());

std::string const final_rate("final_rate.txt");
std::ofstream final_flux(final_rate.c_str());

// Process
dpp::base_module::process_status
TimeDifference::process(datatools::things& data_record_) {
  DT_THROW_IF(! is_initialized(), std::logic_error,
              "Module '" << get_name () << "' is not initialized !");

  //Counting the number of simulated events
  _number_event_++;
  // std::cout << "Number event = " << _number_event_ << std::endl;

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

  //Particle track data base
  // DT_THROW_IF(! data_record_.has("PTD"), std::logic_error,
  //             "Data has no PTD !");
  // const snemo::datamodel::particle_track_data & a_ptd
  //   =  data_record_.get<snemo::datamodel::particle_track_data>("PTD");

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
  double my_energy_sum = 0;
  double my_internal_probability = 0;
  double my_length_Emin = 0;
  double my_length_Emax = 0;
  double my_minimal_energy = 0;
  double my_maximal_energy = 0;
  std::string my_maximal_energy_name;
  std::string my_minimal_energy_name;
  double my_time_Emin = 0;
  double my_time_Emax = 0;
  double my_sigma_time_Emin = 0;
  double my_sigma_time_Emax = 0;
  // double electron_mass_energy = 0.511;
  // std::cout << my_energy_sum << std::endl;
  if (a_td.has_pattern()
      && a_td.has_pattern_as<snemo::datamodel::topology_2e_pattern>()) {
    const snemo::datamodel::topology_2e_pattern & a_2e_topology
      = a_td.get_pattern_as<snemo::datamodel::topology_2e_pattern>();
    const double & a_energy_sum
      = a_2e_topology.get_electrons_energy_sum();
    const double & a_minimal_energy
      = a_2e_topology.get_electron_minimal_energy();
    const std::string & a_minimal_energy_name
      = a_2e_topology.get_minimal_energy_electron_name();
    const std::string & a_maximal_energy_name
      = a_2e_topology.get_maximal_energy_electron_name();
    const double & a_maximal_energy
      = a_2e_topology.get_electron_maximal_energy();
    const double & a_internal_probability
      = a_2e_topology.get_electrons_internal_probability();

    double time_Emin;
    double time_Emax;
    double sigma_time_Emin;
    double sigma_time_Emax;
    datatools::invalidate(time_Emin);
    datatools::invalidate(time_Emax);
    datatools::invalidate(sigma_time_Emin);
    datatools::invalidate(sigma_time_Emax);
    if (a_td.get_pattern().get_particle_track(a_minimal_energy_name).has_associated_calorimeter_hits()
        && a_td.get_pattern().get_particle_track(a_maximal_energy_name).has_associated_calorimeter_hits()) {
      const snemo::datamodel::calibrated_calorimeter_hit::collection_type & a_calorimeter_min =
        a_td.get_pattern().get_particle_track(a_minimal_energy_name).get_associated_calorimeter_hits();
      const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit_min = a_calorimeter_min.front().get();
      const snemo::datamodel::calibrated_calorimeter_hit::collection_type & a_calorimeter_max =
        a_td.get_pattern().get_particle_track(a_maximal_energy_name).get_associated_calorimeter_hits();
      const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit_max = a_calorimeter_max.front().get();
      time_Emin = a_calo_hit_min.get_time();
      time_Emax = a_calo_hit_max.get_time();
      sigma_time_Emin = a_calo_hit_min.get_sigma_time();
      sigma_time_Emax = a_calo_hit_max.get_sigma_time();
    }
    else {
      DT_THROW_IF(true,std::logic_error,"Particle track is not associated to any calorimeter block !");
    }

    double length_Emin = datatools::invalid_real();
    double length_Emax = datatools::invalid_real();
    if (a_td.get_pattern().get_particle_track(a_2e_topology.get_minimal_energy_electron_name()).has_trajectory()
        && a_td.get_pattern().get_particle_track(a_2e_topology.get_maximal_energy_electron_name()).has_trajectory()) {
      const snemo::datamodel::tracker_trajectory & a_trajectory_min =
        a_td.get_pattern().get_particle_track(a_2e_topology.get_minimal_energy_electron_name()).get_trajectory();
      const snemo::datamodel::base_trajectory_pattern & a_track_pattern_min = a_trajectory_min.get_pattern();
      length_Emin = a_track_pattern_min.get_shape().get_length();
      const snemo::datamodel::tracker_trajectory & a_trajectory_max =
        a_td.get_pattern().get_particle_track(a_2e_topology.get_maximal_energy_electron_name()).get_trajectory();
      const snemo::datamodel::base_trajectory_pattern & a_track_pattern_max = a_trajectory_max.get_pattern();
      length_Emax = a_track_pattern_max.get_shape().get_length();
    }
    else {
      DT_THROW_IF(true,std::logic_error,"Electron of minimal energy has no attached trajectory !");
    }

    if (a_energy_sum/CLHEP::MeV >= 2.7 && a_energy_sum/CLHEP::MeV <= 3.2) {
      _nb_2e_topology_++;
      my_energy_sum = a_energy_sum;
      my_minimal_energy = a_minimal_energy;
      my_maximal_energy = a_maximal_energy;
      my_minimal_energy_name = a_minimal_energy_name;
      my_maximal_energy_name = a_maximal_energy_name;
      my_internal_probability = a_internal_probability;
      my_length_Emin = length_Emin;
      my_length_Emax = length_Emax;
      my_time_Emin = time_Emin;
      my_time_Emax = time_Emax;
      my_sigma_time_Emin = sigma_time_Emin;
      my_sigma_time_Emax = sigma_time_Emax;
      // std::cout << "Internal probability = " << my_internal_probability << std::endl;
      // std::cout << "2e topology = " << _nb_2e_topology_ << std::endl;
      // std::cout << "Cut on TD base OK" << std::endl;
    }
  }

  //Cut on SD base
  int nb_electron = 0;
  double a_time_difference = 0;
  for (const auto & iparticle : the_primary_particles) {
    if (iparticle.is_electron()) {
      nb_electron++;
      // std::cout << "Number of electron = " << nb_electron << std::endl;
    }

    const std::string & a_particle_label
      = iparticle.get_particle_label();
    _particle_label_ = a_particle_label;
    // std::cout << "Particle type = " << _particle_label_ << std::endl;
    const double a_time
      = iparticle.get_time();
    _particle_time_ = a_time;
    if (nb_electron == 2) {
      // std::cout << "Cut on SD base OK" << std::endl;
      a_time_difference = fabs(a_time_difference - a_time);
      // std::cout << "Time difference = " << a_time_difference/CLHEP::picosecond << " ps" << std::endl;
    }
  }

  // std::cout << "Energy sum = " << my_energy_sum << std::endl;
  ////Storing data
  if (my_energy_sum != 0 && my_internal_probability != 0) {//Guarantee we entered in the TD cut loop
    //Keep interesting events in a root tree
    if (nb_electron == 2) {
      _nb_internal_conversion_++;
      if (a_time_difference != 0 && my_time_Emin != 0 && my_time_Emax != 0) {// To remove if working with 0nubb simulations
        _sd_output_file_->cd();
        _time_= a_time_difference/CLHEP::picosecond;
        _internal_probability_ = my_internal_probability;
        _length_Emin_ = my_length_Emin;
        _length_Emax_ = my_length_Emax;
        _energy_ = my_energy_sum;
        _minimal_energy_ = my_minimal_energy;
        _maximal_energy_ = my_maximal_energy;
        _time_Emin_ = my_time_Emin;
        _time_Emax_ = my_time_Emax;
        _time_difference_E_ = fabs(my_time_Emax - my_time_Emin);
        _sigma_time_Emin_ = my_sigma_time_Emin;
        _sigma_time_Emax_ = my_sigma_time_Emax;
        _sd_tree_->Fill();
        // std::cout << "Internal probability = " << _internal_probability_ << std::endl;
        // std::cout << "Energy sum = " << my_energy_sum << std::endl;
        // std::cout << "Calo--minimal energy : " << my_minimal_energy_name << std::endl;
        // std::cout << "Calo--maximal energy : " << my_maximal_energy_name << std::endl;
        // std::cout << "Sigma time of e- of min energy = " << my_sigma_time_Emin << std::endl;
        // std::cout << "sigma time of e- of max energy = " << my_sigma_time_Emax << std::endl;
        // std::cout << "Time of e- of min energy = " << my_time_Emin << std::endl;
        // std::cout << "Time of e- of max energy = " << my_time_Emax << std::endl;
      }
    }

    //Keep other events in a .txt file
    if (nb_electron != 2) {
      _nb_other_process_++;
      DT_THROW_IF(! other_events_flux, std::logic_error,
                  "ERROR: cannot open the other_events.txt file!");
      other_events_flux << "Event # " << _number_event_-1 << std::endl;
      other_events_flux << "Particle type = " << _particle_label_ << std::endl;
      other_events_flux << "Event generation time = " << _particle_time_ << "\n" <<std::endl;
      // std::cout << "There is other events" << std::endl;
    }
  }

  //Final rates stored in a .txt file
  DT_THROW_IF(! final_flux, std::logic_error,
              "ERROR: cannot open the final_rate.txt file!");
  final_flux << "Event # " << _number_event_-1 << std::endl;
  final_flux << "2e topology = " << ((double)_nb_2e_topology_/(double)_number_event_)*100 << "% : " << std::endl;
  final_flux << "   -Internal conversion = " << ((double)_nb_internal_conversion_/(double)_number_event_)*100 << "%" << std::endl;
  final_flux << "   -Other processes = " << ((double)_nb_other_process_/(double)_number_event_)*100 << "%" << std::endl;
  // std::cout << "Final flux" << std::endl;
  // MUST return a status, see ref dpp::base_module::process_status
  return PROCESS_OK;
}

// Reset
void TimeDifference::reset() {
  // Root tree
  _sd_output_file_->cd();
  // _sd_tree_->Write();
  _sd_output_file_->Write();
  _sd_output_file_->Close();

  this->_set_initialized(false);
}

void TimeDifference::_set_defaults() {
}
