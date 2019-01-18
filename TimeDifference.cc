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
#include <bayeux/datatools/event_id.h>
#include <falaise/snemo/datamodels/vertex_measurement.h>
#include <falaise/snemo/datamodels/pid_utils.h>

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
  _sd_output_file_ = new TFile ("sd_tree.root", "recreate", "Output file of Simulation data");
  _sd_output_file_->cd();
  _sd_tree_= new TTree("calorimeter_hit", "calorimeter_hit");
  _sd_tree_->Branch("gen_time", &_gen_time_,"gen_time/D");
  _sd_tree_->Branch("event_counter", &_event_counter_,"event_counter/I");
  _sd_tree_->Branch("event_counter_ytrue", &_event_counter_ytrue_,"event_counter_ytrue/I");
  _sd_tree_->Branch("time_difference_E", &_time_difference_E_,"time_difference_E/D");
  _sd_tree_->Branch("probability", &_internal_probability_,"probability/D");
  _sd_tree_->Branch("length_Emin", &_length_Emin_,"length_Emin/D");
  _sd_tree_->Branch("length_Emax", &_length_Emax_,"length_Emax/D");
  _sd_tree_->Branch("energy_sum", &_energy_,"energy/D");
  _sd_tree_->Branch("minimal_energy", &_minimal_energy_,"minimal_energy/D");
  _sd_tree_->Branch("maximal_energy", &_maximal_energy_,"maximal_energy/D");
  _sd_tree_->Branch("time_Emin", &_time_Emin_,"time_Emin/D");
  _sd_tree_->Branch("time_Emax", &_time_Emax_,"time_Emax/D");
  _sd_tree_->Branch("sigma_time_Emin", &_sigma_time_Emin_,"sigma_time_Emin/D");
  _sd_tree_->Branch("sigma_time_Emax", &_sigma_time_Emax_,"sigma_time_Emax/D");
  _sd_tree_->Branch("event_number", &_event_number_,"event_number/I");

  this->_set_initialized(true);
}

//output file
std::string const final_rate("final_rate.txt");
std::ofstream final_flux(final_rate.c_str());

// Process
dpp::base_module::process_status
TimeDifference::process(datatools::things& data_record_) {
  DT_THROW_IF(! is_initialized(), std::logic_error,
              "Module '" << get_name () << "' is not initialized !");

  //Counting the number of simulated events
  _number_event_++;
  
  ////Storing data bases
  //Simulated data base
  const std::string & sd_label = snemo::datamodel::data_info::default_simulated_data_label();
  DT_THROW_IF(! data_record_.has(sd_label), std::logic_error,
              "Data has no SD !");
  const mctools::simulated_data &a_sd
    = data_record_.get<mctools::simulated_data>(sd_label);
  const mctools::simulated_data::primary_event_type & a_primary_event
    = a_sd.get_primary_event();
  const genbb::primary_event::particles_col_type & the_primary_particles
    = a_primary_event.get_particles();


  //Topology data base
  DT_THROW_IF(! data_record_.has("TD"), std::logic_error,
              "Data has no TD !");
  const snemo::datamodel::topology_data &a_td
    = data_record_.get<snemo::datamodel::topology_data>("TD");

  //Particle track data base
  DT_THROW_IF(! data_record_.has("PTD"), std::logic_error,
              "Data has no PTD !");
  const snemo::datamodel::particle_track_data &a_ptd
    = data_record_.get<snemo::datamodel::particle_track_data>("PTD");


  ////Applying cuts on data banks
  //Cut on PTD bank
  bool vertex_on_source_foil = false;
  bool vertex_on_internal_pads_bulk = false;
  const snemo::datamodel::particle_track_data::particle_collection_type & the_particles
     = a_ptd.get_particles();

  // std::cout << "nbr particles in PTD bank = " << the_particles.size() << std::endl;

  if (abs(a_sd.get_vertex().y())<2371.5){
   _number_event_ytrue_++;
   vertex_on_internal_pads_bulk = true;
   for (unsigned part_i=0; part_i<the_particles.size(); part_i++){
    const snemo::datamodel::particle_track &part_track = a_ptd.get_particle(part_i);
    const snemo::datamodel::particle_track::vertex_collection_type &part_vertices = part_track.get_vertices();
    // std::cout << "nbr vertices per particle in PTD bank = " << part_vertices.size() << std::endl;
    for (unsigned vtx_i=0; vtx_i<part_vertices.size(); vtx_i++){

     if (snemo::datamodel::particle_track::vertex_is_on_source_foil(part_vertices[vtx_i].get())){
     // std::cout << "vertex on source foil !!" << std::endl;
     vertex_on_source_foil = true;
     // const geomtools::blur_spot &foil_vertex = part_vertices[vtx_i].get();
     // if (foil_vertex.get_position()[1] < 2371.5){
     }
    }
   }
  }

  //Cut on TD bank
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
    bool part_hits_are_on_same_calo = false;
    bool hit_calo_double_count = false;

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

     for (unsigned i=0; i < a_calorimeter_max.size(); i++) {

      for (unsigned j=0; j < a_calorimeter_min.size(); j++) {

        const geomtools::geom_id &part_hit_geom_id_min = a_calorimeter_min[j].get().get_geom_id();

        const geomtools::geom_id &part_hit_geom_id_max = a_calorimeter_max[i].get().get_geom_id();
        // std::cout << "geom id max = " << part_hit_geom_id_max << std::endl;
        if (part_hit_geom_id_min == part_hit_geom_id_max) {
          part_hits_are_on_same_calo = true;
        }
        if (a_calorimeter_max.size() > 1 || a_calorimeter_min.size() > 1){
        hit_calo_double_count = true;
        }
       }
      }

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

    //SD bank
    int nb_gen_electron = 0;
    double a_gen_time_difference = 0;
    for (const auto & iparticle : the_primary_particles) {

      if (iparticle.is_electron()) {
	nb_gen_electron++;
      }

      const double a_gen_time
	= iparticle.get_time();
      _gen_particle_time_ = a_gen_time;
      if (nb_gen_electron == 2) {
	a_gen_time_difference = fabs(a_gen_time_difference - a_gen_time);
      }
    }


    //if (a_energy_sum/CLHEP::MeV >= 2.7 && a_energy_sum/CLHEP::MeV <= 3.2) {  
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
    //}

    ////Storing data
    if (my_energy_sum != 0 && my_internal_probability != 0 && !part_hits_are_on_same_calo && vertex_on_source_foil && vertex_on_internal_pads_bulk && !hit_calo_double_count) {//Guarantee we entered in the TD cut loop
      //Keep interesting events in a root tree
      if (a_td.has_pattern_as<snemo::datamodel::topology_2e_pattern>()) {
	_sd_output_file_->cd();
	_gen_time_= a_gen_time_difference/CLHEP::nanosecond;
	_event_counter_ = _number_event_-1;
	_event_counter_ytrue_ = _number_event_ytrue_-1;
	_internal_probability_ = my_internal_probability;
	_length_Emin_ = my_length_Emin/CLHEP::millimeter;
	_length_Emax_ = my_length_Emax/CLHEP::millimeter;
	_energy_ = my_energy_sum/CLHEP::MeV;
	_minimal_energy_ = my_minimal_energy/CLHEP::MeV;
	_maximal_energy_ = my_maximal_energy/CLHEP::MeV;
	_time_Emin_ = my_time_Emin/CLHEP::nanosecond;
	_time_Emax_ = my_time_Emax/CLHEP::nanosecond;
	_time_difference_E_ = fabs(my_time_Emax - my_time_Emin)/CLHEP::nanosecond;
	_sigma_time_Emin_ = my_sigma_time_Emin/CLHEP::nanosecond;
	_sigma_time_Emax_ = my_sigma_time_Emax/CLHEP::nanosecond;
	_sd_tree_->Fill();
      }

      DT_THROW_IF(! final_flux, std::logic_error,
              "ERROR: cannot open the final_rate.txt file!");
       final_flux << "Event # " << _event_counter_ << std::endl;
   }
  }

  return PROCESS_OK;
}

// Reset
void TimeDifference::reset() {
  // Root tree
  _sd_output_file_->cd();
  _sd_output_file_->Write();
  _sd_output_file_->Close();

  this->_set_initialized(false);
}

void TimeDifference::_set_defaults() {
}
