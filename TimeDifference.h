//! \file    TimeDifference.h
//! \brief   study of time difference for 2 internal electrons for 208Tl
//! \details Process a things object
#ifndef TIMEDIFFERENCE_HH
#define TIMEDIFFERENCE_HH
// Standard Library
#include <vector>
// Third Party
// - Bayeux
#include <TFile.h>
#include <TTree.h>

#include "bayeux/dpp/base_module.h"
#include <bayeux/geomtools/blur_spot.h>

class TimeDifference : public dpp::base_module {

public:

  //! Construct module
  TimeDifference();

  //! Destructor
  virtual ~TimeDifference();

  //! Configure the module
  virtual void initialize(const datatools::properties& myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);

  //! Process supplied data record
  virtual dpp::base_module::process_status process(datatools::things& workItem);

  //! Reset the module
  virtual void reset();

protected:

  /// Give default values to specific class members
  void _set_defaults();

private:

  // The root output file:
TFile * _sd_output_file_;
std::string _output_file_name_;

  // The root tree:
  TTree * _sd_tree_;
  double _time_;
  double _energy_;
  double _minimal_energy_;
  double _maximal_energy_;
  double _internal_probability_;
  double _length_Emin_;
  double _length_Emax_;
  double _time_Emin_;
  double _time_Emax_;
  double _time_difference_E_;
  double _sigma_time_Emin_;
  double _sigma_time_Emax_;
  int _event_counter_;
  int _event_counter_ytrue_;

  std::string _particle_label_;
  double _particle_time_;
  int _number_event_ = 0;
  int _number_event_ytrue_ = 0;
  int _nb_internal_conversion_ = 0;
  int _nb_other_process_ = 0;
  int _nb_2e_topology_ = 0;
  int _counter_test_ = 0;
  int _event_number_ = 0;
   
  geomtools::blur_spot _vertex_; 
  // Macro which automatically creates the interface needed
  // to enable the module to be loaded at runtime
  DPP_MODULE_REGISTRATION_INTERFACE(TimeDifference)
};
#endif // TIMEDIFFERENCE_HH
