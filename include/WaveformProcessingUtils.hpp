#ifndef WAVEFORMPROCESSOR_H
#define WAVEFORMPROCESSOR_H

#include "PlottingUtils.hpp"
#include <TArrayF.h>
#include <TArrayS.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

struct WaveformFeatures {
  Short_t raw_pulse_height;
  Float_t pulse_height;
  Int_t peak_position;
  Int_t trigger_position;
  Float_t short_integral;
  Float_t long_integral;
  Float_t negative_fraction;
  Bool_t passes_cuts;
  ULong64_t timestamp;
};

struct ProcessingStats {
  Int_t total_processed = 0;
  Int_t accepted = 0;
  Int_t rejected_no_trigger = 0;
  Int_t rejected_insufficient_samples = 0;
  Int_t rejected_negative_integral = 0;
  Int_t rejected_baseline = 0;
  Int_t rejected_clipped = 0;
  Double_t sum_baseline_rms = 0.0;
  Int_t baseline_rms_count = 0;
  Double_t sum_baseline_rms_accepted = 0.0;
  Int_t baseline_rms_count_accepted = 0;
};

enum class InputFormat { kCOMPASS, kWAVEDUMP };

struct FileProcessingConfig {
  Int_t polarity = -1;
  Float_t trigger_threshold = 0.15;
  Int_t num_samples_baseline = 10;
  Int_t pre_samples = 17;
  Int_t post_samples = 190;
  Int_t pre_gate = 10;
  Int_t short_gate = 10;
  Int_t long_gate = 200;
  Int_t sample_waveforms_to_save = 5;
  Int_t max_events = -1;
  Bool_t verbose = kTRUE;
  Bool_t store_waveforms = kTRUE;
  InputFormat input_format = InputFormat::kCOMPASS;
};

class WaveformProcessingUtils {
private:
  Int_t polarity_;
  Double_t trigger_threshold_;
  Int_t num_samples_baseline_;
  Int_t pre_samples_;
  Int_t post_samples_;
  Int_t pre_gate_;
  Int_t short_gate_;
  Int_t long_gate_;
  Int_t max_events_;
  Bool_t verbose_;

  static std::mutex canvas_mutex_;
  Int_t sample_waveforms_to_save_;
  Int_t sample_waveforms_saved_;
  TString current_output_name_;

  ProcessingStats stats_;

  Float_t current_baseline_rms_ = 0.0f;
  Bool_t current_baseline_rms_valid_ = kFALSE;

  TFile *output_file_;
  TTree *output_tree_;
  WaveformFeatures current_features_;
  Bool_t store_waveforms_;
  TArrayF *save_waveform_;
  ULong64_t current_timestamp_;
  InputFormat input_format_;

public:
  WaveformProcessingUtils();
  WaveformProcessingUtils(const FileProcessingConfig &config);
  ~WaveformProcessingUtils();

  void SetPolarity(const Int_t polarity) { polarity_ = polarity; }
  void SetTriggerThreshold(Double_t threshold) {
    trigger_threshold_ = threshold;
  }
  void SetNumberOfSamplesForBaseline(Int_t num_samples_baseline) {
    num_samples_baseline_ = num_samples_baseline;
  }
  void SetSampleWindows(Int_t pre_samples, Int_t post_samples) {
    pre_samples_ = pre_samples;
    post_samples_ = post_samples;
  }
  void SetGates(Int_t pre_gate, Int_t short_gate, Int_t long_gate) {
    pre_gate_ = pre_gate;
    short_gate_ = short_gate;
    long_gate_ = long_gate;
  }
  void SetMaxEvents(Int_t max_events) { max_events_ = max_events; }
  void SetVerbose(Bool_t verbose) { verbose_ = verbose; }
  void SetStoreWaveforms(Bool_t store = kTRUE) { store_waveforms_ = store; }
  void SetSaveSampleWaveforms(Int_t count) {
    sample_waveforms_to_save_ = count;
  }

  Bool_t ProcessWaveform(const TArrayS &samples);

  void SubtractBaseline(const TArrayS &samples);
  Float_t FindTrigger(const TArrayF &waveform);
  void CropWaveform(const TArrayF &waveform, Int_t trigger_pos);
  WaveformFeatures ExtractFeatures(const TArrayF &cropped_wf);
  Bool_t ApplyQualityCuts(const WaveformFeatures &features);
  void SaveSampleWaveform(const TArrayF &waveform);

  void PrintAllStatistics() const;
  ProcessingStats GetStats() const { return stats_; };

  Bool_t ProcessFile(const TString filepath, const TString output_name);

  static void ProcessFilesParallel(const std::vector<TString> &filepaths,
                                   const std::vector<TString> &output_names,
                                   const FileProcessingConfig &config,
                                   Int_t max_workers = 4);
};

#endif
