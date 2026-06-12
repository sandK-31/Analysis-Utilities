#include "WaveformProcessingUtils.hpp"
#include "IOUtils.hpp"

WaveformProcessingUtils::WaveformProcessingUtils()
    : polarity_(1), trigger_threshold_(0.15), num_samples_baseline_(10),
      pre_samples_(10), post_samples_(100), max_events_(-1), verbose_(kFALSE),
      sample_waveforms_to_save_(0), sample_waveforms_saved_(0),
      output_file_(nullptr), output_tree_(nullptr), store_waveforms_(kTRUE),
      save_waveform_(new TArrayF()), input_format_(InputFormat::kCOMPASS) {}

WaveformProcessingUtils::WaveformProcessingUtils(
    const FileProcessingConfig &config)
    : polarity_(config.polarity), trigger_threshold_(config.trigger_threshold),
      num_samples_baseline_(config.num_samples_baseline),
      pre_samples_(config.pre_samples), post_samples_(config.post_samples),
      pre_gate_(config.pre_gate), short_gate_(config.short_gate),
      long_gate_(config.long_gate), max_events_(config.max_events),
      verbose_(config.verbose),
      sample_waveforms_to_save_(config.sample_waveforms_to_save),
      sample_waveforms_saved_(0), output_file_(nullptr), output_tree_(nullptr),
      store_waveforms_(config.store_waveforms), save_waveform_(new TArrayF()),
      input_format_(config.input_format) {}

WaveformProcessingUtils::~WaveformProcessingUtils() {
  if (output_file_) {
    if (output_file_->IsOpen()) {
      output_file_->Close();
    }
    delete output_file_;
    output_file_ = nullptr;
    if (store_waveforms_) {
      save_waveform_ = nullptr;
    }
  }
  delete save_waveform_;
  save_waveform_ = nullptr;
}

Bool_t WaveformProcessingUtils::ProcessWaveform(const TArrayS &samples) {
  Int_t n = samples.GetSize();

  if (n == 0)
    return kFALSE;

  Short_t raw_max = samples.At(0);

  if (polarity_ == 1) {
    for (Int_t i = 1; i < n; ++i) {
      if (samples[i] > raw_max)
        raw_max = samples[i];
    }
  } else {
    for (Int_t i = 1; i < n; ++i) {
      if (samples[i] < raw_max)
        raw_max = samples[i];
    }
  }

  SubtractBaseline(samples);

  Float_t trigger_pos = FindTrigger(*save_waveform_);
  if (trigger_pos < 0) {
    stats_.rejected_no_trigger++;
    return kFALSE;
  }

  if (trigger_pos < pre_samples_ ||
      (save_waveform_->GetSize() - trigger_pos) <= post_samples_) {
    stats_.rejected_insufficient_samples++;
    return kFALSE;
  }

  CropWaveform(*save_waveform_, trigger_pos);

  WaveformFeatures features = ExtractFeatures(*save_waveform_);
  features.raw_pulse_height = std::abs(raw_max);
  features.trigger_position = FindTrigger(*save_waveform_);
  Bool_t passes_cuts = ApplyQualityCuts(features);
  features.passes_cuts = passes_cuts;

  if (!passes_cuts) {
    return kFALSE;
  }

  if (sample_waveforms_saved_ < sample_waveforms_to_save_) {
    SaveSampleWaveform(*save_waveform_);
  }

  current_features_ = features;
  output_tree_->Fill();

  stats_.accepted++;
  if (current_baseline_rms_valid_) {
    stats_.sum_baseline_rms_accepted += current_baseline_rms_;
    stats_.baseline_rms_count_accepted++;
  }
  return kTRUE;
}

std::mutex WaveformProcessingUtils::canvas_mutex_;

void WaveformProcessingUtils::SaveSampleWaveform(const TArrayF &waveform) {

  std::lock_guard<std::mutex> lock(canvas_mutex_);

  Int_t n = waveform.GetSize();
  const Float_t *arr = waveform.GetArray();
  std::vector<Double_t> x(n), y(n);
  for (Int_t i = 0; i < n; ++i) {
    x[i] = i;
    y[i] = arr[i];
  }

  PlottingUtils::SetStylePreferences(PlotSaveFormat::kPNG);

  TGraph *graph = new TGraph(n, x.data(), y.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);

  PlottingUtils::ConfigureGraph(graph, kBlue + 1, ";Sample;Amplitude [ADC]");
  graph->Draw("AL");

  TString output_name = Form("%s_waveform_%04d", current_output_name_.Data(),
                             sample_waveforms_saved_);

  PlottingUtils::SaveFigure(canvas, output_name, "samplewaveforms",
                            PlotSaveOptions::kLINEAR);
  delete graph;
  delete canvas;

  sample_waveforms_saved_++;
}

void WaveformProcessingUtils::SubtractBaseline(const TArrayS &samples) {
  Int_t n = samples.GetSize();

  Float_t baseline = 0;
  Int_t baseline_samples = TMath::Min(num_samples_baseline_, n);
  for (Int_t i = 0; i < baseline_samples; ++i) {
    baseline += samples.GetAt(i);
  }
  baseline /= baseline_samples;

  current_baseline_rms_valid_ = kFALSE;
  if (baseline_samples > 1) {
    Double_t sum_sq = 0.0;
    for (Int_t i = 0; i < baseline_samples; ++i) {
      Double_t d = samples.GetAt(i) - baseline;
      sum_sq += d * d;
    }
    Float_t rms = TMath::Sqrt(sum_sq / (baseline_samples - 1));
    stats_.sum_baseline_rms += rms;
    stats_.baseline_rms_count++;
    current_baseline_rms_ = rms;
    current_baseline_rms_valid_ = kTRUE;
  }

  save_waveform_->Set(n);
  if (polarity_ == -1) {
    for (Int_t i = 0; i < n; ++i) {
      save_waveform_->SetAt(baseline - samples.GetAt(i), i);
    }
  } else {
    for (Int_t i = 0; i < n; ++i) {
      save_waveform_->SetAt(samples.GetAt(i) - baseline, i);
    }
  }
}

Float_t WaveformProcessingUtils::FindTrigger(const TArrayF &waveform) {
  Int_t n = waveform.GetSize();
  const Float_t *arr = waveform.GetArray();

  Float_t peak_value = *std::max_element(arr, arr + n);
  Float_t trigger_level = peak_value * trigger_threshold_;

  for (Int_t i = 0; i < n; ++i) {
    if (arr[i] >= trigger_level) {
      return Float_t(i);
    }
  }

  return -1.0;
}

void WaveformProcessingUtils::CropWaveform(const TArrayF &waveform,
                                           Int_t trigger_pos) {
  Int_t start = trigger_pos - pre_samples_;
  Int_t end = TMath::Min(trigger_pos + post_samples_, waveform.GetSize());
  Int_t crop_size = end - start;

  TArrayF cropped(crop_size);
  const Float_t *src = waveform.GetArray();
  for (Int_t i = 0; i < crop_size; ++i) {
    cropped[i] = src[start + i];
  }

  *save_waveform_ = cropped;
}

WaveformFeatures
WaveformProcessingUtils::ExtractFeatures(const TArrayF &cropped_wf) {
  WaveformFeatures features;
  Int_t integration_start = pre_samples_ - pre_gate_;

  Int_t n = cropped_wf.GetSize();
  const Float_t *arr = cropped_wf.GetArray();

  const Float_t *max_it = std::max_element(arr, arr + n);
  features.pulse_height = *max_it;
  features.peak_position = std::distance(arr, max_it);

  features.short_integral = 0;
  features.long_integral = 0;

  Int_t negative_samples = 0;
  Int_t short_end = TMath::Min(integration_start + short_gate_, n);
  Int_t long_end = TMath::Min(integration_start + long_gate_, n);

  for (Int_t i = integration_start; i < long_end; ++i) {
    Float_t sample_value = arr[i];
    features.long_integral += sample_value;
    if (i < short_end) {
      features.short_integral += sample_value;
    }
    if (sample_value < 0)
      negative_samples++;
  }
  features.timestamp = current_timestamp_;

  features.passes_cuts = kTRUE;
  features.negative_fraction =
      Float_t(negative_samples) / Float_t(long_end - integration_start);

  return features;
}

Bool_t
WaveformProcessingUtils::ApplyQualityCuts(const WaveformFeatures &features) {

  if (((features.raw_pulse_height == 16384) && (polarity_ == 1)) ||
      ((features.raw_pulse_height == 0) && (polarity_ == -1))) {
    stats_.rejected_clipped++;
    return kFALSE;
  }

  if (features.negative_fraction > 0.50) {
    stats_.rejected_baseline++;
    return kFALSE;
  }

  if (features.long_integral <= 0) {
    stats_.rejected_negative_integral++;
    return kFALSE;
  }

  return kTRUE;
}

void WaveformProcessingUtils::PrintAllStatistics() const {
  std::cout << "Waveform processing statistics..." << std::endl;
  std::cout << "Total processed: " << stats_.total_processed << std::endl;
  std::cout << std::endl;
  std::cout << "Accepted: " << stats_.accepted << std::endl;
  std::cout << std::endl;
  std::cout << "Rejected no trigger: " << stats_.rejected_no_trigger
            << std::endl;
  std::cout << "Rejected clipped ADC: " << stats_.rejected_clipped << std::endl;
  std::cout << "Rejected insufficient samples: "
            << stats_.rejected_insufficient_samples << std::endl;
  std::cout << "Rejected negative integral: "
            << stats_.rejected_negative_integral << std::endl;
  std::cout << "Rejected bad baseline: " << stats_.rejected_baseline
            << std::endl;
  std::cout << std::endl;

  if (stats_.total_processed > 0) {
    std::cout << "Acceptance rate: "
              << 100 * Float_t(stats_.accepted) /
                     Float_t(stats_.total_processed)
              << "%" << std::endl;
  }
  if (stats_.baseline_rms_count > 0) {
    std::cout << "Mean baseline RMS (all processed): "
              << stats_.sum_baseline_rms / stats_.baseline_rms_count
              << " ADC counts (over " << stats_.baseline_rms_count
              << " waveforms)" << std::endl;
  }
  if (stats_.baseline_rms_count_accepted > 0) {
    std::cout << "Mean baseline RMS (accepted only): "
              << stats_.sum_baseline_rms_accepted /
                     stats_.baseline_rms_count_accepted
              << " ADC counts (over " << stats_.baseline_rms_count_accepted
              << " waveforms)" << std::endl;
  }
  std::cout << std::endl;

  TString stats_path =
      IO::GetRootFilesBaseDir() + "/" + current_output_name_ + ".stats";
  std::ofstream stats_file(stats_path.Data(), std::ios::app);
  if (stats_file.is_open()) {
    stats_file << "Waveform processing statistics..." << std::endl;
    stats_file << "Total processed: " << stats_.total_processed << std::endl;
    stats_file << std::endl;
    stats_file << "Accepted: " << stats_.accepted << std::endl;
    stats_file << std::endl;
    stats_file << "Rejected no trigger: " << stats_.rejected_no_trigger
               << std::endl;
    stats_file << "Rejected clipped ADC: " << stats_.rejected_clipped
               << std::endl;
    stats_file << "Rejected insufficient samples: "
               << stats_.rejected_insufficient_samples << std::endl;
    stats_file << "Rejected negative integral: "
               << stats_.rejected_negative_integral << std::endl;
    stats_file << "Rejected bad baseline: " << stats_.rejected_baseline
               << std::endl;
    stats_file << std::endl;

    if (stats_.total_processed > 0) {
      stats_file << "Acceptance rate: "
                 << 100 * Float_t(stats_.accepted) /
                        Float_t(stats_.total_processed)
                 << "%" << std::endl;
    }
    if (stats_.baseline_rms_count > 0) {
      stats_file << "Mean baseline RMS (all processed): "
                 << stats_.sum_baseline_rms / stats_.baseline_rms_count
                 << " ADC counts (over " << stats_.baseline_rms_count
                 << " waveforms)" << std::endl;
    }
    if (stats_.baseline_rms_count_accepted > 0) {
      stats_file << "Mean baseline RMS (accepted only): "
                 << stats_.sum_baseline_rms_accepted /
                        stats_.baseline_rms_count_accepted
                 << " ADC counts (over " << stats_.baseline_rms_count_accepted
                 << " waveforms)" << std::endl;
    }
    stats_file << std::endl;
  }
}

Bool_t WaveformProcessingUtils::ProcessFile(const TString filepath,
                                            const TString output_name) {
  current_output_name_ = output_name;
  sample_waveforms_saved_ = 0;
  if (!save_waveform_) {
    save_waveform_ = new TArrayF();
  }

  const TString base_dir = IO::GetRootFilesBaseDir();
  if (gSystem->AccessPathName(base_dir)) {
    gSystem->mkdir(base_dir, kTRUE);
  }

  // clear file
  TString clear_path = base_dir + "/" + output_name + ".stats";
  std::ofstream(clear_path.Data(), std::ios::trunc);

  TString output_subpath = output_name + ".root";
  TString output_filename = base_dir + "/" + output_subpath;
  output_file_ = IO::OpenForWriting(output_subpath);
  if (!output_file_ || output_file_->IsZombie()) {
    std::cout << "ERROR: Could not create output file " << output_filename
              << std::endl;
    return kFALSE;
  }

  output_tree_ = new TTree("features", "Waveform Features");

  output_tree_->Branch("pulse_height", &current_features_.pulse_height,
                       "pulse_height/F");
  output_tree_->Branch("trigger_position", &current_features_.trigger_position,
                       "trigger_position/I");
  output_tree_->Branch("short_integral", &current_features_.short_integral,
                       "short_integral/F");
  output_tree_->Branch("long_integral", &current_features_.long_integral,
                       "long_integral/F");
  output_tree_->Branch("timestamp", &current_features_.timestamp,
                       "timestamp/l");

  if (store_waveforms_) {
    output_tree_->Branch("Samples", &save_waveform_);
    std::cout << "Storing events that pass cuts." << std::endl;
  }

  TFile *file = TFile::Open(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "ERROR opening file: " << filepath << std::endl;
    return kFALSE;
  }

  TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
  if (!tree) {
    std::cout << "ERROR: TTree 'Data_R' not found in " << filepath << std::endl;
    file->Close();
    return kFALSE;
  }

  TArrayS *samples = new TArrayS();
  tree->SetBranchAddress("Samples", &samples);

  UInt_t trigger_time_tag = 0;
  if (input_format_ == InputFormat::kCOMPASS) {
    tree->SetBranchAddress("Timestamp", &current_timestamp_);
  } else if (input_format_ == InputFormat::kWAVEDUMP) {
    tree->SetBranchAddress("TriggerTimeTag", &trigger_time_tag);
  }

  Long64_t n_entries = tree->GetEntries();

  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    if (max_events_ > 0 && stats_.accepted >= max_events_) {
      break;
    }
    tree->GetEntry(entry);
    if (input_format_ == InputFormat::kWAVEDUMP) {
      current_timestamp_ = static_cast<ULong64_t>(trigger_time_tag);
    }
    stats_.total_processed++;
    ProcessWaveform(*samples);
  }

  delete samples;
  file->Close();
  delete file;

  output_file_->cd();
  output_tree_->Write("", TObject::kOverwrite);
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;
  output_tree_ = nullptr;
  if (store_waveforms_) {
    save_waveform_ = nullptr;
  }

  if (verbose_) {
    PrintAllStatistics();
  }

  return kTRUE;
}

void WaveformProcessingUtils::ProcessFilesParallel(
    const std::vector<TString> &filepaths,
    const std::vector<TString> &output_names,
    const FileProcessingConfig &config, Int_t max_workers) {

  ROOT::EnableThreadSafety();
  IO::SetThreadSafe(kTRUE);

  Int_t n_files = Int_t(filepaths.size());
  Int_t n_workers = max_workers > 0
                        ? max_workers
                        : Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Processing " << n_files << " files with " << n_workers
            << " workers." << std::endl;

  std::function<Bool_t(const TString &, const TString &)> process_one =
      [&config](const TString &filepath, const TString &output_name) -> Bool_t {
    WaveformProcessingUtils *processor = new WaveformProcessingUtils(config);
    Bool_t result = processor->ProcessFile(filepath, output_name);
    delete processor;
    return result;
  };

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);

    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(std::async(std::launch::async, process_one,
                                   std::cref(filepaths[j]),
                                   std::cref(output_names[j])));
    }

    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      std::cout << "Finished: " << output_names[i + j]
                << (result ? " [OK]" : " [FAILED]") << std::endl;
    }
  }
}
