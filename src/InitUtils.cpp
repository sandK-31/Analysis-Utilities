#include "InitUtils.hpp"

namespace {

// Compute the fixed per-event record size (in bytes) for a CoMPASS file given
// the global header's control bits. Only valid when no waveform is present.
Int_t ComputeCoMPASSRecordSize(UShort_t global_header) {
  Int_t record_size = 16; // board(2) + channel(2) + timestamp(8) + flags(4)
  if (global_header & 0x0001)
    record_size += 2; // energy_ch
  if (global_header & 0x0002)
    record_size += 8; // energy_cal
  if (global_header & 0x0004)
    record_size += 2; // energy_short
  return record_size;
}

// Compute the per-field byte offsets within a fixed-stride CoMPASS record.
// Offsets to optional fields are set to -1 when not present.
void ComputeCoMPASSFieldOffsets(UShort_t global_header, Long64_t &off_board,
                                Long64_t &off_channel, Long64_t &off_timestamp,
                                Long64_t &off_energy_ch,
                                Long64_t &off_energy_cal,
                                Long64_t &off_energy_short,
                                Long64_t &off_flags) {
  Bool_t has_energy_ch = (global_header & 0x0001);
  Bool_t has_energy_cal = (global_header & 0x0002);
  Bool_t has_energy_short = (global_header & 0x0004);

  off_board = 0;
  off_channel = 2;
  off_timestamp = 4;
  off_energy_ch = -1;
  off_energy_cal = -1;
  off_energy_short = -1;

  Long64_t cursor = 12;
  if (has_energy_ch) {
    off_energy_ch = cursor;
    cursor += 2;
  }
  if (has_energy_cal) {
    off_energy_cal = cursor;
    cursor += 8;
  }
  if (has_energy_short) {
    off_energy_short = cursor;
    cursor += 2;
  }
  off_flags = cursor;
}

// Load the data section of a fixed-stride CoMPASS file into a buffer using
// a single fread. Returns kFALSE and an empty buffer on error.
Bool_t LoadCoMPASSBulkBuffer(const TString &input_filename,
                             Long64_t header_bytes_to_skip, Int_t record_size,
                             std::vector<char> &buf, Long64_t &n_events,
                             Long64_t &file_size) {
  buf.clear();
  n_events = 0;
  file_size = 0;

  FILE *fp = std::fopen(input_filename.Data(), "rb");
  if (!fp) {
    std::cout << "ERROR: Failed to open file for bulk read: " << input_filename
              << std::endl;
    return kFALSE;
  }

  if (std::fseek(fp, 0, SEEK_END) != 0) {
    std::cout << "ERROR: fseek to end failed on " << input_filename
              << std::endl;
    std::fclose(fp);
    return kFALSE;
  }
  Long64_t fsize = Long64_t(std::ftell(fp));
  if (fsize < 0) {
    std::cout << "ERROR: ftell failed on " << input_filename << std::endl;
    std::fclose(fp);
    return kFALSE;
  }
  file_size = fsize;

  Long64_t data_bytes = file_size - header_bytes_to_skip;
  if (data_bytes < 0)
    data_bytes = 0;

  Long64_t record_size_l = Long64_t(record_size);
  if (data_bytes % record_size_l != 0) {
    std::cout << "WARNING: File data section (" << data_bytes
              << " bytes) is not a multiple of record size " << record_size
              << "; truncating to floor." << std::endl;
  }
  n_events = data_bytes / record_size_l;
  Long64_t want = n_events * record_size_l;

  if (std::fseek(fp, long(header_bytes_to_skip), SEEK_SET) != 0) {
    std::cout << "ERROR: fseek to data start failed on " << input_filename
              << std::endl;
    std::fclose(fp);
    return kFALSE;
  }

  buf.resize(size_t(want));
  size_t got = (want > 0) ? std::fread(buf.data(), 1, size_t(want), fp) : 0;
  std::fclose(fp);

  if (Long64_t(got) != want) {
    std::cout << "ERROR: Short read on " << input_filename << " (got "
              << Long64_t(got) << " of " << want << " bytes)" << std::endl;
    buf.clear();
    n_events = 0;
    return kFALSE;
  }

  return kTRUE;
}

// Print the same per-file header summary that CoMPASSData::PrintHeader produces
// for the first event in the slow path. Called once from the bulk-read paths.
void PrintCoMPASSHeaderSummary(UShort_t global_header) {
  Bool_t has_energy_ch = (global_header & 0x0001);
  Bool_t has_energy_cal = (global_header & 0x0002);
  Bool_t has_energy_short = (global_header & 0x0004);

  std::cout << "CoMPASS event header..." << std::endl;
  std::cout << "Header:    0x" << std::hex << global_header << std::dec
            << " (Binary: " << std::bitset<16>(global_header) << ")"
            << std::endl;
  std::cout << std::endl;
  std::cout << "Control bits..." << std::endl;
  std::cout << "Energy (ch):    " << (has_energy_ch ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (cal):   " << (has_energy_cal ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (short): " << (has_energy_short ? "YES" : "NO")
            << std::endl;
  std::cout << "Waveform:       NO" << std::endl;
}

// Print the summary lines (event counts, warning aggregations, bytes read) that
// terminate both ConvertCoMPASS* paths. Kept here so the slow and fast paths
// emit identical output.
void PrintCoMPASSConversionSummary(
    Long64_t event_count, Long64_t bytes_read, Bool_t skip_bad_events,
    Long64_t warning_fake, Long64_t warning_saturated, Long64_t warning_pileup,
    Long64_t warning_memory_full, Long64_t warning_trigger_lost,
    Long64_t warning_pll_loss, Long64_t warning_over_temp,
    Long64_t warning_adc_shutdown) {
  std::cout << "Conversion complete." << std::endl;
  std::cout << "Total events processed: " << event_count << std::endl;

  if (warning_fake > 0 || warning_saturated > 0 || warning_pileup > 0) {
    std::cout << "Events with rejection-quality flags:" << std::endl;
    if (warning_fake > 0) {
      std::cout << "  Fake events: " << warning_fake;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_saturated > 0) {
      std::cout << "  Saturated: " << warning_saturated;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_pileup > 0) {
      std::cout << "  Pileup: " << warning_pileup;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  if (warning_memory_full > 0) {
    std::cout << "WARNING: " << warning_memory_full
              << " events with memory full flag" << std::endl;
  }
  if (warning_trigger_lost > 0) {
    std::cout << "WARNING: " << warning_trigger_lost
              << " events with trigger lost flag" << std::endl;
  }
  if (warning_pll_loss > 0) {
    std::cout << "WARNING: " << warning_pll_loss << " events with PLL lock loss"
              << std::endl;
  }
  if (warning_over_temp > 0) {
    std::cout << "WARNING: " << warning_over_temp
              << " events with over temperature" << std::endl;
  }
  if (warning_adc_shutdown > 0) {
    std::cout << "WARNING: " << warning_adc_shutdown
              << " events with ADC shutdown" << std::endl;
  }

  std::cout << "Total bytes read: " << bytes_read << std::endl;
}

// Bulk-read fast path for ConvertCoMPASSBinToHits. Returns the parsed hits
// vector. The caller must guarantee the file has no waveform records.
std::pair<std::vector<RawHit>, UShort_t>
BulkReadCoMPASSHits(const TString &input_filename, UShort_t global_header,
                    UShort_t global_header_override, Bool_t skip_bad_events) {
  std::vector<RawHit> hits;

  Int_t record_size = ComputeCoMPASSRecordSize(global_header);
  Long64_t header_bytes_to_skip = (global_header_override != 0) ? 0 : 2;

  std::vector<char> buf;
  Long64_t n_events = 0;
  Long64_t file_size = 0;
  if (!LoadCoMPASSBulkBuffer(input_filename, header_bytes_to_skip, record_size,
                             buf, n_events, file_size)) {
    return std::make_pair(hits, static_cast<UShort_t>(0));
  }

  Long64_t off_board = 0, off_channel = 0, off_timestamp = 0;
  Long64_t off_energy_ch = -1, off_energy_cal = -1, off_energy_short = -1;
  Long64_t off_flags = 0;
  ComputeCoMPASSFieldOffsets(global_header, off_board, off_channel,
                             off_timestamp, off_energy_ch, off_energy_cal,
                             off_energy_short, off_flags);

  Long64_t event_count = 0;
  Long64_t warning_fake = 0;
  Long64_t warning_saturated = 0;
  Long64_t warning_pileup = 0;
  Long64_t warning_memory_full = 0;
  Long64_t warning_trigger_lost = 0;
  Long64_t warning_pll_loss = 0;
  Long64_t warning_over_temp = 0;
  Long64_t warning_adc_shutdown = 0;

  std::cout << "Reading events..." << std::endl;
  if (skip_bad_events) {
    std::cout
        << "Filtering enabled: skipping fake, saturated, and pileup events"
        << std::endl;
  }

  Bool_t header_printed = kFALSE;

  hits.reserve(size_t(n_events));

  const char *ptr = buf.data();
  Long64_t record_size_l = Long64_t(record_size);

  for (Long64_t i = 0; i < n_events; i++) {
    if (!header_printed) {
      PrintCoMPASSHeaderSummary(global_header);
      header_printed = kTRUE;
    }

    const char *rec = ptr + i * record_size_l;

    UInt_t flags_val;
    std::memcpy(&flags_val, rec + off_flags, 4);

    Bool_t is_fake = (flags_val & CoMPASSData::FAKE_EVENT) != 0;
    Bool_t is_saturated = ((flags_val & CoMPASSData::INPUT_SATURATING) ||
                           (flags_val & CoMPASSData::SATURATION_IN_GATE)) != 0;
    Bool_t is_pileup = (flags_val & CoMPASSData::PILEUP) != 0;

    // The slow path increments via cascading early-continue when skipping,
    // so each rejected event contributes to only the first matching counter.
    // Replicate that exactly for output parity.
    if (is_fake) {
      warning_fake++;
      if (skip_bad_events)
        continue;
    }
    if (is_saturated) {
      warning_saturated++;
      if (skip_bad_events)
        continue;
    }
    if (is_pileup) {
      warning_pileup++;
      if (skip_bad_events)
        continue;
    }

    if (flags_val & CoMPASSData::MEMORY_FULL)
      warning_memory_full++;
    if (flags_val & CoMPASSData::TRIGGER_LOST)
      warning_trigger_lost++;
    if (flags_val & CoMPASSData::PLL_LOCK_LOSS)
      warning_pll_loss++;
    if (flags_val & CoMPASSData::OVER_TEMPERATURE)
      warning_over_temp++;
    if (flags_val & CoMPASSData::ADC_SHUTDOWN)
      warning_adc_shutdown++;

    RawHit hit;
    std::memcpy(&hit.board, rec + off_board, 2);
    std::memcpy(&hit.channel, rec + off_channel, 2);
    std::memcpy(&hit.timestamp, rec + off_timestamp, 8);
    hit.energy = 0;
    if (off_energy_ch >= 0)
      std::memcpy(&hit.energy, rec + off_energy_ch, 2);
    hit.flags = flags_val;
    hits.push_back(hit);

    event_count++;
  }

  PrintCoMPASSConversionSummary(
      event_count, file_size, skip_bad_events, warning_fake, warning_saturated,
      warning_pileup, warning_memory_full, warning_trigger_lost,
      warning_pll_loss, warning_over_temp, warning_adc_shutdown);

  return std::make_pair(hits, global_header);
}

} // namespace

void InitUtils::SetROOTPreferences(PlotSaveFormat save_format,
                                   const TString &plots_dir,
                                   const TString &root_files_dir,
                                   Bool_t enable_mt) {
  PlottingUtils::SetStylePreferences(save_format);
  gROOT->ForceStyle(kTRUE);
  gROOT->SetBatch(kTRUE);

  if (enable_mt) {
    IO::SetThreadSafe(kTRUE);
    // Detach new histograms from gDirectory so concurrent TFile openings on
    // other threads don't race on the directory's child list.
    TH1::AddDirectory(kFALSE);
  }

  TString resolved_plots_dir = plots_dir;
  if (resolved_plots_dir.Length() == 0) {
    std::cout
        << "WARNING: InitUtils::SetROOTPreferences called without plots_dir; "
           "defaulting to CWD-relative \"plots\". Pass an absolute path to "
           "drive output into a project root."
        << std::endl;
    resolved_plots_dir = "plots";
  }
  PlottingUtils::SetPlotsBaseDir(resolved_plots_dir);

  TString resolved_root_files_dir = root_files_dir;
  if (resolved_root_files_dir.Length() == 0) {
    std::cout << "WARNING: InitUtils::SetROOTPreferences called without "
                 "root_files_dir; defaulting to CWD-relative \"root_files\"."
              << std::endl;
    resolved_root_files_dir = "root_files";
  } else {
    IO::SetRootFilesBaseDir(resolved_root_files_dir);
  }

  if (gSystem->AccessPathName(PlottingUtils::GetPlotsBaseDir())) {
    gSystem->mkdir(PlottingUtils::GetPlotsBaseDir(), kTRUE);
  }
  if (gSystem->AccessPathName(IO::GetRootFilesBaseDir())) {
    gSystem->mkdir(IO::GetRootFilesBaseDir(), kTRUE);
  }
}

UShort_t InitUtils::ConvertCoMPASSBinToROOT(const TString input_filename,
                                            const TString output_name,
                                            UShort_t global_header_override,
                                            Bool_t skip_bad_events) {
  const TString base_dir = IO::GetRootFilesBaseDir();
  if (gSystem->AccessPathName(base_dir)) {
    gSystem->mkdir(base_dir, kTRUE);
  }

  if (gSystem->AccessPathName(input_filename)) {
    std::cout << "ERROR: Input file does not exist: " << input_filename
              << std::endl;
    return 0;
  }

  TString output_subpath = output_name + ".root";
  TString output_filename = base_dir + "/" + output_subpath;

  CoMPASSReader reader;
  Bool_t open_success =
      (global_header_override != 0)
          ? reader.Open(input_filename.Data(), global_header_override)
          : reader.Open(input_filename.Data());

  if (!open_success) {
    std::cout << "ERROR: Failed to open CoMPASS binary file" << std::endl;
    return 0;
  }

  UShort_t global_header = reader.GetGlobalHeader();

  Bool_t has_energy_ch = (global_header & 0x0001);
  Bool_t has_energy_cal = (global_header & 0x0002);
  Bool_t has_energy_short = (global_header & 0x0004);
  Bool_t has_waveform = (global_header & 0x0008);

  TFile *outfile = nullptr;
  TTree *tree = nullptr;
  UShort_t board = 0, channel = 0, energy = 0, energy_short = 0;
  ULong64_t timestamp = 0;
  Double_t energy_cal = 0.0;
  UInt_t flags = 0, num_samples = 0;
  UChar_t waveform_code = 0;
  TArrayS *samples = nullptr;

  {
    IO::ScopedRootLock setup_guard;

    outfile = IO::OpenForWriting(output_subpath);
    if (!outfile || outfile->IsZombie()) {
      std::cout << "ERROR: Could not create output file " << output_filename
                << std::endl;
      reader.Close();
      return 0;
    }

    tree = new TTree("Data_R", "CoMPASS Binary Data");

    tree->Branch("Board", &board, "Board/s");
    tree->Branch("Channel", &channel, "Channel/s");
    tree->Branch("Timestamp", &timestamp, "Timestamp/l");

    if (has_energy_ch) {
      tree->Branch("Energy", &energy, "Energy/s");
      std::cout << "Energy type: Channel (ADC counts)" << std::endl;
    } else if (has_energy_cal) {
      tree->Branch("Energy", &energy_cal, "Energy/D");
      std::cout << "Energy type: Calibrated (keV/MeV)" << std::endl;
    }

    if (has_energy_short) {
      tree->Branch("EnergyShort", &energy_short, "EnergyShort/s");
    }

    tree->Branch("Flags", &flags, "Flags/i");

    if (has_waveform) {
      samples = new TArrayS();
      tree->Branch("WaveformCode", &waveform_code, "WaveformCode/b");
      tree->Branch("NumSamples", &num_samples, "NumSamples/i");
      tree->Branch("Samples", &samples);
    }
  }

  Long64_t event_count = 0;
  Long64_t warning_fake = 0;
  Long64_t warning_saturated = 0;
  Long64_t warning_pileup = 0;
  Long64_t warning_memory_full = 0;
  Long64_t warning_trigger_lost = 0;
  Long64_t warning_pll_loss = 0;
  Long64_t warning_over_temp = 0;
  Long64_t warning_adc_shutdown = 0;
  Long64_t bytes_read_for_summary = 0;

  std::cout << "Reading events..." << std::endl;
  if (skip_bad_events) {
    std::cout
        << "Filtering enabled: skipping fake, saturated, and pileup events"
        << std::endl;
  }

  if (!has_waveform) {
    // Fast bulk-read path: fixed-stride records, one fread + memcpy parse.
    std::cout
        << "Using bulk-read fast path (fixed-stride records, no waveform)."
        << std::endl;
    reader.Close();

    Int_t record_size = ComputeCoMPASSRecordSize(global_header);
    Long64_t header_bytes_to_skip = (global_header_override != 0) ? 0 : 2;

    std::vector<char> buf;
    Long64_t n_events = 0;
    Long64_t file_size = 0;
    Bool_t load_ok =
        LoadCoMPASSBulkBuffer(input_filename, header_bytes_to_skip, record_size,
                              buf, n_events, file_size);
    if (!load_ok) {
      IO::ScopedRootLock teardown_guard;
      outfile->Close();
      delete outfile;
      return 0;
    }
    bytes_read_for_summary = file_size;

    Long64_t off_board = 0, off_channel = 0, off_timestamp = 0;
    Long64_t off_energy_ch = -1, off_energy_cal = -1, off_energy_short = -1;
    Long64_t off_flags = 0;
    ComputeCoMPASSFieldOffsets(global_header, off_board, off_channel,
                               off_timestamp, off_energy_ch, off_energy_cal,
                               off_energy_short, off_flags);

    Bool_t header_printed = kFALSE;

    const char *ptr = buf.data();
    Long64_t record_size_l = Long64_t(record_size);

    for (Long64_t i = 0; i < n_events; i++) {
      if (!header_printed) {
        PrintCoMPASSHeaderSummary(global_header);
        header_printed = kTRUE;
      }

      const char *rec = ptr + i * record_size_l;

      UInt_t flags_val;
      std::memcpy(&flags_val, rec + off_flags, 4);

      Bool_t is_fake = (flags_val & CoMPASSData::FAKE_EVENT) != 0;
      Bool_t is_saturated =
          ((flags_val & CoMPASSData::INPUT_SATURATING) ||
           (flags_val & CoMPASSData::SATURATION_IN_GATE)) != 0;
      Bool_t is_pileup = (flags_val & CoMPASSData::PILEUP) != 0;

      if (is_fake) {
        warning_fake++;
        if (skip_bad_events)
          continue;
      }
      if (is_saturated) {
        warning_saturated++;
        if (skip_bad_events)
          continue;
      }
      if (is_pileup) {
        warning_pileup++;
        if (skip_bad_events)
          continue;
      }

      if (flags_val & CoMPASSData::MEMORY_FULL)
        warning_memory_full++;
      if (flags_val & CoMPASSData::TRIGGER_LOST)
        warning_trigger_lost++;
      if (flags_val & CoMPASSData::PLL_LOCK_LOSS)
        warning_pll_loss++;
      if (flags_val & CoMPASSData::OVER_TEMPERATURE)
        warning_over_temp++;
      if (flags_val & CoMPASSData::ADC_SHUTDOWN)
        warning_adc_shutdown++;

      std::memcpy(&board, rec + off_board, 2);
      std::memcpy(&channel, rec + off_channel, 2);
      std::memcpy(&timestamp, rec + off_timestamp, 8);
      if (off_energy_ch >= 0)
        std::memcpy(&energy, rec + off_energy_ch, 2);
      if (off_energy_cal >= 0)
        std::memcpy(&energy_cal, rec + off_energy_cal, 8);
      if (off_energy_short >= 0)
        std::memcpy(&energy_short, rec + off_energy_short, 2);
      flags = flags_val;

      tree->Fill();
      event_count++;
    }
  } else {
    // Slow path: waveform records are variable-length; keep the existing
    // CoMPASSReader-based loop.
    while (reader.ReadEvent()) {
      const CoMPASSData &event = reader.GetCurrentEvent();
      if (event_count == 0) {
        event.PrintHeader();
      }

      if (event.isFakeEvent()) {
        warning_fake++;
        if (skip_bad_events)
          continue;
      }
      if (event.isInputSaturating() || event.hasSaturation()) {
        warning_saturated++;
        if (skip_bad_events)
          continue;
      }
      if (event.isPileup()) {
        warning_pileup++;
        if (skip_bad_events)
          continue;
      }

      if (event.hasMemoryFull())
        warning_memory_full++;
      if (event.hasTriggerLost())
        warning_trigger_lost++;
      if (event.hasPLLLockLoss())
        warning_pll_loss++;
      if (event.isOverTemperature())
        warning_over_temp++;
      if (event.isADCShutdown())
        warning_adc_shutdown++;

      board = event.board;
      channel = event.channel;
      timestamp = event.timestamp;
      flags = event.flags;

      if (has_energy_ch) {
        energy = event.energy_ch;
      }
      if (has_energy_cal) {
        energy_cal = event.energy_cal;
      }
      if (has_energy_short) {
        energy_short = event.energy_short_ch;
      }
      if (has_waveform) {
        waveform_code = event.waveform_code;
        num_samples = event.num_samples;
        *samples = event.samples;
      }

      tree->Fill();
      event_count++;
    }
    bytes_read_for_summary = reader.GetBytesRead();
  }

  PrintCoMPASSConversionSummary(
      event_count, bytes_read_for_summary, skip_bad_events, warning_fake,
      warning_saturated, warning_pileup, warning_memory_full,
      warning_trigger_lost, warning_pll_loss, warning_over_temp,
      warning_adc_shutdown);

  {
    IO::ScopedRootLock teardown_guard;
    outfile->cd();
    tree->Write("", TObject::kOverwrite);
    outfile->Close();
    reader.Close();
    delete outfile;
  }

  std::cout << "Output saved to: " << output_filename << std::endl;

  return global_header;
}

std::pair<std::vector<RawHit>, UShort_t>
InitUtils::ConvertCoMPASSBinToHits(const TString input_filename,
                                   UShort_t global_header_override,
                                   Bool_t skip_bad_events) {
  std::vector<RawHit> hits;

  if (gSystem->AccessPathName(input_filename)) {
    std::cout << "ERROR: Input file does not exist: " << input_filename
              << std::endl;
    return std::make_pair(hits, static_cast<UShort_t>(0));
  }

  CoMPASSReader reader;
  Bool_t open_success =
      (global_header_override != 0)
          ? reader.Open(input_filename.Data(), global_header_override)
          : reader.Open(input_filename.Data());

  if (!open_success) {
    std::cout << "ERROR: Failed to open CoMPASS binary file" << std::endl;
    return std::make_pair(hits, static_cast<UShort_t>(0));
  }

  UShort_t global_header = reader.GetGlobalHeader();
  Bool_t has_energy_ch = (global_header & 0x0001);
  Bool_t has_waveform = (global_header & 0x0008);

  if (!has_energy_ch) {
    std::cout << "WARNING: File has no channel-energy field; RawHit.energy "
                 "will be zero for every hit."
              << std::endl;
  }

  if (!has_waveform) {
    // Fast bulk-read path: single fread of the data section, then memcpy-parse.
    std::cout
        << "Using bulk-read fast path (fixed-stride records, no waveform)."
        << std::endl;
    reader.Close();
    return BulkReadCoMPASSHits(input_filename, global_header,
                               global_header_override, skip_bad_events);
  }

  Long64_t event_count = 0;
  Long64_t warning_fake = 0;
  Long64_t warning_saturated = 0;
  Long64_t warning_pileup = 0;
  Long64_t warning_memory_full = 0;
  Long64_t warning_trigger_lost = 0;
  Long64_t warning_pll_loss = 0;
  Long64_t warning_over_temp = 0;
  Long64_t warning_adc_shutdown = 0;

  std::cout << "Reading events..." << std::endl;
  if (skip_bad_events) {
    std::cout
        << "Filtering enabled: skipping fake, saturated, and pileup events"
        << std::endl;
  }

  while (reader.ReadEvent()) {
    const CoMPASSData &event = reader.GetCurrentEvent();
    if (event_count == 0) {
      event.PrintHeader();
    }

    if (event.isFakeEvent()) {
      warning_fake++;
      if (skip_bad_events)
        continue;
    }
    if (event.isInputSaturating() || event.hasSaturation()) {
      warning_saturated++;
      if (skip_bad_events)
        continue;
    }
    if (event.isPileup()) {
      warning_pileup++;
      if (skip_bad_events)
        continue;
    }

    if (event.hasMemoryFull())
      warning_memory_full++;
    if (event.hasTriggerLost())
      warning_trigger_lost++;
    if (event.hasPLLLockLoss())
      warning_pll_loss++;
    if (event.isOverTemperature())
      warning_over_temp++;
    if (event.isADCShutdown())
      warning_adc_shutdown++;

    RawHit hit;
    hit.board = event.board;
    hit.channel = event.channel;
    hit.energy = event.energy_ch;
    hit.timestamp = event.timestamp;
    hit.flags = event.flags;
    hits.push_back(hit);

    event_count++;
  }

  std::cout << "Conversion complete." << std::endl;
  std::cout << "Total events processed: " << event_count << std::endl;

  if (warning_fake > 0 || warning_saturated > 0 || warning_pileup > 0) {
    std::cout << "Events with rejection-quality flags:" << std::endl;
    if (warning_fake > 0) {
      std::cout << "  Fake events: " << warning_fake;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_saturated > 0) {
      std::cout << "  Saturated: " << warning_saturated;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_pileup > 0) {
      std::cout << "  Pileup: " << warning_pileup;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  if (warning_memory_full > 0) {
    std::cout << "WARNING: " << warning_memory_full
              << " events with memory full flag" << std::endl;
  }
  if (warning_trigger_lost > 0) {
    std::cout << "WARNING: " << warning_trigger_lost
              << " events with trigger lost flag" << std::endl;
  }
  if (warning_pll_loss > 0) {
    std::cout << "WARNING: " << warning_pll_loss << " events with PLL lock loss"
              << std::endl;
  }
  if (warning_over_temp > 0) {
    std::cout << "WARNING: " << warning_over_temp
              << " events with over temperature" << std::endl;
  }
  if (warning_adc_shutdown > 0) {
    std::cout << "WARNING: " << warning_adc_shutdown
              << " events with ADC shutdown" << std::endl;
  }

  std::cout << "Total bytes read: " << reader.GetBytesRead() << std::endl;

  reader.Close();

  return std::make_pair(hits, global_header);
}

Bool_t InitUtils::ConvertWavedumpBinToROOT(const TString input_filename,
                                           const TString output_name,
                                           Bool_t corrections_enabled) {
  const TString base_dir = IO::GetRootFilesBaseDir();
  if (gSystem->AccessPathName(base_dir)) {
    gSystem->mkdir(base_dir, kTRUE);
  }

  if (gSystem->AccessPathName(input_filename)) {
    std::cout << "ERROR: Input file does not exist: " << input_filename
              << std::endl;
    return kFALSE;
  }

  TString output_subpath = output_name + "_raw.root";
  TString output_filename = base_dir + "/" + output_subpath;

  WaveDump742Reader reader(corrections_enabled);

  if (!reader.Open(input_filename.Data())) {
    std::cout << "ERROR: Failed to open WaveDump binary file" << std::endl;
    return kFALSE;
  }

  std::cout << "Corrections: " << (corrections_enabled ? "enabled" : "disabled")
            << std::endl;

  TFile *outfile = nullptr;
  TTree *tree = nullptr;
  UInt_t channel_br = 0, event_counter = 0, trigger_time_tag = 0;
  TArrayS *samples = nullptr;

  {
    IO::ScopedRootLock setup_guard;

    outfile = IO::OpenForWriting(output_subpath);
    if (!outfile || outfile->IsZombie()) {
      std::cout << "ERROR: Could not create output file " << output_filename
                << std::endl;
      reader.Close();
      return kFALSE;
    }

    tree = new TTree("Data_R", "WaveDump 742 Binary Data");

    tree->Branch("Channel", &channel_br, "Channel/i");
    tree->Branch("EventCounter", &event_counter, "EventCounter/i");
    tree->Branch("TriggerTimeTag", &trigger_time_tag, "TriggerTimeTag/i");
    samples = new TArrayS();
    tree->Branch("Samples", &samples);
  }

  Long64_t event_count = 0;

  std::cout << "Reading events..." << std::endl;

  while (reader.ReadEvent()) {
    const WaveDump742Data &event = reader.GetCurrentEvent();

    channel_br = event.channel;
    event_counter = event.event_counter;
    trigger_time_tag = event.group_trigger_time_tag;
    *samples = event.samples;

    tree->Fill();
    event_count++;
  }

  std::cout << "Conversion complete." << std::endl;
  std::cout << "Total events processed: " << event_count << std::endl;
  std::cout << "Samples per event: "
            << (event_count > 0 ? samples->GetSize() : 0) << std::endl;
  std::cout << "Total bytes read: " << reader.GetBytesRead() << std::endl;

  {
    IO::ScopedRootLock teardown_guard;
    outfile->cd();
    tree->Write("", TObject::kOverwrite);
    outfile->Close();
    reader.Close();
    delete outfile;
  }

  std::cout << "Output saved to: " << output_filename << std::endl;

  return kTRUE;
}
