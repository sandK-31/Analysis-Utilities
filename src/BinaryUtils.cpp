#include "BinaryUtils.hpp"

CoMPASSData::CoMPASSData()
    : header(0), board(0), channel(0), timestamp(0), energy_ch(0),
      energy_cal(0.0), energy_short_ch(0), flags(0), waveform_code(0),
      num_samples(0) {
  std::cout << "This version of CoMPASS binary conversion is based on manual "
               "revision 25.1..."
            << std::endl;
}

TString CoMPASSData::getWaveformCodeName() const {
  switch (waveform_code) {
  case INPUT:
    return "Input";
  case RC_CR:
    return "RC-CR (DPP-PHA)";
  case RC_CR2:
    return "RC-CR2 (DPP-PHA)";
  case TRAPEZOID:
    return "Trapezoid (DPP-PHA)";
  case BASELINE:
    return "Baseline";
  case THRESHOLD:
    return "Threshold";
  case CFD:
    return "CFD (DPP-PSD)";
  case TRAPEZOID_BASELINE:
    return "Trapezoid-Baseline (DPP-PHA)";
  case FAST_TRIANGLE:
    return "Fast Triangle (x27xx DPP-PHA)";
  case SMOOTHED_INPUT:
    return "Smoothed Input (DPP-PSD)";
  default:
    return "Unknown";
  }
}

std::vector<TString> CoMPASSData::getActiveFlags() const {
  std::vector<TString> active_flags;

  if (hasDeadtime())
    active_flags.push_back("Deadtime");
  if (hasTimestampRollover())
    active_flags.push_back("Timestamp Rollover");
  if (hasTimestampResetExt())
    active_flags.push_back("Timestamp Reset (Ext)");
  if (isFakeEvent())
    active_flags.push_back("Fake Event");
  if (hasMemoryFull())
    active_flags.push_back("Memory Full");
  if (hasTriggerLost())
    active_flags.push_back("Trigger Lost");
  if (hasNTriggersLost())
    active_flags.push_back("N Triggers Lost");
  if (hasSaturation())
    active_flags.push_back("Saturation");
  if (has1024Triggers())
    active_flags.push_back("1024 Triggers");
  if (isFirstAfterBusy())
    active_flags.push_back("First After Busy");
  if (isInputSaturating())
    active_flags.push_back("Input Saturating");
  if (hasNTriggersCounted())
    active_flags.push_back("N Triggers Counted");
  if (isNotMatchedTimeFilter())
    active_flags.push_back("Not Matched Time Filter");
  if (hasFineTimestamp())
    active_flags.push_back("Fine Timestamp");
  if (isPileup())
    active_flags.push_back("Pile-up");
  if (hasPLLLockLoss())
    active_flags.push_back("PLL Lock Loss");
  if (isOverTemperature())
    active_flags.push_back("Over Temperature");
  if (isADCShutdown())
    active_flags.push_back("ADC Shutdown");

  return active_flags;
}

void CoMPASSData::PrintHeader() const {
  std::cout << "CoMPASS event header..." << std::endl;
  std::cout << "Header:    0x" << std::hex << header << std::dec
            << " (Binary: " << std::bitset<16>(header) << ")" << std::endl;

  std::cout << std::endl;
  std::cout << "Control bits..." << std::endl;
  std::cout << "Energy (ch):    " << (hasEnergyCh() ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (cal):   " << (hasEnergyCal() ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (short): " << (hasEnergyShort() ? "YES" : "NO")
            << std::endl;
  std::cout << "Waveform:       " << (hasWaveform() ? "YES" : "NO")
            << std::endl;
}

void CoMPASSData::PrintFlags() const {
  std::cout << "Flags (0x" << std::hex << flags << std::dec << ")" << std::endl;
  std::cout << "Binary: " << std::bitset<32>(flags) << std::endl;

  std::vector<TString> active = getActiveFlags();
  if (active.empty()) {
    std::cout << "No flags set" << std::endl;
  } else {
    std::cout << "Active flags:" << std::endl;
    Int_t n_flags = active.size();
    for (Int_t flag = 0; flag < n_flags; flag++) {
      std::cout << "  - " << active.at(flag).Data() << std::endl;
    }
  }
}

void CoMPASSData::PrintWaveform() const {
  if (hasWaveform()) {
    std::cout << "Waveform..." << std::endl;
    std::cout << "Code:    " << static_cast<Int_t>(waveform_code) << " ("
              << getWaveformCodeName().Data() << ")" << std::endl;
    std::cout << "Samples: " << num_samples << std::endl;
    if (samples.GetSize() > 0) {
      std::cout << "First 5 samples: ";
      for (Int_t i = 0; i < TMath::Min(5, samples.GetSize()); i++) {
        std::cout << samples[i] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void CoMPASSData::Print() const {
  PrintHeader();

  if (hasEnergyCh()) {
    std::cout << "Energy (ch):    " << energy_ch << std::endl;
  }
  if (hasEnergyCal()) {
    std::cout << "Energy (cal):   " << energy_cal << " keV/MeV" << std::endl;
  }
  if (hasEnergyShort()) {
    std::cout << "Energy (short): " << energy_short_ch << std::endl;
  }

  PrintFlags();
  PrintWaveform();
}

Bool_t CoMPASSReader::Open(const char *fname) {
  if (!BinaryReader::Open(fname)) {
    return kFALSE;
  }

  file.read(reinterpret_cast<char *>(&global_header), sizeof(UShort_t));
  bytes_read += sizeof(UShort_t);

  if ((global_header & 0xCAE0) != 0xCAE0) {
    std::cerr << "WARNING: Header does not match 0xCAEx pattern: 0x" << std::hex
              << global_header << std::dec << std::endl;
  }

  std::cout << "CoMPASS file opened: " << fname << std::endl;

  return kTRUE;
}

Bool_t CoMPASSReader::Open(const char *fname, UShort_t header_override) {
  if (!BinaryReader::Open(fname)) {
    return kFALSE;
  }

  if (header_override != 0) {
    global_header = header_override;
    std::cout << "CoMPASS continuation file opened: " << fname << std::endl;
    std::cout << "Using provided global header: 0x" << std::hex << global_header
              << std::dec << std::endl;
    std::cout << "Control bits: " << std::bitset<4>(global_header & 0x000F)
              << std::endl;
  } else {
    file.read(reinterpret_cast<char *>(&global_header), sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);

    if ((global_header & 0xCAE0) != 0xCAE0) {
      std::cerr << "WARNING: Header does not match 0xCAEx pattern: 0x"
                << std::hex << global_header << std::dec << std::endl;
    }

    std::cout << "CoMPASS file opened: " << fname << std::endl;
    std::cout << "Global header: 0x" << std::hex << global_header << std::dec
              << std::endl;
    std::cout << "Control bits: " << std::bitset<4>(global_header & 0x000F)
              << std::endl;
  }

  return kTRUE;
}

Bool_t CoMPASSReader::ReadEvent() {
  if (IsEOF()) {
    return kFALSE;
  }

  current_event.header = global_header;

  file.read(reinterpret_cast<char *>(&current_event.board), sizeof(UShort_t));
  file.read(reinterpret_cast<char *>(&current_event.channel), sizeof(UShort_t));
  file.read(reinterpret_cast<char *>(&current_event.timestamp),
            sizeof(ULong64_t));
  bytes_read += sizeof(UShort_t) * 2 + sizeof(ULong64_t);

  if (current_event.hasEnergyCh()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_ch),
              sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);
  }

  if (current_event.hasEnergyCal()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_cal),
              sizeof(Double_t));
    bytes_read += sizeof(Double_t);
  }

  if (current_event.hasEnergyShort()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_short_ch),
              sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);
  }

  file.read(reinterpret_cast<char *>(&current_event.flags), sizeof(UInt_t));
  bytes_read += sizeof(UInt_t);

  if (current_event.hasWaveform()) {
    file.read(reinterpret_cast<char *>(&current_event.waveform_code),
              sizeof(UChar_t));
    file.read(reinterpret_cast<char *>(&current_event.num_samples),
              sizeof(UInt_t));
    bytes_read += sizeof(UChar_t) + sizeof(UInt_t);

    // Sanity-cap num_samples; a corrupt header can read garbage and trigger
    // a multi-GB allocation in TArrayS::Set.
    const UInt_t kMaxSamples = 1u << 20;
    if (current_event.num_samples > kMaxSamples) {
      std::cerr << "ERROR: implausible num_samples "
                << current_event.num_samples << " at byte " << bytes_read
                << " (treating as truncated)" << std::endl;
      return kFALSE;
    }

    current_event.samples.Set(current_event.num_samples);
    std::vector<UShort_t> sample_buf(current_event.num_samples);
    file.read(reinterpret_cast<char *>(sample_buf.data()),
              current_event.num_samples * sizeof(UShort_t));
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete waveform at byte " << bytes_read
                << " (truncated file, event discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < current_event.num_samples; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(sample_buf[i]), i);
    }
    bytes_read += current_event.num_samples * sizeof(UShort_t);
  }

  return kTRUE;
}

WaveDump742Data::WaveDump742Data()
    : event_size(0), board_id(0), pattern(0), channel(0), event_counter(0),
      group_trigger_time_tag(0), dc_offset(0), start_index_cell(0) {
  std::cout << "This version of wavedump binary conversion for 742 family "
               "digitizers is based on manual "
               "revision 21"
            << std::endl;
}

void WaveDump742Data::Print() const {
  std::cout << "WaveDump 742 Event..." << std::endl;
  std::cout << "Event size:        " << event_size << std::endl;
  std::cout << "Board ID:          " << board_id << std::endl;
  std::cout << "Channel:           " << channel << std::endl;
  std::cout << "Event counter:     " << event_counter << std::endl;
  std::cout << "Group trigger tag: " << group_trigger_time_tag << std::endl;
  std::cout << "DC offset:         " << dc_offset << std::endl;
  std::cout << "Start cell:        " << start_index_cell << std::endl;
  std::cout << "Samples:           " << samples.GetSize() << std::endl;
}

Bool_t WaveDump742Reader::ReadEvent() {
  if (IsEOF()) {
    return kFALSE;
  }

  UInt_t headers[8];
  file.read(reinterpret_cast<char *>(headers), 8 * sizeof(UInt_t));
  if (file.fail()) {
    if (file.gcount() > 0) {
      std::cerr << "WARNING: Incomplete event header at byte " << bytes_read
                << " (" << file.gcount() << " of " << 8 * sizeof(UInt_t)
                << " header bytes read, truncated file)" << std::endl;
    }
    return kFALSE;
  }
  bytes_read += 8 * sizeof(UInt_t);

  current_event.event_size = headers[0];
  current_event.board_id = headers[1];
  current_event.pattern = headers[2];
  current_event.channel = headers[3];
  current_event.event_counter = headers[4];
  current_event.group_trigger_time_tag = headers[5];
  current_event.dc_offset = headers[6];
  current_event.start_index_cell = headers[7];

  // event_size is in bytes and includes the 8-word (32-byte) header
  UInt_t header_bytes = 8 * sizeof(UInt_t);
  if (current_event.event_size <= header_bytes) {
    std::cerr << "ERROR: Invalid event_size " << current_event.event_size
              << " at event " << current_event.event_counter << std::endl;
    return kFALSE;
  }

  UInt_t sample_bytes = current_event.event_size - header_bytes;
  UInt_t sample_size = sample_bytes / sizeof(UInt_t);

  current_event.samples.Set(sample_size);

  if (corrections_enabled) {
    std::vector<Float_t> float_samples(sample_size);
    file.read(reinterpret_cast<char *>(float_samples.data()), sample_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete event at byte " << bytes_read
                << " (truncated file, event " << current_event.event_counter
                << " discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < sample_size; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(float_samples[i]), i);
    }
  } else {
    std::vector<UInt_t> int_samples(sample_size);
    file.read(reinterpret_cast<char *>(int_samples.data()), sample_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete event at byte " << bytes_read
                << " (truncated file, event " << current_event.event_counter
                << " discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < sample_size; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(int_samples[i] & 0xFFF),
                                  i);
    }
  }

  bytes_read += sample_bytes;

  return kTRUE;
}
