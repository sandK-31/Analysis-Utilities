#ifndef INITUTILS_H
#define INITUTILS_H

#include "BinaryUtils.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstdio>
#include <cstring>
#include <utility>

class InitUtils {
public:
  static void
  SetROOTPreferences(PlotSaveFormat save_format = PlotSaveFormat::kPNG,
                     const TString &plots_dir = "",
                     const TString &root_files_dir = "",
                     Bool_t enable_mt = kTRUE);
  static Bool_t ConvertWavedumpBinToROOT(const TString input_filename,
                                         const TString output_name,
                                         Bool_t corrections_enabled = kTRUE);
  static UShort_t ConvertCoMPASSBinToROOT(const TString input_filename,
                                          const TString output_name,
                                          UShort_t global_header_override,
                                          Bool_t skip_bad_events = kFALSE);
  static std::pair<std::vector<RawHit>, UShort_t>
  ConvertCoMPASSBinToHits(const TString input_filename,
                          UShort_t global_header_override = 0,
                          Bool_t skip_bad_events = kFALSE);
  static Bool_t ConvertCoMPASSCSVToROOT();
};

#endif
