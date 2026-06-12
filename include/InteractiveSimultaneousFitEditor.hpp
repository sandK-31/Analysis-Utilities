#ifndef INTERACTIVESIMULTANEOUSFITEDITOR_H
#define INTERACTIVESIMULTANEOUSFITEDITOR_H

#include "RooFitUtils.hpp"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <TCanvas.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGDoubleSlider.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGSlider.h>
#include <TGTab.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TTimer.h>
#include <map>
#include <vector>

struct SimEditorChannelView {
  TString name;
  TH1 *hist;
  const std::vector<Double_t> *events;
  Float_t display_bin_width_kev;
  RooAbsPdf *pdf;
  RooAbsData *data;
  std::vector<RooFitPeakModel> *peaks;
  RooFitBackgroundModel *bkg;
  Int_t num_peaks;

  std::vector<RooRealVar *> params;
  std::vector<Double_t> original_params;
  std::vector<Double_t> original_bounds_low;
  std::vector<Double_t> original_bounds_high;
  std::vector<Bool_t> original_fixed;
  std::vector<Double_t> current_bounds_low;
  std::vector<Double_t> current_bounds_high;

  std::vector<TGHSlider *> sliders;
  std::vector<TGNumberEntry *> value_entries;
  std::vector<TGCheckButton *> fix_checks;
  std::vector<TGNumberEntry *> lo_bound_entries;
  std::vector<TGNumberEntry *> hi_bound_entries;

  TRootEmbeddedCanvas *embedded_canvas;
  TPad *main_pad;
  TPad *residual_pad;
  TH1 *hist_draw;
  TGraph *total_graph;
  TGraph *bkg_graph;
  TGraph *comp_graphs[3][4];
  TGraph *res_graph;
  TF1 *zero_line;
  TF1 *plus3_line;
  TF1 *minus3_line;
  TLatex *chi2_label;
  Int_t n_res_points;
};

class InteractiveSimultaneousFitEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 10000;
  static const Int_t kNDrawPts = 500;

  static const Int_t kBtnRefit = 1000;
  static const Int_t kBtnAccept = 1001;
  static const Int_t kBtnCancel = 1002;
  static const Int_t kBtnReset = 1003;

  static const Int_t kRangeSlider = 7000;
  static const Int_t kRangeLoEntry = 7001;
  static const Int_t kRangeHiEntry = 7002;

  static const Int_t kSliderStride = 100;
  static const Int_t kSliderBase = 10000;
  static const Int_t kEntryBase = 20000;
  static const Int_t kFixBase = 30000;
  static const Int_t kLoBoundBase = 40000;
  static const Int_t kHiBoundBase = 50000;

  RooSimultaneous *sim_pdf_;
  RooAbsData *combined_data_;
  RooRealVar *x_;
  std::vector<SimEditorChannelView> channels_;
  TString info_label_text_;
  Double_t range_low_;
  Double_t range_high_;
  Double_t original_range_low_;
  Double_t original_range_high_;
  Double_t hist_x_min_;
  Double_t hist_x_max_;

  TGDoubleHSlider *range_slider_;
  TGNumberEntry *range_lo_entry_;
  TGNumberEntry *range_hi_entry_;

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t done_;
  Bool_t syncing_;
  Bool_t fit_debug_;
  TTimer *redraw_timer_;

  Int_t ChannelIndexFromWidgetId(Int_t base, Int_t parm1, Int_t &local_idx);

  void BuildGUI();
  void BuildChannelTab(TGCompositeFrame *parent, Int_t ch_idx);
  void BuildPeakSubTab(TGCompositeFrame *parent, Int_t ch_idx, Int_t peak_idx);
  void BuildBackgroundSubTab(TGCompositeFrame *parent, Int_t ch_idx);
  void AddParamRow(TGCompositeFrame *parent, Int_t ch_idx, Int_t param_idx,
                   const char *name);

  void InitDrawing();
  void InitChannelDrawing(SimEditorChannelView &cv);
  void UpdateCanvases();
  void UpdateChannelGraphs(SimEditorChannelView &cv);
  void UpdateChannelResiduals(SimEditorChannelView &cv);
  void UpdateChannelChi2(SimEditorChannelView &cv);

  void SyncAllWidgets();
  void SyncChannelWidget(Int_t ch_idx, Int_t param_idx);

  void OnSliderMoved(Int_t ch_idx, Int_t param_idx);
  void OnEntryChanged(Int_t ch_idx, Int_t param_idx);
  void OnBoundsChanged(Int_t ch_idx, Int_t param_idx);
  void OnFixToggled(Int_t ch_idx, Int_t param_idx);
  void OnRangeChanged();

  void DoRefit();
  void DoAccept();
  void DoCancel();
  void DoReset();

  // Re-tie each channel's linear-background slope lower bound to the current
  // shared fit range so 1+slope*x stays positive across the whole window.
  void ApplyBackgroundSlopeBounds();

  // Debug aid: scan each channel's components for a non-finite normalization
  // integral or a non-finite raw value across the range, to name the pdf
  // behind an invalid NLL. Only called when fit_debug_ is set.
  void DiagnoseInvalidComponents();

  Int_t ValToSlider(Int_t ch_idx, Int_t param_idx, Double_t val);
  Double_t SliderToVal(Int_t ch_idx, Int_t param_idx, Int_t pos);
  Bool_t IsFixed(Int_t ch_idx, Int_t param_idx);
  void GetDefaultBounds(Int_t ch_idx, Int_t param_idx, Double_t &lo,
                        Double_t &hi);
  Int_t BkgConstIdx(Int_t ch_idx) { return channels_[ch_idx].num_peaks * 10; }
  Int_t BkgSlopeIdx(Int_t ch_idx) {
    return channels_[ch_idx].num_peaks * 10 + 1;
  }
  static Int_t PeakStyle(Int_t peak_idx);

public:
  InteractiveSimultaneousFitEditor(
      const TGWindow *parent, RooSimultaneous *sim_pdf,
      RooAbsData *combined_data, RooRealVar *x,
      const std::vector<SimEditorChannelView> &channel_views,
      Double_t range_low, Double_t range_high, const TString &info_label = "",
      Bool_t fit_debug = kFALSE);
  virtual ~InteractiveSimultaneousFitEditor();

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  virtual Bool_t HandleTimer(TTimer *timer);
  virtual void CloseWindow();

  Bool_t WasAccepted() const { return accepted_; }
  Bool_t IsDone() const { return done_; }
  TTimer *GetRedrawTimer() { return redraw_timer_; }
  std::vector<SimEditorChannelView> *GetChannels() { return &channels_; }
};

Bool_t LaunchInteractiveSimultaneousFitEditor(
    RooSimultaneous *sim_pdf, RooAbsData *combined_data, RooRealVar *x,
    std::vector<SimEditorChannelView> &channel_views, Double_t range_low,
    Double_t range_high, const TString &info_label, Bool_t fit_debug = kFALSE);

#endif
