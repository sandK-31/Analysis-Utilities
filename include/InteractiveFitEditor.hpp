#ifndef INTERACTIVEFITEDITOR_H
#define INTERACTIVEFITEDITOR_H

#include "FittingUtils.hpp"
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
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TTimer.h>
#include <iostream>

class InteractiveFitEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 10000;
  static const Int_t kNDrawPts = 500;

  static const Int_t kBtnRefit = 1000;
  static const Int_t kBtnAccept = 1001;
  static const Int_t kBtnCancel = 1002;
  static const Int_t kBtnReset = 1003;
  static const Int_t kSliderBase = 2000;
  static const Int_t kEntryBase = 3000;
  static const Int_t kFixBase = 4000;
  static const Int_t kLoBoundBase = 5000;
  static const Int_t kHiBoundBase = 6000;
  static const Int_t kRangeSlider = 7000;
  static const Int_t kRangeLoEntry = 7001;
  static const Int_t kRangeHiEntry = 7002;

  TH1 *hist_;
  TF1 *fit_func_;
  TString info_label_text_;
  Double_t range_low_;
  Double_t range_high_;
  Double_t original_range_low_;
  Double_t original_range_high_;
  Double_t hist_x_min_;
  Double_t hist_x_max_;
  Int_t num_peaks_;
  Int_t num_params_;

  Double_t *original_params_;
  Double_t *original_bounds_low_;
  Double_t *original_bounds_high_;
  Bool_t *original_fixed_;

  Double_t *current_bounds_low_;
  Double_t *current_bounds_high_;

  TRootEmbeddedCanvas *embedded_canvas_;
  TPad *main_pad_;
  TPad *residual_pad_;

  TGHSlider **sliders_;
  TGNumberEntry **value_entries_;
  TGCheckButton **fix_checks_;
  TGNumberEntry **lo_bound_entries_;
  TGNumberEntry **hi_bound_entries_;

  TGDoubleHSlider *range_slider_;
  TGNumberEntry *range_lo_entry_;
  TGNumberEntry *range_hi_entry_;

  TH1 *hist_draw_;
  TGraph *comp_graphs_[3][4];
  TF1 *bkg_draw_;
  TGraph *res_graph_;
  TF1 *zero_line_;
  TF1 *plus3_line_;
  TF1 *minus3_line_;
  TLatex *chi2_label_;
  Int_t n_res_points_;

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t done_;
  Bool_t syncing_;
  TTimer *redraw_timer_;

  void BuildGUI();
  void BuildPeakTab(TGCompositeFrame *parent, Int_t peak_idx);
  void BuildBackgroundTab(TGCompositeFrame *parent);
  void AddParamRow(TGCompositeFrame *parent, Int_t param_idx, const char *name);

  void InitDrawing();
  void UpdateCanvas();
  void UpdateCompPoints();
  void UpdateResPoints();

  void SyncAllWidgets();
  void SyncWidget(Int_t param_idx);

  void OnSliderMoved(Int_t param_idx);
  void OnEntryChanged(Int_t param_idx);
  void OnBoundsChanged(Int_t param_idx);
  void OnFixToggled(Int_t param_idx);
  void OnRangeChanged();

  void DoRefit();
  void DoAccept();
  void DoCancel();
  void DoReset();

  Int_t ValToSlider(Int_t param_idx, Double_t val);
  Double_t SliderToVal(Int_t param_idx, Int_t pos);
  Bool_t IsFixed(Int_t param_idx);
  void GetDefaultBounds(Int_t param_idx, Double_t &lo, Double_t &hi);
  Int_t BkgConstIdx() { return num_peaks_ * 10; }
  Int_t BkgSlopeIdx() { return num_peaks_ * 10 + 1; }
  static Int_t PeakStyle(Int_t peak_idx);

public:
  InteractiveFitEditor(const TGWindow *parent, TH1 *hist, TF1 *fit_func,
                       Double_t range_low, Double_t range_high, Int_t num_peaks,
                       const TString &info_label = "");
  virtual ~InteractiveFitEditor();

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  virtual Bool_t HandleTimer(TTimer *timer);
  virtual void CloseWindow();

  Bool_t WasAccepted() const { return accepted_; }
  Bool_t IsDone() const { return done_; }
  TTimer *GetRedrawTimer() { return redraw_timer_; }
  TRootEmbeddedCanvas *GetEmbeddedCanvas() { return embedded_canvas_; }
};

Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func, Double_t range_low,
                                  Double_t range_high, Int_t num_peaks,
                                  const TString &info_label);

#endif
