#include "InteractiveRooFitEditor.hpp"

#include "InteractiveEditorX11Guard.hpp"
#include "PlottingUtils.hpp"

#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <TROOT.h>
#include <iostream>

InteractiveRooFitEditor::InteractiveRooFitEditor(
    const TGWindow *parent, TH1 *hist, const std::vector<Double_t> *events,
    Float_t display_bin_width_kev, RooAbsPdf *total_pdf, RooRealVar *x,
    RooAbsData *data, std::vector<RooFitPeakModel> *peaks,
    RooFitBackgroundModel *bkg, Double_t range_low, Double_t range_high,
    const TString &info_label)
    : TGMainFrame(parent, 1400, 900) {

  hist_ = hist;
  events_ = events;
  display_bin_width_kev_ = display_bin_width_kev;
  total_pdf_ = total_pdf;
  x_ = x;
  data_ = data;
  peaks_ = peaks;
  bkg_ = bkg;
  info_label_text_ = info_label;
  range_low_ = range_low;
  range_high_ = range_high;
  original_range_low_ = range_low;
  original_range_high_ = range_high;
  hist_x_min_ = hist->GetXaxis()->GetXmin();
  hist_x_max_ = hist->GetXaxis()->GetXmax();
  num_peaks_ = (Int_t)peaks->size();
  num_params_ = num_peaks_ * 10 + 2;

  accepted_ = kFALSE;
  done_ = kFALSE;
  needs_redraw_ = kFALSE;
  syncing_ = kFALSE;

  for (Int_t pi = 0; pi < num_peaks_; pi++) {
    RooFitPeakModel &p = (*peaks_)[pi];
    params_.push_back(p.mu);
    params_.push_back(p.sigma);
    params_.push_back(p.gaus_yield);
    params_.push_back(p.ratio_step);
    params_.push_back(p.ratio_low_exp);
    params_.push_back(p.tau_low_exp);
    params_.push_back(p.ratio_low_lin);
    params_.push_back(p.slope_low_lin);
    params_.push_back(p.ratio_high_exp);
    params_.push_back(p.tau_high_exp);
  }
  params_.push_back(bkg_->bkg_yield);
  params_.push_back(bkg_->bkg_slope);

  original_params_.resize(num_params_);
  original_bounds_low_.resize(num_params_);
  original_bounds_high_.resize(num_params_);
  original_fixed_.resize(num_params_);
  current_bounds_low_.resize(num_params_);
  current_bounds_high_.resize(num_params_);

  for (Int_t i = 0; i < num_params_; i++) {
    original_params_[i] = params_[i]->getVal();
    original_fixed_[i] = params_[i]->isConstant();

    Double_t lo, hi;
    if (params_[i]->hasMin() && params_[i]->hasMax()) {
      lo = params_[i]->getMin();
      hi = params_[i]->getMax();
    } else {
      GetDefaultBounds(i, lo, hi);
    }
    original_bounds_low_[i] = lo;
    original_bounds_high_[i] = hi;
    current_bounds_low_[i] = lo;
    current_bounds_high_[i] = hi;
  }

  sliders_.resize(num_params_);
  value_entries_.resize(num_params_);
  fix_checks_.resize(num_params_);
  lo_bound_entries_.resize(num_params_);
  hi_bound_entries_.resize(num_params_);

  hist_draw_ = nullptr;
  total_graph_ = nullptr;
  bkg_graph_ = nullptr;
  res_graph_ = nullptr;
  zero_line_ = nullptr;
  plus3_line_ = nullptr;
  minus3_line_ = nullptr;
  chi2_label_ = nullptr;
  n_res_points_ = 0;
  for (Int_t p = 0; p < 3; p++) {
    for (Int_t c = 0; c < 4; c++) {
      comp_graphs_[p][c] = nullptr;
    }
  }

  BuildGUI();
  InitDrawing();

  redraw_timer_ = new TTimer(this, 50);
  redraw_timer_->TurnOn();

  SetWindowName("Interactive RooFit Editor");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  UpdateCanvas();
}

InteractiveRooFitEditor::~InteractiveRooFitEditor() {
  if (redraw_timer_) {
    redraw_timer_->TurnOff();
    delete redraw_timer_;
  }

  // Drawing primitives drawn on the embedded canvas — TPad does not
  // auto-delete user primitives, so they have to go by hand.
  //
  // Pre-clear each pad's primitives list with "nodelete" before destroying
  // the primitives themselves. The launcher removes the embedded canvas
  // from gROOT->GetListOfCanvases() to suppress X11 errors during teardown,
  // which also disables RecursiveRemove's cleanup of pad lists; without
  // this manual Clear, ~TPad later walks stale pointers and segfaults.
  if (main_pad_ && main_pad_->GetListOfPrimitives()) {
    main_pad_->GetListOfPrimitives()->Clear("nodelete");
  }
  if (residual_pad_ && residual_pad_->GetListOfPrimitives()) {
    residual_pad_->GetListOfPrimitives()->Clear("nodelete");
  }
  for (Int_t p = 0; p < 3; p++) {
    for (Int_t c = 0; c < 4; c++) {
      delete comp_graphs_[p][c];
    }
  }
  delete chi2_label_;
  delete zero_line_;
  delete plus3_line_;
  delete minus3_line_;
  delete res_graph_;
  delete bkg_graph_;
  delete total_graph_;
  delete hist_draw_;
  delete residual_pad_;
  delete main_pad_;
}

void InteractiveRooFitEditor::BuildGUI() {
  TGHorizontalFrame *main_frame = new TGHorizontalFrame(this, 1400, 850);
  AddFrame(main_frame,
           new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  embedded_canvas_ =
      new TRootEmbeddedCanvas("RooFitEditorCanvas", main_frame, 850, 800);
  main_frame->AddFrame(
      embedded_canvas_,
      new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  TGVerticalFrame *controls = new TGVerticalFrame(main_frame, 520, 800);
  main_frame->AddFrame(
      controls, new TGLayoutHints(kLHintsExpandY | kLHintsRight, 2, 2, 2, 2));

  if (info_label_text_.Length() > 0) {
    TGGroupFrame *info_grp = new TGGroupFrame(controls, "Info", kVerticalFrame);
    controls->AddFrame(info_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    TGLabel *info_label = new TGLabel(info_grp, info_label_text_);
    info_grp->AddFrame(info_label,
                       new TGLayoutHints(kLHintsCenterX, 4, 4, 4, 4));
  }

  TGGroupFrame *range_grp =
      new TGGroupFrame(controls, "Fit Range", kVerticalFrame);
  controls->AddFrame(range_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));

  range_slider_ =
      new TGDoubleHSlider(range_grp, 200, kDoubleScaleNo, kRangeSlider);
  range_slider_->SetRange(hist_x_min_, hist_x_max_);
  range_slider_->SetPosition(range_low_, range_high_);
  range_slider_->Associate(this);
  range_grp->AddFrame(range_slider_,
                      new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

  TGHorizontalFrame *range_entries = new TGHorizontalFrame(range_grp, 200, 28);
  range_grp->AddFrame(range_entries,
                      new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

  TGLabel *range_lo_label = new TGLabel(range_entries, "Low:");
  range_entries->AddFrame(range_lo_label,
                          new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  range_lo_entry_ = new TGNumberEntry(
      range_entries, range_low_, 8, kRangeLoEntry, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  range_lo_entry_->GetNumberEntry()->Associate(this);
  range_entries->AddFrame(
      range_lo_entry_,
      new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));

  TGLabel *range_hi_label = new TGLabel(range_entries, "High:");
  range_entries->AddFrame(range_hi_label,
                          new TGLayoutHints(kLHintsCenterY, 8, 2, 2, 2));
  range_hi_entry_ = new TGNumberEntry(
      range_entries, range_high_, 8, kRangeHiEntry, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  range_hi_entry_->GetNumberEntry()->Associate(this);
  range_entries->AddFrame(
      range_hi_entry_,
      new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));

  TGTab *tabs = new TGTab(controls, 510, 700);
  controls->AddFrame(
      tabs, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  for (Int_t p = 0; p < num_peaks_; p++) {
    TString tab_name = TString::Format("Peak %d", p + 1);
    TGCompositeFrame *tab_frame = tabs->AddTab(tab_name);
    BuildPeakTab(tab_frame, p);
  }

  TGCompositeFrame *bkg_tab = tabs->AddTab("Background");
  BuildBackgroundTab(bkg_tab);

  TGHorizontalFrame *btn_frame = new TGHorizontalFrame(controls, 510, 40);
  controls->AddFrame(
      btn_frame, new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 5, 5));

  TGTextButton *refit_btn = new TGTextButton(btn_frame, "Refit", kBtnRefit);
  refit_btn->Associate(this);
  TGTextButton *accept_btn = new TGTextButton(btn_frame, "Accept", kBtnAccept);
  accept_btn->Associate(this);
  TGTextButton *cancel_btn = new TGTextButton(btn_frame, "Cancel", kBtnCancel);
  cancel_btn->Associate(this);
  TGTextButton *reset_btn = new TGTextButton(btn_frame, "Reset", kBtnReset);
  reset_btn->Associate(this);

  TGLayoutHints *btn_hints = new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2);
  btn_frame->AddFrame(refit_btn, btn_hints);
  btn_frame->AddFrame(accept_btn, btn_hints);
  btn_frame->AddFrame(cancel_btn, btn_hints);
  btn_frame->AddFrame(reset_btn, btn_hints);
}

void InteractiveRooFitEditor::BuildPeakTab(TGCompositeFrame *parent,
                                           Int_t peak_idx) {
  Int_t offset = peak_idx * 10;

  TGGroupFrame *gaus_grp = new TGGroupFrame(parent, "Gaussian", kVerticalFrame);
  parent->AddFrame(gaus_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 1));
  AddParamRow(gaus_grp, offset + 0, "Mu");
  AddParamRow(gaus_grp, offset + 1, "Sigma");
  AddParamRow(gaus_grp, offset + 2, "Yield");

  TGGroupFrame *step_grp = new TGGroupFrame(parent, "Step", kVerticalFrame);
  parent->AddFrame(step_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(step_grp, offset + 3, "Step Ratio");

  TGGroupFrame *ltail_grp =
      new TGGroupFrame(parent, "Low-Side Tails", kVerticalFrame);
  parent->AddFrame(ltail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(ltail_grp, offset + 4, "Exp Amp");
  AddParamRow(ltail_grp, offset + 5, "Exp Decay/Sigma");
  AddParamRow(ltail_grp, offset + 6, "Lin Ratio");
  AddParamRow(ltail_grp, offset + 7, "Lin Slope");

  TGGroupFrame *htail_grp =
      new TGGroupFrame(parent, "High-Side Tail", kVerticalFrame);
  parent->AddFrame(htail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 3));
  AddParamRow(htail_grp, offset + 8, "Exp Amp");
  AddParamRow(htail_grp, offset + 9, "Exp Decay/Sigma");
}

void InteractiveRooFitEditor::BuildBackgroundTab(TGCompositeFrame *parent) {
  TGGroupFrame *bkg_grp =
      new TGGroupFrame(parent, "Background", kVerticalFrame);
  parent->AddFrame(bkg_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
  AddParamRow(bkg_grp, BkgConstIdx(), "Yield");
  AddParamRow(bkg_grp, BkgSlopeIdx(), "Slope");
}

void InteractiveRooFitEditor::AddParamRow(TGCompositeFrame *parent,
                                          Int_t param_idx, const char *name) {
  TGHorizontalFrame *row = new TGHorizontalFrame(parent, 490, 28);
  parent->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

  Double_t val = params_[param_idx]->getVal();
  Bool_t fixed = original_fixed_[param_idx];

  TGLabel *label = new TGLabel(row, name);
  row->AddFrame(label, new TGLayoutHints(kLHintsCenterY, 2, 4, 2, 2));
  label->Resize(90, 20);

  TGHSlider *slider =
      new TGHSlider(row, 140, kSlider1, kSliderBase + param_idx);
  slider->SetRange(0, kSliderRes);
  slider->SetPosition(ValToSlider(param_idx, val));
  slider->Associate(this);
  row->AddFrame(slider,
                new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));
  sliders_[param_idx] = slider;

  TGNumberEntry *entry = new TGNumberEntry(
      row, val, 8, kEntryBase + param_idx, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  entry->GetNumberEntry()->Associate(this);
  row->AddFrame(entry, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  value_entries_[param_idx] = entry;

  TGCheckButton *fix_cb = new TGCheckButton(row, "Fix", kFixBase + param_idx);
  fix_cb->Associate(this);
  if (fixed) {
    fix_cb->SetState(kButtonDown);
  }
  row->AddFrame(fix_cb, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  fix_checks_[param_idx] = fix_cb;

  TGNumberEntry *lo_entry = new TGNumberEntry(
      row, current_bounds_low_[param_idx], 6, kLoBoundBase + param_idx,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  lo_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(lo_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  lo_bound_entries_[param_idx] = lo_entry;

  TGNumberEntry *hi_entry = new TGNumberEntry(
      row, current_bounds_high_[param_idx], 6, kHiBoundBase + param_idx,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  hi_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(hi_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  hi_bound_entries_[param_idx] = hi_entry;

  if (fixed) {
    slider->SetEnabled(kFALSE);
    entry->GetNumberEntry()->SetEnabled(kFALSE);
  }
}

Double_t InteractiveRooFitEditor::EvalPdfDensity(RooAbsPdf *pdf, Double_t xv) {
  Double_t saved = x_->getVal();
  x_->setVal(xv);
  RooArgSet nset(*x_);
  Double_t v = pdf->getVal(&nset);
  x_->setVal(saved);
  return v;
}

Double_t InteractiveRooFitEditor::ComponentExpected(RooAbsPdf *pdf,
                                                    Double_t yield,
                                                    Double_t xv) {
  Double_t bin_width = hist_->GetBinWidth(1);
  return yield * EvalPdfDensity(pdf, xv) * bin_width;
}

void InteractiveRooFitEditor::InitDrawing() {
  TCanvas *canvas = embedded_canvas_->GetCanvas();
  canvas->cd();

  main_pad_ = new TPad("rfeditor_main", "rfeditor_main", 0, 0.3, 1, 1.0);
  main_pad_->SetBottomMargin(0.04);
  main_pad_->SetGridx(1);
  main_pad_->SetGridy(1);
  main_pad_->SetTopMargin(0.08);
  main_pad_->Draw();

  residual_pad_ = new TPad("rfeditor_res", "rfeditor_res", 0, 0, 1, 0.3);
  residual_pad_->SetTopMargin(0.04);
  residual_pad_->SetBottomMargin(0.35);
  residual_pad_->SetGridx(1);
  residual_pad_->SetGridy(1);
  residual_pad_->Draw();

  main_pad_->cd();

  hist_draw_ = static_cast<TH1 *>(hist_->Clone("hist_rfeditor_draw"));
  hist_draw_->SetDirectory(0);
  Double_t min_disp = 0.9 * range_low_;
  Double_t max_disp = 1.1 * range_high_;
  hist_draw_->GetXaxis()->SetRangeUser(min_disp, max_disp);
  hist_draw_->GetXaxis()->SetLabelSize(0);
  hist_draw_->GetXaxis()->SetTitleSize(0);
  hist_draw_->SetLineColor(kViolet);
  hist_draw_->SetLineWidth(2);
  hist_draw_->Draw();

  total_graph_ = new TGraph(kNDrawPts);
  total_graph_->SetLineColor(kAzure);
  total_graph_->SetLineWidth(2);
  total_graph_->Draw("L same");

  bkg_graph_ = new TGraph(kNDrawPts);
  bkg_graph_->SetLineColor(kGreen);
  bkg_graph_->SetLineWidth(2);
  bkg_graph_->Draw("L same");

  Int_t comp_colors[4] = {kBlack, kGray, kRed, kOrange};
  for (Int_t p = 0; p < num_peaks_; p++) {
    for (Int_t c = 0; c < 4; c++) {
      comp_graphs_[p][c] = new TGraph(kNDrawPts);
      comp_graphs_[p][c]->SetLineColor(comp_colors[c]);
      comp_graphs_[p][c]->SetLineStyle(PeakStyle(p));
      comp_graphs_[p][c]->SetLineWidth(2);
      comp_graphs_[p][c]->Draw("L same");
    }
  }

  chi2_label_ = new TLatex();
  chi2_label_->SetNDC();
  chi2_label_->SetTextSize(0.045);
  chi2_label_->SetTextAlign(31);
  chi2_label_->Draw();

  main_pad_->SetLogy(kTRUE);

  residual_pad_->cd();

  n_res_points_ = 0;
  Int_t nbins = hist_->GetNbinsX();
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t xv = hist_->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t data = hist_->GetBinContent(i);
    Double_t error = hist_->GetBinError(i);
    if (error > 0 && data > 0)
      n_res_points_++;
  }
  if (n_res_points_ < 1)
    n_res_points_ = 1;

  res_graph_ = new TGraph(n_res_points_);
  res_graph_->SetMarkerStyle(20);
  res_graph_->SetMarkerSize(0.8);
  res_graph_->SetMarkerColor(kAzure);
  res_graph_->SetLineColor(kAzure);
  res_graph_->SetTitle("");

  Double_t ax_min =
      hist_draw_->GetXaxis()->GetBinLowEdge(hist_draw_->GetXaxis()->GetFirst());
  Double_t ax_max =
      hist_draw_->GetXaxis()->GetBinUpEdge(hist_draw_->GetXaxis()->GetLast());
  res_graph_->GetXaxis()->SetLimits(ax_min, ax_max);
  res_graph_->GetYaxis()->SetTitle("#delta/#sigma");
  res_graph_->GetXaxis()->SetTitle(hist_->GetXaxis()->GetTitle());
  res_graph_->GetXaxis()->SetTitleSize(0.13);
  res_graph_->GetYaxis()->SetTitleSize(0.13);
  res_graph_->GetXaxis()->SetLabelSize(0.12);
  res_graph_->GetYaxis()->SetLabelSize(0.12);
  res_graph_->GetXaxis()->SetTitleOffset(1.0);
  res_graph_->GetYaxis()->SetTitleOffset(0.3);
  res_graph_->GetYaxis()->SetNdivisions(505);
  res_graph_->GetXaxis()->SetNdivisions(510);
  res_graph_->GetYaxis()->CenterTitle(kTRUE);
  res_graph_->GetYaxis()->SetRangeUser(-5.5, 5.5);
  res_graph_->Draw("AP");

  zero_line_ = new TF1("zero_rfeditor", "0", ax_min, ax_max);
  zero_line_->SetLineColor(kBlack);
  zero_line_->SetLineStyle(2);
  zero_line_->SetLineWidth(2);
  zero_line_->Draw("same");

  plus3_line_ = new TF1("plus3_rfeditor", "3", ax_min, ax_max);
  plus3_line_->SetLineColor(kGray + 2);
  plus3_line_->SetLineStyle(3);
  plus3_line_->SetLineWidth(2);
  plus3_line_->Draw("same");

  minus3_line_ = new TF1("minus3_rfeditor", "-3", ax_min, ax_max);
  minus3_line_->SetLineColor(kGray + 2);
  minus3_line_->SetLineStyle(3);
  minus3_line_->SetLineWidth(2);
  minus3_line_->Draw("same");

  UpdateAllGraphs();
  UpdateResPoints();
}

void InteractiveRooFitEditor::UpdateAllGraphs() {
  Double_t x_step = (range_high_ - range_low_) / (kNDrawPts - 1);
  Double_t bin_width = hist_->GetBinWidth(1);
  Double_t saved = x_->getVal();
  RooArgSet nset(*x_);

  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < kNDrawPts; i++) {
    Double_t xv = range_low_ + i * x_step;
    x_->setVal(xv);
    total_graph_->SetPoint(i, xv,
                           total_exp * total_pdf_->getVal(&nset) * bin_width);
  }

  Double_t bkg_yield = bkg_->bkg_yield->getVal();
  for (Int_t i = 0; i < kNDrawPts; i++) {
    Double_t xv = range_low_ + i * x_step;
    x_->setVal(xv);
    bkg_graph_->SetPoint(i, xv,
                         bkg_yield * bkg_->bkg_pdf->getVal(&nset) * bin_width);
  }

  for (Int_t pi = 0; pi < num_peaks_; pi++) {
    RooFitPeakModel &p = (*peaks_)[pi];
    Double_t gy = p.gaus_yield->getVal();
    Double_t sy = p.step_yield->getVal();
    Double_t lexp_y = p.low_exp_yield->getVal();
    Double_t llin_y = p.low_lin_yield->getVal();
    Double_t hexp_y = p.high_exp_yield->getVal();

    for (Int_t i = 0; i < kNDrawPts; i++) {
      Double_t xv = range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t bkg_v = bkg_yield * bkg_->bkg_pdf->getVal(&nset) * bin_width;

      comp_graphs_[pi][0]->SetPoint(
          i, xv, gy * p.gauss_pdf->getVal(&nset) * bin_width + bkg_v);
      comp_graphs_[pi][1]->SetPoint(
          i, xv, sy * p.step_pdf->getVal(&nset) * bin_width + bkg_v);
      Double_t y_low = lexp_y * p.low_exp_pdf->getVal(&nset) * bin_width +
                       llin_y * p.low_lin_pdf->getVal(&nset) * bin_width;
      comp_graphs_[pi][2]->SetPoint(i, xv, y_low + bkg_v);
      comp_graphs_[pi][3]->SetPoint(
          i, xv, hexp_y * p.high_exp_pdf->getVal(&nset) * bin_width + bkg_v);
    }
  }

  x_->setVal(saved);
}

void InteractiveRooFitEditor::UpdateResPoints() {
  Int_t nbins = hist_->GetNbinsX();
  Double_t bin_width = hist_->GetBinWidth(1);
  RooArgSet nset(*x_);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  Double_t saved = x_->getVal();

  Int_t pt = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t xv = hist_->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t data = hist_->GetBinContent(i);
    Double_t error = hist_->GetBinError(i);
    if (error > 0 && data > 0) {
      x_->setVal(xv);
      Double_t fit_val = total_exp * total_pdf_->getVal(&nset) * bin_width;
      Double_t pull = (data - fit_val) / error;
      res_graph_->SetPoint(pt, xv, pull);
      pt++;
    }
  }
  res_graph_->Set(pt);
  n_res_points_ = pt;

  Double_t disp_lo = 0.9 * range_low_;
  Double_t disp_hi = 1.1 * range_high_;
  res_graph_->GetXaxis()->SetLimits(disp_lo, disp_hi);

  x_->setVal(saved);
}

void InteractiveRooFitEditor::UpdateCanvas() {
  UpdateAllGraphs();
  UpdateResPoints();

  Double_t chi2_sum = 0;
  Int_t ndf_bins = 0;
  Int_t nbins = hist_->GetNbinsX();
  Double_t bin_width = hist_->GetBinWidth(1);
  RooArgSet nset(*x_);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  Double_t saved = x_->getVal();
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t xv = hist_->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t data = hist_->GetBinContent(i);
    Double_t error = hist_->GetBinError(i);
    if (error > 0 && data > 0) {
      x_->setVal(xv);
      Double_t fit_val = total_exp * total_pdf_->getVal(&nset) * bin_width;
      Double_t residual = (data - fit_val) / error;
      chi2_sum += residual * residual;
      ndf_bins++;
    }
  }
  x_->setVal(saved);
  Int_t n_free = 0;
  for (Int_t i = 0; i < num_params_; i++) {
    if (!IsFixed(i))
      n_free++;
  }
  Int_t ndf = ndf_bins - n_free;
  if (ndf > 0) {
    chi2_label_->SetText(0.85, 0.85,
                         Form("#chi^{2}/ndf = %.3f", chi2_sum / ndf));
  } else {
    chi2_label_->SetText(0.85, 0.85, "#chi^{2}/ndf = N/A");
  }

  zero_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);
  plus3_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);
  minus3_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);
  res_graph_->GetYaxis()->SetRangeUser(-5.5, 5.5);

  TCanvas *canvas = embedded_canvas_->GetCanvas();
  main_pad_->Modified();
  residual_pad_->Modified();
  canvas->Modified();
  canvas->Update();

  needs_redraw_ = kFALSE;
}

void InteractiveRooFitEditor::SyncAllWidgets() {
  syncing_ = kTRUE;
  for (Int_t i = 0; i < num_params_; i++) {
    SyncWidget(i);
  }
  syncing_ = kFALSE;
}

void InteractiveRooFitEditor::SyncWidget(Int_t param_idx) {
  Double_t val = params_[param_idx]->getVal();
  Bool_t fixed = params_[param_idx]->isConstant();

  if (!fixed && params_[param_idx]->hasMin() && params_[param_idx]->hasMax()) {
    current_bounds_low_[param_idx] = params_[param_idx]->getMin();
    current_bounds_high_[param_idx] = params_[param_idx]->getMax();
  }

  Double_t clamped = val;
  if (clamped < current_bounds_low_[param_idx])
    clamped = current_bounds_low_[param_idx];
  if (clamped > current_bounds_high_[param_idx])
    clamped = current_bounds_high_[param_idx];

  sliders_[param_idx]->SetPosition(ValToSlider(param_idx, clamped));
  value_entries_[param_idx]->SetNumber(val);

  if (fixed) {
    fix_checks_[param_idx]->SetState(kButtonDown);
    sliders_[param_idx]->SetEnabled(kFALSE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kFALSE);
  } else {
    fix_checks_[param_idx]->SetState(kButtonUp);
    sliders_[param_idx]->SetEnabled(kTRUE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kTRUE);
  }

  lo_bound_entries_[param_idx]->SetNumber(current_bounds_low_[param_idx]);
  hi_bound_entries_[param_idx]->SetNumber(current_bounds_high_[param_idx]);
}

void InteractiveRooFitEditor::OnSliderMoved(Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(param_idx))
    return;

  Int_t pos = sliders_[param_idx]->GetPosition();
  Double_t val = SliderToVal(param_idx, pos);
  params_[param_idx]->setVal(val);

  syncing_ = kTRUE;
  value_entries_[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveRooFitEditor::OnEntryChanged(Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(param_idx))
    return;

  Double_t val = value_entries_[param_idx]->GetNumber();

  Bool_t bounds_changed = kFALSE;
  if (val < current_bounds_low_[param_idx]) {
    current_bounds_low_[param_idx] = val;
    bounds_changed = kTRUE;
  }
  if (val > current_bounds_high_[param_idx]) {
    current_bounds_high_[param_idx] = val;
    bounds_changed = kTRUE;
  }

  params_[param_idx]->setRange(current_bounds_low_[param_idx],
                               current_bounds_high_[param_idx]);
  params_[param_idx]->setVal(val);

  syncing_ = kTRUE;
  sliders_[param_idx]->SetPosition(ValToSlider(param_idx, val));
  if (bounds_changed) {
    lo_bound_entries_[param_idx]->SetNumber(current_bounds_low_[param_idx]);
    hi_bound_entries_[param_idx]->SetNumber(current_bounds_high_[param_idx]);
  }
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveRooFitEditor::OnBoundsChanged(Int_t param_idx) {
  if (syncing_)
    return;

  Double_t new_lo = lo_bound_entries_[param_idx]->GetNumber();
  Double_t new_hi = hi_bound_entries_[param_idx]->GetNumber();

  if (new_lo >= new_hi)
    return;

  current_bounds_low_[param_idx] = new_lo;
  current_bounds_high_[param_idx] = new_hi;

  Double_t val = params_[param_idx]->getVal();
  if (val < new_lo)
    val = new_lo;
  if (val > new_hi)
    val = new_hi;

  if (!IsFixed(param_idx)) {
    params_[param_idx]->setRange(new_lo, new_hi);
  }
  params_[param_idx]->setVal(val);

  syncing_ = kTRUE;
  sliders_[param_idx]->SetPosition(ValToSlider(param_idx, val));
  value_entries_[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveRooFitEditor::OnFixToggled(Int_t param_idx) {
  Bool_t now_fixed = (fix_checks_[param_idx]->GetState() == kButtonDown);

  if (now_fixed) {
    params_[param_idx]->setConstant(kTRUE);
    sliders_[param_idx]->SetEnabled(kFALSE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kFALSE);
  } else {
    params_[param_idx]->setRange(current_bounds_low_[param_idx],
                                 current_bounds_high_[param_idx]);
    params_[param_idx]->setConstant(kFALSE);
    sliders_[param_idx]->SetEnabled(kTRUE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kTRUE);
  }

  needs_redraw_ = kTRUE;
}

void InteractiveRooFitEditor::OnRangeChanged() {
  if (syncing_)
    return;

  Double_t new_lo = range_slider_->GetMinPosition();
  Double_t new_hi = range_slider_->GetMaxPosition();

  if (new_lo >= new_hi)
    return;

  range_low_ = new_lo;
  range_high_ = new_hi;

  x_->setRange(range_low_, range_high_);
  x_->setRange("fitrange", range_low_, range_high_);

  if (events_) {
    RooFitUtils::RefillDisplayHistogram(
        hist_draw_, *events_, (Float_t)range_low_, (Float_t)range_high_,
        display_bin_width_kev_);
  }
  hist_draw_->GetXaxis()->SetRangeUser(0.9 * range_low_, 1.1 * range_high_);

  syncing_ = kTRUE;
  range_lo_entry_->SetNumber(range_low_);
  range_hi_entry_->SetNumber(range_high_);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveRooFitEditor::DoRefit() {
  for (Int_t i = 0; i < num_params_; i++) {
    if (!IsFixed(i)) {
      params_[i]->setRange(current_bounds_low_[i], current_bounds_high_[i]);
    }
  }

  RooFitResult *res = total_pdf_->fitTo(
      *data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE),
      RooFit::Range("fitrange"), RooFit::SumW2Error(kFALSE),
      RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), RooFit::Strategy(2),
      RooFit::Minimizer("Minuit2", "migrad"), BestAvailableBackend());
  if (res)
    delete res;

  SyncAllWidgets();
  needs_redraw_ = kTRUE;
  UpdateCanvas();
}

void InteractiveRooFitEditor::DoAccept() {
  accepted_ = kTRUE;
  done_ = kTRUE;
}

void InteractiveRooFitEditor::DoCancel() {
  for (Int_t i = 0; i < num_params_; i++) {
    params_[i]->setVal(original_params_[i]);
    if (original_fixed_[i]) {
      params_[i]->setConstant(kTRUE);
    } else {
      params_[i]->setRange(original_bounds_low_[i], original_bounds_high_[i]);
      params_[i]->setConstant(kFALSE);
    }
  }
  x_->setRange(original_range_low_, original_range_high_);
  x_->setRange("fitrange", original_range_low_, original_range_high_);
  accepted_ = kFALSE;
  done_ = kTRUE;
}

void InteractiveRooFitEditor::DoReset() {
  for (Int_t i = 0; i < num_params_; i++) {
    params_[i]->setVal(original_params_[i]);
    current_bounds_low_[i] = original_bounds_low_[i];
    current_bounds_high_[i] = original_bounds_high_[i];
    if (original_fixed_[i]) {
      params_[i]->setConstant(kTRUE);
    } else {
      params_[i]->setRange(original_bounds_low_[i], original_bounds_high_[i]);
      params_[i]->setConstant(kFALSE);
    }
  }

  range_low_ = original_range_low_;
  range_high_ = original_range_high_;
  x_->setRange(range_low_, range_high_);
  x_->setRange("fitrange", range_low_, range_high_);
  hist_draw_->GetXaxis()->SetRangeUser(0.9 * range_low_, 1.1 * range_high_);
  syncing_ = kTRUE;
  range_slider_->SetPosition(range_low_, range_high_);
  range_lo_entry_->SetNumber(range_low_);
  range_hi_entry_->SetNumber(range_high_);
  syncing_ = kFALSE;

  SyncAllWidgets();
  needs_redraw_ = kTRUE;
  UpdateCanvas();
}

Bool_t InteractiveRooFitEditor::ProcessMessage(Long_t msg, Long_t parm1,
                                               Long_t parm2) {
  if (syncing_)
    return kTRUE;

  (void)parm2;

  switch (GET_MSG(msg)) {
  case kC_COMMAND:
    switch (GET_SUBMSG(msg)) {
    case kCM_BUTTON:
      if (parm1 == kBtnRefit)
        DoRefit();
      else if (parm1 == kBtnAccept)
        DoAccept();
      else if (parm1 == kBtnCancel)
        DoCancel();
      else if (parm1 == kBtnReset)
        DoReset();
      break;
    case kCM_CHECKBUTTON:
      if (parm1 >= kFixBase && parm1 < kFixBase + num_params_) {
        OnFixToggled(parm1 - kFixBase);
      }
      break;
    }
    break;
  case kC_HSLIDER:
    if (parm1 == kRangeSlider) {
      OnRangeChanged();
    } else if (parm1 >= kSliderBase && parm1 < kSliderBase + num_params_) {
      OnSliderMoved(parm1 - kSliderBase);
    }
    break;
  case kC_TEXTENTRY:
    if (GET_SUBMSG(msg) == kTE_ENTER || GET_SUBMSG(msg) == kTE_TAB) {
      if (parm1 == kRangeLoEntry || parm1 == kRangeHiEntry) {
        Double_t new_lo = range_lo_entry_->GetNumber();
        Double_t new_hi = range_hi_entry_->GetNumber();
        if (new_lo < new_hi) {
          syncing_ = kTRUE;
          range_slider_->SetPosition(new_lo, new_hi);
          syncing_ = kFALSE;
          OnRangeChanged();
        }
      } else if (parm1 >= kEntryBase && parm1 < kEntryBase + num_params_) {
        OnEntryChanged(parm1 - kEntryBase);
      } else if (parm1 >= kLoBoundBase && parm1 < kLoBoundBase + num_params_) {
        OnBoundsChanged(parm1 - kLoBoundBase);
      } else if (parm1 >= kHiBoundBase && parm1 < kHiBoundBase + num_params_) {
        OnBoundsChanged(parm1 - kHiBoundBase);
      }
    }
    break;
  }

  return kTRUE;
}

Bool_t InteractiveRooFitEditor::HandleTimer(TTimer *timer) {
  if (timer == redraw_timer_ && needs_redraw_) {
    UpdateCanvas();
  }
  return kTRUE;
}

void InteractiveRooFitEditor::CloseWindow() { DoCancel(); }

Int_t InteractiveRooFitEditor::ValToSlider(Int_t param_idx, Double_t val) {
  Double_t lo = current_bounds_low_[param_idx];
  Double_t hi = current_bounds_high_[param_idx];
  if (hi <= lo)
    return 0;
  Double_t frac = (val - lo) / (hi - lo);
  if (frac < 0.0)
    frac = 0.0;
  if (frac > 1.0)
    frac = 1.0;
  return static_cast<Int_t>(frac * kSliderRes);
}

Double_t InteractiveRooFitEditor::SliderToVal(Int_t param_idx, Int_t pos) {
  Double_t lo = current_bounds_low_[param_idx];
  Double_t hi = current_bounds_high_[param_idx];
  Double_t frac = static_cast<Double_t>(pos) / kSliderRes;
  return lo + frac * (hi - lo);
}

Bool_t InteractiveRooFitEditor::IsFixed(Int_t param_idx) {
  return (fix_checks_[param_idx]->GetState() == kButtonDown);
}

void InteractiveRooFitEditor::GetDefaultBounds(Int_t param_idx, Double_t &lo,
                                               Double_t &hi) {
  Double_t peak_height = hist_->GetMaximum();
  Double_t range_width = range_high_ - range_low_;

  if (param_idx == BkgConstIdx()) {
    lo = 0;
    hi = peak_height * range_width * 10.0;
    return;
  }
  if (param_idx == BkgSlopeIdx()) {
    lo = -1000.0;
    hi = 1000.0;
    return;
  }

  Int_t local = param_idx % 10;
  switch (local) {
  case 0:
    lo = range_low_;
    hi = range_high_;
    break;
  case 1:
    lo = range_width * 0.001;
    hi = range_width * 0.5;
    break;
  case 2:
    lo = 0;
    hi = peak_height * range_width * 10.0;
    break;
  case 3:
  case 4:
  case 6:
  case 8:
    lo = 0;
    hi = 0.5;
    break;
  case 5:
  case 9:
    lo = 1.0;
    hi = 100.0;
    break;
  case 7:
    lo = -1.0;
    hi = 1.0;
    break;
  default:
    lo = -1000;
    hi = 1000;
    break;
  }
}

Int_t InteractiveRooFitEditor::PeakStyle(Int_t peak_idx) {
  if (peak_idx == 0)
    return 1;
  if (peak_idx == 1)
    return 3;
  return 4;
}

Bool_t LaunchInteractiveRooFitEditor(
    TH1 *hist, const std::vector<Double_t> *events,
    Float_t display_bin_width_kev, RooAbsPdf *total_pdf, RooRealVar *x,
    RooAbsData *data, std::vector<RooFitPeakModel> *peaks,
    RooFitBackgroundModel *bkg, Double_t range_low, Double_t range_high,
    const TString &info_label) {
  if (!gClient) {
    std::cerr << "InteractiveRooFitEditor: GUI not available (gClient is null)."
              << std::endl;
    return kFALSE;
  }

  // Swallow recoverable X protocol errors for the editor's lifetime so a
  // transient bad redraw doesn't trip ROOT's crashing default handler.
  AUXErrorHandlerSave xerr_save = AUInstallTolerantXErrorHandler();

  InteractiveRooFitEditor *editor = new InteractiveRooFitEditor(
      gClient->GetRoot(), hist, events, display_bin_width_kev, total_pdf, x,
      data, peaks, bkg, range_low, range_high, info_label);

  while (!editor->IsDone()) {
    gSystem->ProcessEvents();
    gSystem->Sleep(10);
  }

  Bool_t result = editor->WasAccepted();

  editor->GetRedrawTimer()->TurnOff();

  // Prevent TCanvas::Close() from doing X11 operations during teardown.
  // Removing the canvas from the canvases list also suppresses X11 events
  // for it on the next ProcessEvents (otherwise the next editor's first
  // ProcessEvents picks up stale paint events and crashes in DrawString).
  // Batch mode disables remaining X11 drawing calls.
  // Note: removing the canvas from this list disables RecursiveRemove's
  // pad-primitive cleanup, so the editor destructor explicitly clears
  // pad primitive lists before deleting drawn objects.
  TCanvas *ecanvas = editor->GetEmbeddedCanvas()->GetCanvas();
  if (ecanvas) {
    gROOT->GetListOfCanvases()->Remove(ecanvas);
    ecanvas->SetBatch(kTRUE);
  }

  editor->DontCallClose();
  editor->UnmapWindow();
  gSystem->ProcessEvents();
  delete editor;

  AURestoreXErrorHandler(xerr_save);
  return result;
}
