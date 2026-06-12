#include "InteractiveSimultaneousFitEditor.hpp"

#include "InteractiveEditorX11Guard.hpp"
#include "PlottingUtils.hpp"

#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <TROOT.h>
#include <cmath>
#include <iostream>

InteractiveSimultaneousFitEditor::InteractiveSimultaneousFitEditor(
    const TGWindow *parent, RooSimultaneous *sim_pdf, RooAbsData *combined_data,
    RooRealVar *x, const std::vector<SimEditorChannelView> &channel_views,
    Double_t range_low, Double_t range_high, const TString &info_label,
    Bool_t fit_debug)
    : TGMainFrame(parent, 1500, 950) {

  sim_pdf_ = sim_pdf;
  combined_data_ = combined_data;
  x_ = x;
  channels_ = channel_views;
  info_label_text_ = info_label;
  fit_debug_ = fit_debug;
  range_low_ = range_low;
  range_high_ = range_high;
  original_range_low_ = range_low;
  original_range_high_ = range_high;

  hist_x_min_ = channels_[0].hist->GetXaxis()->GetXmin();
  hist_x_max_ = channels_[0].hist->GetXaxis()->GetXmax();
  for (size_t i = 1; i < channels_.size(); i++) {
    Double_t cmin = channels_[i].hist->GetXaxis()->GetXmin();
    Double_t cmax = channels_[i].hist->GetXaxis()->GetXmax();
    if (cmin < hist_x_min_)
      hist_x_min_ = cmin;
    if (cmax > hist_x_max_)
      hist_x_max_ = cmax;
  }

  accepted_ = kFALSE;
  done_ = kFALSE;
  needs_redraw_ = kFALSE;
  syncing_ = kFALSE;

  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    cv.params.clear();
    for (Int_t pi = 0; pi < cv.num_peaks; pi++) {
      RooFitPeakModel &p = (*cv.peaks)[pi];
      cv.params.push_back(p.mu);
      cv.params.push_back(p.sigma);
      cv.params.push_back(p.gaus_yield);
      cv.params.push_back(p.ratio_step);
      cv.params.push_back(p.ratio_low_exp);
      cv.params.push_back(p.tau_low_exp);
      cv.params.push_back(p.ratio_low_lin);
      cv.params.push_back(p.slope_low_lin);
      cv.params.push_back(p.ratio_high_exp);
      cv.params.push_back(p.tau_high_exp);
    }
    cv.params.push_back(cv.bkg->bkg_yield);
    cv.params.push_back(cv.bkg->bkg_slope);

    Int_t np = (Int_t)cv.params.size();
    cv.original_params.resize(np);
    cv.original_bounds_low.resize(np);
    cv.original_bounds_high.resize(np);
    cv.original_fixed.resize(np);
    cv.current_bounds_low.resize(np);
    cv.current_bounds_high.resize(np);
    cv.sliders.resize(np);
    cv.value_entries.resize(np);
    cv.fix_checks.resize(np);
    cv.lo_bound_entries.resize(np);
    cv.hi_bound_entries.resize(np);

    for (Int_t i = 0; i < np; i++) {
      cv.original_params[i] = cv.params[i]->getVal();
      cv.original_fixed[i] = cv.params[i]->isConstant();
      Double_t lo, hi;
      if (cv.params[i]->hasMin() && cv.params[i]->hasMax()) {
        lo = cv.params[i]->getMin();
        hi = cv.params[i]->getMax();
      } else {
        GetDefaultBounds(ci, i, lo, hi);
      }
      cv.original_bounds_low[i] = lo;
      cv.original_bounds_high[i] = hi;
      cv.current_bounds_low[i] = lo;
      cv.current_bounds_high[i] = hi;
    }

    cv.embedded_canvas = nullptr;
    cv.main_pad = nullptr;
    cv.residual_pad = nullptr;
    cv.hist_draw = nullptr;
    cv.total_graph = nullptr;
    cv.bkg_graph = nullptr;
    cv.res_graph = nullptr;
    cv.zero_line = nullptr;
    cv.plus3_line = nullptr;
    cv.minus3_line = nullptr;
    cv.chi2_label = nullptr;
    cv.n_res_points = 0;
    for (Int_t p = 0; p < 3; p++) {
      for (Int_t c = 0; c < 4; c++) {
        cv.comp_graphs[p][c] = nullptr;
      }
    }
  }

  ApplyBackgroundSlopeBounds();

  BuildGUI();
  InitDrawing();

  redraw_timer_ = new TTimer(this, 50);
  redraw_timer_->TurnOn();

  SetWindowName("Interactive Simultaneous Fit Editor");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  UpdateCanvases();
}

InteractiveSimultaneousFitEditor::~InteractiveSimultaneousFitEditor() {
  if (redraw_timer_) {
    redraw_timer_->TurnOff();
    delete redraw_timer_;
  }

  // Drawing primitives drawn on each channel's embedded canvas — TPad does
  // not auto-delete user primitives, so they have to go by hand.
  //
  // Pre-clear each pad's primitives list with "nodelete" before destroying
  // the primitives themselves. The launcher removes the embedded canvases
  // from gROOT->GetListOfCanvases() to suppress X11 errors during teardown,
  // which also disables RecursiveRemove's cleanup of pad lists; without
  // this manual Clear, ~TPad later walks stale pointers and segfaults.
  for (Int_t ci = 0; ci < (Int_t)channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    if (cv.main_pad && cv.main_pad->GetListOfPrimitives()) {
      cv.main_pad->GetListOfPrimitives()->Clear("nodelete");
    }
    if (cv.residual_pad && cv.residual_pad->GetListOfPrimitives()) {
      cv.residual_pad->GetListOfPrimitives()->Clear("nodelete");
    }
    for (Int_t p = 0; p < 3; p++) {
      for (Int_t c = 0; c < 4; c++) {
        delete cv.comp_graphs[p][c];
      }
    }
    delete cv.chi2_label;
    delete cv.zero_line;
    delete cv.plus3_line;
    delete cv.minus3_line;
    delete cv.res_graph;
    delete cv.bkg_graph;
    delete cv.total_graph;
    delete cv.hist_draw;
    delete cv.residual_pad;
    delete cv.main_pad;
  }
}

Int_t InteractiveSimultaneousFitEditor::ChannelIndexFromWidgetId(
    Int_t base, Int_t parm1, Int_t &local_idx) {
  Int_t offset = parm1 - base;
  Int_t ch_idx = offset / kSliderStride;
  local_idx = offset % kSliderStride;
  if (ch_idx < 0 || ch_idx >= (Int_t)channels_.size())
    return -1;
  if (local_idx < 0 || local_idx >= (Int_t)channels_[ch_idx].params.size())
    return -1;
  return ch_idx;
}

void InteractiveSimultaneousFitEditor::BuildGUI() {
  TGHorizontalFrame *main_frame = new TGHorizontalFrame(this, 1500, 900);
  AddFrame(main_frame,
           new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  TGTab *outer_tabs = new TGTab(main_frame, 1500, 900);
  main_frame->AddFrame(
      outer_tabs,
      new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  for (size_t ci = 0; ci < channels_.size(); ci++) {
    TGCompositeFrame *channel_tab =
        outer_tabs->AddTab(channels_[ci].name.Data());
    BuildChannelTab(channel_tab, (Int_t)ci);
  }

  TGHorizontalFrame *bottom = new TGHorizontalFrame(this, 1500, 80);
  AddFrame(bottom,
           new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 2, 5));

  if (info_label_text_.Length() > 0) {
    TGLabel *info_label = new TGLabel(bottom, info_label_text_);
    bottom->AddFrame(info_label, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                                   8, 8, 4, 4));
  }

  TGGroupFrame *range_grp =
      new TGGroupFrame(bottom, "Fit Range (shared)", kHorizontalFrame);
  bottom->AddFrame(range_grp,
                   new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 4, 0, 0));

  range_slider_ =
      new TGDoubleHSlider(range_grp, 220, kDoubleScaleNo, kRangeSlider);
  range_slider_->SetRange(hist_x_min_, hist_x_max_);
  range_slider_->SetPosition(range_low_, range_high_);
  range_slider_->Associate(this);
  range_grp->AddFrame(range_slider_,
                      new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

  range_lo_entry_ = new TGNumberEntry(
      range_grp, range_low_, 7, kRangeLoEntry, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  range_lo_entry_->GetNumberEntry()->Associate(this);
  range_grp->AddFrame(range_lo_entry_,
                      new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

  range_hi_entry_ = new TGNumberEntry(
      range_grp, range_high_, 7, kRangeHiEntry, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  range_hi_entry_->GetNumberEntry()->Associate(this);
  range_grp->AddFrame(range_hi_entry_,
                      new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

  TGTextButton *refit_btn = new TGTextButton(bottom, "Refit", kBtnRefit);
  refit_btn->Associate(this);
  TGTextButton *accept_btn = new TGTextButton(bottom, "Accept", kBtnAccept);
  accept_btn->Associate(this);
  TGTextButton *cancel_btn = new TGTextButton(bottom, "Cancel", kBtnCancel);
  cancel_btn->Associate(this);
  TGTextButton *reset_btn = new TGTextButton(bottom, "Reset", kBtnReset);
  reset_btn->Associate(this);

  TGLayoutHints *btn_hints =
      new TGLayoutHints(kLHintsRight | kLHintsCenterY, 4, 4, 4, 4);
  bottom->AddFrame(reset_btn, btn_hints);
  bottom->AddFrame(cancel_btn, btn_hints);
  bottom->AddFrame(accept_btn, btn_hints);
  bottom->AddFrame(refit_btn, btn_hints);
}

void InteractiveSimultaneousFitEditor::BuildChannelTab(TGCompositeFrame *parent,
                                                       Int_t ch_idx) {
  SimEditorChannelView &cv = channels_[ch_idx];

  TGHorizontalFrame *row = new TGHorizontalFrame(parent, 1490, 850);
  parent->AddFrame(
      row, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  TString canvas_name = "SimEditorCanvas_" + cv.name;
  cv.embedded_canvas =
      new TRootEmbeddedCanvas(canvas_name.Data(), row, 850, 800);
  row->AddFrame(cv.embedded_canvas,
                new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  TGVerticalFrame *controls = new TGVerticalFrame(row, 540, 800);
  row->AddFrame(controls,
                new TGLayoutHints(kLHintsExpandY | kLHintsRight, 2, 2, 2, 2));

  TGTab *peak_tabs = new TGTab(controls, 530, 800);
  controls->AddFrame(
      peak_tabs,
      new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  for (Int_t pi = 0; pi < cv.num_peaks; pi++) {
    TString tab_name = TString::Format("Peak %d", pi + 1);
    TGCompositeFrame *pf = peak_tabs->AddTab(tab_name);
    BuildPeakSubTab(pf, ch_idx, pi);
  }
  TGCompositeFrame *bf = peak_tabs->AddTab("Background");
  BuildBackgroundSubTab(bf, ch_idx);
}

void InteractiveSimultaneousFitEditor::BuildPeakSubTab(TGCompositeFrame *parent,
                                                       Int_t ch_idx,
                                                       Int_t peak_idx) {
  Int_t offset = peak_idx * 10;

  TGGroupFrame *gaus_grp = new TGGroupFrame(parent, "Gaussian", kVerticalFrame);
  parent->AddFrame(gaus_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 1));
  AddParamRow(gaus_grp, ch_idx, offset + 0, "Mu");
  AddParamRow(gaus_grp, ch_idx, offset + 1, "Sigma");
  AddParamRow(gaus_grp, ch_idx, offset + 2, "Yield");

  TGGroupFrame *step_grp = new TGGroupFrame(parent, "Step", kVerticalFrame);
  parent->AddFrame(step_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(step_grp, ch_idx, offset + 3, "Step Ratio");

  TGGroupFrame *ltail_grp =
      new TGGroupFrame(parent, "Low-Side Tails", kVerticalFrame);
  parent->AddFrame(ltail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(ltail_grp, ch_idx, offset + 4, "Exp Amp");
  AddParamRow(ltail_grp, ch_idx, offset + 5, "Exp Decay/Sigma");
  AddParamRow(ltail_grp, ch_idx, offset + 6, "Lin Ratio");
  AddParamRow(ltail_grp, ch_idx, offset + 7, "Lin Slope");

  TGGroupFrame *htail_grp =
      new TGGroupFrame(parent, "High-Side Tail", kVerticalFrame);
  parent->AddFrame(htail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 3));
  AddParamRow(htail_grp, ch_idx, offset + 8, "Exp Amp");
  AddParamRow(htail_grp, ch_idx, offset + 9, "Exp Decay/Sigma");
}

void InteractiveSimultaneousFitEditor::BuildBackgroundSubTab(
    TGCompositeFrame *parent, Int_t ch_idx) {
  TGGroupFrame *bkg_grp =
      new TGGroupFrame(parent, "Background", kVerticalFrame);
  parent->AddFrame(bkg_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
  AddParamRow(bkg_grp, ch_idx, BkgConstIdx(ch_idx), "Yield");
  AddParamRow(bkg_grp, ch_idx, BkgSlopeIdx(ch_idx), "Slope");
}

void InteractiveSimultaneousFitEditor::AddParamRow(TGCompositeFrame *parent,
                                                   Int_t ch_idx,
                                                   Int_t param_idx,
                                                   const char *name) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Int_t global_id = ch_idx * kSliderStride + param_idx;

  TGHorizontalFrame *row = new TGHorizontalFrame(parent, 510, 28);
  parent->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

  Double_t val = cv.params[param_idx]->getVal();
  Bool_t fixed = cv.original_fixed[param_idx];

  TGLabel *label = new TGLabel(row, name);
  row->AddFrame(label, new TGLayoutHints(kLHintsCenterY, 2, 4, 2, 2));
  label->Resize(90, 20);

  TGHSlider *slider =
      new TGHSlider(row, 140, kSlider1, kSliderBase + global_id);
  slider->SetRange(0, kSliderRes);
  slider->SetPosition(ValToSlider(ch_idx, param_idx, val));
  slider->Associate(this);
  row->AddFrame(slider,
                new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));
  cv.sliders[param_idx] = slider;

  TGNumberEntry *entry = new TGNumberEntry(
      row, val, 8, kEntryBase + global_id, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  entry->GetNumberEntry()->Associate(this);
  row->AddFrame(entry, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  cv.value_entries[param_idx] = entry;

  TGCheckButton *fix_cb = new TGCheckButton(row, "Fix", kFixBase + global_id);
  fix_cb->Associate(this);
  if (fixed)
    fix_cb->SetState(kButtonDown);
  row->AddFrame(fix_cb, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  cv.fix_checks[param_idx] = fix_cb;

  TGNumberEntry *lo_entry = new TGNumberEntry(
      row, cv.current_bounds_low[param_idx], 6, kLoBoundBase + global_id,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  lo_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(lo_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  cv.lo_bound_entries[param_idx] = lo_entry;

  TGNumberEntry *hi_entry = new TGNumberEntry(
      row, cv.current_bounds_high[param_idx], 6, kHiBoundBase + global_id,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  hi_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(hi_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  cv.hi_bound_entries[param_idx] = hi_entry;

  if (fixed) {
    slider->SetEnabled(kFALSE);
    entry->GetNumberEntry()->SetEnabled(kFALSE);
  }
}

void InteractiveSimultaneousFitEditor::InitDrawing() {
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    InitChannelDrawing(channels_[ci]);
  }
}

void InteractiveSimultaneousFitEditor::InitChannelDrawing(
    SimEditorChannelView &cv) {
  TCanvas *canvas = cv.embedded_canvas->GetCanvas();
  canvas->cd();

  TString main_name = "simed_main_" + cv.name;
  cv.main_pad = new TPad(main_name.Data(), main_name.Data(), 0, 0.3, 1, 1.0);
  cv.main_pad->SetBottomMargin(0.04);
  cv.main_pad->SetGridx(1);
  cv.main_pad->SetGridy(1);
  cv.main_pad->SetTopMargin(0.08);
  cv.main_pad->Draw();

  TString res_name = "simed_res_" + cv.name;
  cv.residual_pad = new TPad(res_name.Data(), res_name.Data(), 0, 0, 1, 0.3);
  cv.residual_pad->SetTopMargin(0.04);
  cv.residual_pad->SetBottomMargin(0.35);
  cv.residual_pad->SetGridx(1);
  cv.residual_pad->SetGridy(1);
  cv.residual_pad->Draw();

  cv.main_pad->cd();
  TString draw_name = "hist_simed_" + cv.name;
  cv.hist_draw = static_cast<TH1 *>(cv.hist->Clone(draw_name.Data()));
  cv.hist_draw->SetDirectory(0);
  cv.hist_draw->GetXaxis()->SetRangeUser(0.9 * range_low_, 1.1 * range_high_);
  cv.hist_draw->GetXaxis()->SetLabelSize(0);
  cv.hist_draw->GetXaxis()->SetTitleSize(0);
  cv.hist_draw->SetLineColor(kViolet);
  cv.hist_draw->SetLineWidth(2);
  cv.hist_draw->Draw();

  cv.total_graph = new TGraph(kNDrawPts);
  cv.total_graph->SetLineColor(kAzure);
  cv.total_graph->SetLineWidth(2);
  cv.total_graph->Draw("L same");

  cv.bkg_graph = new TGraph(kNDrawPts);
  cv.bkg_graph->SetLineColor(kGreen);
  cv.bkg_graph->SetLineWidth(2);
  cv.bkg_graph->Draw("L same");

  Int_t comp_colors[4] = {kBlack, kGray, kRed, kOrange};
  for (Int_t p = 0; p < cv.num_peaks; p++) {
    for (Int_t c = 0; c < 4; c++) {
      cv.comp_graphs[p][c] = new TGraph(kNDrawPts);
      cv.comp_graphs[p][c]->SetLineColor(comp_colors[c]);
      cv.comp_graphs[p][c]->SetLineStyle(PeakStyle(p));
      cv.comp_graphs[p][c]->SetLineWidth(2);
      cv.comp_graphs[p][c]->Draw("L same");
    }
  }

  cv.chi2_label = new TLatex();
  cv.chi2_label->SetNDC();
  cv.chi2_label->SetTextSize(0.045);
  cv.chi2_label->SetTextAlign(31);
  cv.chi2_label->Draw();

  cv.main_pad->SetLogy(kTRUE);

  cv.residual_pad->cd();
  Int_t nbh = cv.hist->GetNbinsX();
  cv.n_res_points = 0;
  for (Int_t i = 1; i <= nbh; i++) {
    Double_t xv = cv.hist->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t d = cv.hist->GetBinContent(i);
    Double_t e = cv.hist->GetBinError(i);
    if (e > 0 && d > 0)
      cv.n_res_points++;
  }
  if (cv.n_res_points < 1)
    cv.n_res_points = 1;

  cv.res_graph = new TGraph(cv.n_res_points);
  cv.res_graph->SetMarkerStyle(20);
  cv.res_graph->SetMarkerSize(0.8);
  cv.res_graph->SetMarkerColor(kAzure);
  cv.res_graph->SetLineColor(kAzure);
  cv.res_graph->SetTitle("");

  Double_t ax_min = 0.9 * range_low_;
  Double_t ax_max = 1.1 * range_high_;
  cv.res_graph->GetXaxis()->SetLimits(ax_min, ax_max);
  cv.res_graph->GetYaxis()->SetTitle("#delta/#sigma");
  cv.res_graph->GetXaxis()->SetTitle(cv.hist->GetXaxis()->GetTitle());
  cv.res_graph->GetXaxis()->SetTitleSize(0.13);
  cv.res_graph->GetYaxis()->SetTitleSize(0.13);
  cv.res_graph->GetXaxis()->SetLabelSize(0.12);
  cv.res_graph->GetYaxis()->SetLabelSize(0.12);
  cv.res_graph->GetXaxis()->SetTitleOffset(1.0);
  cv.res_graph->GetYaxis()->SetTitleOffset(0.3);
  cv.res_graph->GetYaxis()->SetNdivisions(505);
  cv.res_graph->GetXaxis()->SetNdivisions(510);
  cv.res_graph->GetYaxis()->CenterTitle(kTRUE);
  cv.res_graph->GetYaxis()->SetRangeUser(-5.5, 5.5);
  cv.res_graph->Draw("AP");

  TString zero_name = "zero_simed_" + cv.name;
  cv.zero_line = new TF1(zero_name.Data(), "0", ax_min, ax_max);
  cv.zero_line->SetLineColor(kBlack);
  cv.zero_line->SetLineStyle(2);
  cv.zero_line->SetLineWidth(2);
  cv.zero_line->Draw("same");

  TString plus3_name = "plus3_simed_" + cv.name;
  cv.plus3_line = new TF1(plus3_name.Data(), "3", ax_min, ax_max);
  cv.plus3_line->SetLineColor(kGray + 2);
  cv.plus3_line->SetLineStyle(3);
  cv.plus3_line->SetLineWidth(2);
  cv.plus3_line->Draw("same");

  TString minus3_name = "minus3_simed_" + cv.name;
  cv.minus3_line = new TF1(minus3_name.Data(), "-3", ax_min, ax_max);
  cv.minus3_line->SetLineColor(kGray + 2);
  cv.minus3_line->SetLineStyle(3);
  cv.minus3_line->SetLineWidth(2);
  cv.minus3_line->Draw("same");
}

void InteractiveSimultaneousFitEditor::UpdateCanvases() {
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    UpdateChannelGraphs(channels_[ci]);
    UpdateChannelResiduals(channels_[ci]);
    UpdateChannelChi2(channels_[ci]);

    TCanvas *canvas = channels_[ci].embedded_canvas->GetCanvas();
    channels_[ci].main_pad->Modified();
    channels_[ci].residual_pad->Modified();
    canvas->Modified();
    canvas->Update();
  }
  needs_redraw_ = kFALSE;
}

void InteractiveSimultaneousFitEditor::UpdateChannelGraphs(
    SimEditorChannelView &cv) {
  Double_t x_step = (range_high_ - range_low_) / (kNDrawPts - 1);
  Double_t bin_width = cv.hist->GetBinWidth(1);
  Double_t saved = x_->getVal();
  RooArgSet nset(*x_);

  Double_t total_exp = cv.pdf->expectedEvents(&nset);
  for (Int_t i = 0; i < kNDrawPts; i++) {
    Double_t xv = range_low_ + i * x_step;
    x_->setVal(xv);
    cv.total_graph->SetPoint(i, xv,
                             total_exp * cv.pdf->getVal(&nset) * bin_width);
  }

  Double_t bkg_yield = cv.bkg->bkg_yield->getVal();
  for (Int_t i = 0; i < kNDrawPts; i++) {
    Double_t xv = range_low_ + i * x_step;
    x_->setVal(xv);
    cv.bkg_graph->SetPoint(
        i, xv, bkg_yield * cv.bkg->bkg_pdf->getVal(&nset) * bin_width);
  }

  for (Int_t pi = 0; pi < cv.num_peaks; pi++) {
    RooFitPeakModel &p = (*cv.peaks)[pi];
    Double_t gy = p.gaus_yield->getVal();
    Double_t sy = p.step_yield->getVal();
    Double_t lexp_y = p.low_exp_yield->getVal();
    Double_t llin_y = p.low_lin_yield->getVal();
    Double_t hexp_y = p.high_exp_yield->getVal();

    for (Int_t i = 0; i < kNDrawPts; i++) {
      Double_t xv = range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t bkg_v = bkg_yield * cv.bkg->bkg_pdf->getVal(&nset) * bin_width;

      cv.comp_graphs[pi][0]->SetPoint(
          i, xv, gy * p.gauss_pdf->getVal(&nset) * bin_width + bkg_v);
      cv.comp_graphs[pi][1]->SetPoint(
          i, xv, sy * p.step_pdf->getVal(&nset) * bin_width + bkg_v);
      Double_t y_low = lexp_y * p.low_exp_pdf->getVal(&nset) * bin_width +
                       llin_y * p.low_lin_pdf->getVal(&nset) * bin_width;
      cv.comp_graphs[pi][2]->SetPoint(i, xv, y_low + bkg_v);
      cv.comp_graphs[pi][3]->SetPoint(
          i, xv, hexp_y * p.high_exp_pdf->getVal(&nset) * bin_width + bkg_v);
    }
  }

  x_->setVal(saved);
}

void InteractiveSimultaneousFitEditor::UpdateChannelResiduals(
    SimEditorChannelView &cv) {
  Int_t nbh = cv.hist->GetNbinsX();
  Double_t bin_width = cv.hist->GetBinWidth(1);
  RooArgSet nset(*x_);
  Double_t total_exp = cv.pdf->expectedEvents(&nset);
  Double_t saved = x_->getVal();

  Int_t pt = 0;
  for (Int_t i = 1; i <= nbh; i++) {
    Double_t xv = cv.hist->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t d = cv.hist->GetBinContent(i);
    Double_t e = cv.hist->GetBinError(i);
    if (e > 0 && d > 0) {
      x_->setVal(xv);
      Double_t fit_val = total_exp * cv.pdf->getVal(&nset) * bin_width;
      Double_t pull = (d - fit_val) / e;
      cv.res_graph->SetPoint(pt, xv, pull);
      pt++;
    }
  }
  cv.res_graph->Set(pt);
  cv.n_res_points = pt;

  Double_t ax_min = 0.9 * range_low_;
  Double_t ax_max = 1.1 * range_high_;
  cv.res_graph->GetXaxis()->SetLimits(ax_min, ax_max);
  cv.zero_line->SetRange(ax_min, ax_max);
  cv.plus3_line->SetRange(ax_min, ax_max);
  cv.minus3_line->SetRange(ax_min, ax_max);
  cv.res_graph->GetYaxis()->SetRangeUser(-5.5, 5.5);

  x_->setVal(saved);
}

void InteractiveSimultaneousFitEditor::UpdateChannelChi2(
    SimEditorChannelView &cv) {
  Double_t chi2_sum = 0;
  Int_t ndf_bins = 0;
  Int_t nbh = cv.hist->GetNbinsX();
  Double_t bin_width = cv.hist->GetBinWidth(1);
  RooArgSet nset(*x_);
  Double_t total_exp = cv.pdf->expectedEvents(&nset);
  Double_t saved = x_->getVal();
  for (Int_t i = 1; i <= nbh; i++) {
    Double_t xv = cv.hist->GetBinCenter(i);
    if (xv < range_low_ || xv > range_high_)
      continue;
    Double_t d = cv.hist->GetBinContent(i);
    Double_t e = cv.hist->GetBinError(i);
    if (e > 0 && d > 0) {
      x_->setVal(xv);
      Double_t fit_val = total_exp * cv.pdf->getVal(&nset) * bin_width;
      Double_t r = (d - fit_val) / e;
      chi2_sum += r * r;
      ndf_bins++;
    }
  }
  x_->setVal(saved);

  Int_t n_free = 0;
  for (Int_t i = 0; i < (Int_t)cv.params.size(); i++) {
    if (!IsFixed((Int_t)(&cv - channels_.data()), i))
      n_free++;
  }
  Int_t ndf = ndf_bins - n_free;
  if (ndf > 0)
    cv.chi2_label->SetText(0.85, 0.85,
                           Form("#chi^{2}/ndf = %.3f", chi2_sum / ndf));
  else
    cv.chi2_label->SetText(0.85, 0.85, "#chi^{2}/ndf = N/A");
}

void InteractiveSimultaneousFitEditor::SyncAllWidgets() {
  syncing_ = kTRUE;
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    for (Int_t i = 0; i < (Int_t)channels_[ci].params.size(); i++) {
      SyncChannelWidget((Int_t)ci, i);
    }
  }
  syncing_ = kFALSE;
}

void InteractiveSimultaneousFitEditor::SyncChannelWidget(Int_t ch_idx,
                                                         Int_t param_idx) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t val = cv.params[param_idx]->getVal();
  Bool_t fixed = cv.params[param_idx]->isConstant();

  if (!fixed && cv.params[param_idx]->hasMin() &&
      cv.params[param_idx]->hasMax()) {
    cv.current_bounds_low[param_idx] = cv.params[param_idx]->getMin();
    cv.current_bounds_high[param_idx] = cv.params[param_idx]->getMax();
  }

  Double_t clamped = val;
  if (clamped < cv.current_bounds_low[param_idx])
    clamped = cv.current_bounds_low[param_idx];
  if (clamped > cv.current_bounds_high[param_idx])
    clamped = cv.current_bounds_high[param_idx];

  cv.sliders[param_idx]->SetPosition(ValToSlider(ch_idx, param_idx, clamped));
  cv.value_entries[param_idx]->SetNumber(val);

  if (fixed) {
    cv.fix_checks[param_idx]->SetState(kButtonDown);
    cv.sliders[param_idx]->SetEnabled(kFALSE);
    cv.value_entries[param_idx]->GetNumberEntry()->SetEnabled(kFALSE);
  } else {
    cv.fix_checks[param_idx]->SetState(kButtonUp);
    cv.sliders[param_idx]->SetEnabled(kTRUE);
    cv.value_entries[param_idx]->GetNumberEntry()->SetEnabled(kTRUE);
  }

  cv.lo_bound_entries[param_idx]->SetNumber(cv.current_bounds_low[param_idx]);
  cv.hi_bound_entries[param_idx]->SetNumber(cv.current_bounds_high[param_idx]);
}

void InteractiveSimultaneousFitEditor::OnSliderMoved(Int_t ch_idx,
                                                     Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(ch_idx, param_idx))
    return;

  SimEditorChannelView &cv = channels_[ch_idx];
  Int_t pos = cv.sliders[param_idx]->GetPosition();
  Double_t val = SliderToVal(ch_idx, param_idx, pos);
  cv.params[param_idx]->setVal(val);

  syncing_ = kTRUE;
  cv.value_entries[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::OnEntryChanged(Int_t ch_idx,
                                                      Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(ch_idx, param_idx))
    return;

  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t val = cv.value_entries[param_idx]->GetNumber();

  Bool_t bounds_changed = kFALSE;
  if (val < cv.current_bounds_low[param_idx]) {
    cv.current_bounds_low[param_idx] = val;
    bounds_changed = kTRUE;
  }
  if (val > cv.current_bounds_high[param_idx]) {
    cv.current_bounds_high[param_idx] = val;
    bounds_changed = kTRUE;
  }

  cv.params[param_idx]->setRange(cv.current_bounds_low[param_idx],
                                 cv.current_bounds_high[param_idx]);
  cv.params[param_idx]->setVal(val);

  syncing_ = kTRUE;
  cv.sliders[param_idx]->SetPosition(ValToSlider(ch_idx, param_idx, val));
  if (bounds_changed) {
    cv.lo_bound_entries[param_idx]->SetNumber(cv.current_bounds_low[param_idx]);
    cv.hi_bound_entries[param_idx]->SetNumber(
        cv.current_bounds_high[param_idx]);
  }
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::OnBoundsChanged(Int_t ch_idx,
                                                       Int_t param_idx) {
  if (syncing_)
    return;

  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t new_lo = cv.lo_bound_entries[param_idx]->GetNumber();
  Double_t new_hi = cv.hi_bound_entries[param_idx]->GetNumber();
  if (new_lo >= new_hi)
    return;

  cv.current_bounds_low[param_idx] = new_lo;
  cv.current_bounds_high[param_idx] = new_hi;

  Double_t val = cv.params[param_idx]->getVal();
  if (val < new_lo)
    val = new_lo;
  if (val > new_hi)
    val = new_hi;

  if (!IsFixed(ch_idx, param_idx)) {
    cv.params[param_idx]->setRange(new_lo, new_hi);
  }
  cv.params[param_idx]->setVal(val);

  syncing_ = kTRUE;
  cv.sliders[param_idx]->SetPosition(ValToSlider(ch_idx, param_idx, val));
  cv.value_entries[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::OnFixToggled(Int_t ch_idx,
                                                    Int_t param_idx) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Bool_t now_fixed = (cv.fix_checks[param_idx]->GetState() == kButtonDown);

  if (now_fixed) {
    cv.params[param_idx]->setConstant(kTRUE);
    cv.sliders[param_idx]->SetEnabled(kFALSE);
    cv.value_entries[param_idx]->GetNumberEntry()->SetEnabled(kFALSE);
  } else {
    cv.params[param_idx]->setRange(cv.current_bounds_low[param_idx],
                                   cv.current_bounds_high[param_idx]);
    cv.params[param_idx]->setConstant(kFALSE);
    cv.sliders[param_idx]->SetEnabled(kTRUE);
    cv.value_entries[param_idx]->GetNumberEntry()->SetEnabled(kTRUE);
  }
  needs_redraw_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::OnRangeChanged() {
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
  ApplyBackgroundSlopeBounds();

  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    if (cv.events) {
      RooFitUtils::RefillDisplayHistogram(
          cv.hist_draw, *cv.events, (Float_t)range_low_, (Float_t)range_high_,
          cv.display_bin_width_kev);
    }
    cv.hist_draw->GetXaxis()->SetRangeUser(0.9 * range_low_, 1.1 * range_high_);
  }

  syncing_ = kTRUE;
  range_lo_entry_->SetNumber(range_low_);
  range_hi_entry_->SetNumber(range_high_);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::DoRefit() {
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    for (Int_t i = 0; i < (Int_t)cv.params.size(); i++) {
      if (!IsFixed((Int_t)ci, i)) {
        cv.params[i]->setRange(cv.current_bounds_low[i],
                               cv.current_bounds_high[i]);
      }
    }
  }

  if (fit_debug_) {
    // Un-suppress RooFit eval errors so the pdf producing an invalid NLL is
    // named, and report the pre-fit NLL at the current parameter values.
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
    RooAbsReal *nll = sim_pdf_->createNLL(
        *combined_data_, RooFit::Extended(kTRUE), RooFit::Range("fitrange"));
    Double_t nll_seed = (nll != 0) ? nll->getVal() : 0.0;
    std::cout << "=== AU_ROOFIT_FIT_DEBUG: pre-refit NLL = " << nll_seed
              << (std::isfinite(nll_seed) ? "" : "   <== NON-FINITE")
              << " ===" << std::endl;
    if (nll != 0)
      delete nll;
    DiagnoseInvalidComponents();
  }

  Int_t print_level = fit_debug_ ? 1 : -1;
  Int_t eval_errors = fit_debug_ ? 10 : -1;
  RooFitResult *res = sim_pdf_->fitTo(
      *combined_data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE),
      RooFit::Range("fitrange"), RooFit::SumW2Error(kFALSE),
      RooFit::PrintLevel(print_level), RooFit::PrintEvalErrors(eval_errors),
      RooFit::Strategy(1), RooFit::Minimizer("Minuit2", "migrad"),
      BestAvailableBackend());
  if (res)
    delete res;

  SyncAllWidgets();
  needs_redraw_ = kTRUE;
  UpdateCanvases();
}

void InteractiveSimultaneousFitEditor::DoAccept() {
  accepted_ = kTRUE;
  done_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::DoCancel() {
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    for (Int_t i = 0; i < (Int_t)cv.params.size(); i++) {
      cv.params[i]->setVal(cv.original_params[i]);
      if (cv.original_fixed[i]) {
        cv.params[i]->setConstant(kTRUE);
      } else {
        cv.params[i]->setRange(cv.original_bounds_low[i],
                               cv.original_bounds_high[i]);
        cv.params[i]->setConstant(kFALSE);
      }
    }
  }
  x_->setRange(original_range_low_, original_range_high_);
  x_->setRange("fitrange", original_range_low_, original_range_high_);
  accepted_ = kFALSE;
  done_ = kTRUE;
}

void InteractiveSimultaneousFitEditor::DoReset() {
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    for (Int_t i = 0; i < (Int_t)cv.params.size(); i++) {
      cv.params[i]->setVal(cv.original_params[i]);
      cv.current_bounds_low[i] = cv.original_bounds_low[i];
      cv.current_bounds_high[i] = cv.original_bounds_high[i];
      if (cv.original_fixed[i]) {
        cv.params[i]->setConstant(kTRUE);
      } else {
        cv.params[i]->setRange(cv.original_bounds_low[i],
                               cv.original_bounds_high[i]);
        cv.params[i]->setConstant(kFALSE);
      }
    }
  }
  range_low_ = original_range_low_;
  range_high_ = original_range_high_;
  x_->setRange(range_low_, range_high_);
  x_->setRange("fitrange", range_low_, range_high_);
  ApplyBackgroundSlopeBounds();
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    channels_[ci].hist_draw->GetXaxis()->SetRangeUser(0.9 * range_low_,
                                                      1.1 * range_high_);
  }
  syncing_ = kTRUE;
  range_slider_->SetPosition(range_low_, range_high_);
  range_lo_entry_->SetNumber(range_low_);
  range_hi_entry_->SetNumber(range_high_);
  syncing_ = kFALSE;

  SyncAllWidgets();
  needs_redraw_ = kTRUE;
  UpdateCanvases();
}

void InteractiveSimultaneousFitEditor::ApplyBackgroundSlopeBounds() {
  // The linear background 1+slope*x is only non-negative across [0, x_hi] when
  // slope >= -1/x_hi. The build-time bound used each channel's original fit
  // high edge; once the shared range is dragged wider, that bound lets the
  // polynomial dip negative inside the normalization window and poison the NLL.
  // Re-tie the lower bound to the current range high (with the same 0.9 margin
  // BuildChannelModel uses).
  if (range_high_ <= 0)
    return;
  Double_t slope_lo = -0.9 / range_high_;
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    Int_t si = BkgSlopeIdx((Int_t)ci);
    if (si < 0 || si >= (Int_t)cv.params.size())
      continue;
    RooRealVar *slope = cv.params[si];
    if (slope->isConstant())
      continue;
    Double_t hi = cv.current_bounds_high[si];
    if (slope_lo >= hi)
      continue;
    cv.current_bounds_low[si] = slope_lo;
    if (slope->getVal() < slope_lo)
      slope->setVal(slope_lo);
    slope->setRange(slope_lo, hi);
  }
}

void InteractiveSimultaneousFitEditor::DiagnoseInvalidComponents() {
  RooArgSet nset(*x_);
  for (size_t ci = 0; ci < channels_.size(); ci++) {
    SimEditorChannelView &cv = channels_[ci];
    std::cout << "  channel '" << cv.name << "' components:" << std::endl;
    // Per-component normalization integral (one value per parameter set): a NaN
    // or collapsed integral makes pdf/Int non-finite for the whole channel.
    for (size_t pi = 0; pi < cv.peaks->size(); pi++) {
      RooFitPeakModel &p = (*cv.peaks)[pi];
      RooAbsPdf *comps[5] = {p.gauss_pdf, p.step_pdf, p.low_exp_pdf,
                             p.low_lin_pdf, p.high_exp_pdf};
      const char *names[5] = {"gauss", "step", "lowExp", "lowLin", "highExp"};
      for (Int_t k = 0; k < 5; k++) {
        if (comps[k] == 0)
          continue;
        RooAbsReal *integ = comps[k]->createIntegral(
            nset, RooFit::NormSet(nset), RooFit::Range("fitrange"));
        Double_t nrm = (integ != 0) ? integ->getVal() : 0.0;
        if (integ != 0)
          delete integ;
        if (!std::isfinite(nrm) || nrm <= 0.0)
          std::cout << "    [BAD NORM] peak" << (pi + 1) << "." << names[k]
                    << " integral=" << nrm << std::endl;
      }
    }
    // Scan the range for any raw component value that is non-finite (e.g. a
    // tail overflowing to Inf*0=NaN), reporting the first few hits.
    const Int_t n_samp = 4000;
    Double_t step_x = (range_high_ - range_low_) / (Double_t)(n_samp - 1);
    Int_t reported = 0;
    for (Int_t i = 0; i < n_samp && reported < 8; i++) {
      Double_t xv = range_low_ + i * step_x;
      x_->setVal(xv);
      for (size_t pi = 0; pi < cv.peaks->size() && reported < 8; pi++) {
        RooFitPeakModel &p = (*cv.peaks)[pi];
        RooAbsPdf *comps[5] = {p.gauss_pdf, p.step_pdf, p.low_exp_pdf,
                               p.low_lin_pdf, p.high_exp_pdf};
        const char *names[5] = {"gauss", "step", "lowExp", "lowLin", "highExp"};
        for (Int_t k = 0; k < 5; k++) {
          if (comps[k] == 0)
            continue;
          Double_t r = comps[k]->getVal();
          if (!std::isfinite(r)) {
            std::cout << "    [BAD RAW] peak" << (pi + 1) << "." << names[k]
                      << " raw=" << r << " at x=" << xv << std::endl;
            reported++;
          }
        }
      }
    }
  }
}

Bool_t InteractiveSimultaneousFitEditor::ProcessMessage(Long_t msg,
                                                        Long_t parm1,
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
    case kCM_CHECKBUTTON: {
      Int_t local;
      Int_t ch = ChannelIndexFromWidgetId(kFixBase, parm1, local);
      if (ch >= 0)
        OnFixToggled(ch, local);
      break;
    }
    }
    break;
  case kC_HSLIDER:
    if (parm1 == kRangeSlider) {
      OnRangeChanged();
    } else {
      Int_t local;
      Int_t ch = ChannelIndexFromWidgetId(kSliderBase, parm1, local);
      if (ch >= 0)
        OnSliderMoved(ch, local);
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
      } else {
        Int_t local;
        Int_t ch = ChannelIndexFromWidgetId(kEntryBase, parm1, local);
        if (ch >= 0) {
          OnEntryChanged(ch, local);
        } else {
          ch = ChannelIndexFromWidgetId(kLoBoundBase, parm1, local);
          if (ch >= 0) {
            OnBoundsChanged(ch, local);
          } else {
            ch = ChannelIndexFromWidgetId(kHiBoundBase, parm1, local);
            if (ch >= 0)
              OnBoundsChanged(ch, local);
          }
        }
      }
    }
    break;
  }
  return kTRUE;
}

Bool_t InteractiveSimultaneousFitEditor::HandleTimer(TTimer *timer) {
  if (timer == redraw_timer_ && needs_redraw_) {
    UpdateCanvases();
  }
  return kTRUE;
}

void InteractiveSimultaneousFitEditor::CloseWindow() { DoCancel(); }

Int_t InteractiveSimultaneousFitEditor::ValToSlider(Int_t ch_idx,
                                                    Int_t param_idx,
                                                    Double_t val) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t lo = cv.current_bounds_low[param_idx];
  Double_t hi = cv.current_bounds_high[param_idx];
  if (hi <= lo)
    return 0;
  Double_t frac = (val - lo) / (hi - lo);
  if (frac < 0.0)
    frac = 0.0;
  if (frac > 1.0)
    frac = 1.0;
  return static_cast<Int_t>(frac * kSliderRes);
}

Double_t InteractiveSimultaneousFitEditor::SliderToVal(Int_t ch_idx,
                                                       Int_t param_idx,
                                                       Int_t pos) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t lo = cv.current_bounds_low[param_idx];
  Double_t hi = cv.current_bounds_high[param_idx];
  Double_t frac = static_cast<Double_t>(pos) / kSliderRes;
  return lo + frac * (hi - lo);
}

Bool_t InteractiveSimultaneousFitEditor::IsFixed(Int_t ch_idx,
                                                 Int_t param_idx) {
  return (channels_[ch_idx].fix_checks[param_idx]->GetState() == kButtonDown);
}

void InteractiveSimultaneousFitEditor::GetDefaultBounds(Int_t ch_idx,
                                                        Int_t param_idx,
                                                        Double_t &lo,
                                                        Double_t &hi) {
  SimEditorChannelView &cv = channels_[ch_idx];
  Double_t peak_height = cv.hist->GetMaximum();
  Double_t range_width = range_high_ - range_low_;

  if (param_idx == BkgConstIdx(ch_idx)) {
    lo = 0;
    hi = peak_height * range_width * 10.0;
    return;
  }
  if (param_idx == BkgSlopeIdx(ch_idx)) {
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
    lo = cv.hist->GetBinWidth(1);
    hi = range_width * 0.1;
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

Int_t InteractiveSimultaneousFitEditor::PeakStyle(Int_t peak_idx) {
  if (peak_idx == 0)
    return 1;
  if (peak_idx == 1)
    return 3;
  return 4;
}

Bool_t LaunchInteractiveSimultaneousFitEditor(
    RooSimultaneous *sim_pdf, RooAbsData *combined_data, RooRealVar *x,
    std::vector<SimEditorChannelView> &channel_views, Double_t range_low,
    Double_t range_high, const TString &info_label, Bool_t fit_debug) {
  if (!gClient) {
    std::cerr << "InteractiveSimultaneousFitEditor: GUI not available"
              << std::endl;
    return kFALSE;
  }

  // Swallow recoverable X protocol errors for the editor's lifetime so a
  // transient bad redraw doesn't trip ROOT's crashing default handler.
  AUXErrorHandlerSave xerr_save = AUInstallTolerantXErrorHandler();

  InteractiveSimultaneousFitEditor *editor =
      new InteractiveSimultaneousFitEditor(
          gClient->GetRoot(), sim_pdf, combined_data, x, channel_views,
          range_low, range_high, info_label, fit_debug);

  while (!editor->IsDone()) {
    gSystem->ProcessEvents();
    gSystem->Sleep(10);
  }

  Bool_t result = editor->WasAccepted();
  editor->GetRedrawTimer()->TurnOff();

  // Prevent TCanvas::Close() from doing X11 operations during teardown.
  // Removing the canvases from the canvases list also suppresses X11
  // events for them on the next ProcessEvents (otherwise the next editor's
  // first ProcessEvents picks up stale paint events and crashes in
  // DrawString). Batch mode disables remaining X11 drawing calls.
  // Note: removing the canvas from this list disables RecursiveRemove's
  // pad-primitive cleanup, so the editor destructor explicitly clears
  // pad primitive lists before deleting drawn objects.
  std::vector<SimEditorChannelView> *cvs = editor->GetChannels();
  for (size_t ci = 0; ci < cvs->size(); ci++) {
    TCanvas *c = (*cvs)[ci].embedded_canvas->GetCanvas();
    if (c) {
      gROOT->GetListOfCanvases()->Remove(c);
      c->SetBatch(kTRUE);
    }
  }

  editor->DontCallClose();
  editor->UnmapWindow();
  gSystem->ProcessEvents();
  delete editor;

  AURestoreXErrorHandler(xerr_save);
  return result;
}
