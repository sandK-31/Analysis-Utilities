#include "InteractiveFitEditor.hpp"

#include "InteractiveEditorX11Guard.hpp"

InteractiveFitEditor::InteractiveFitEditor(const TGWindow *parent, TH1 *hist,
                                           TF1 *fit_func, Double_t range_low,
                                           Double_t range_high, Int_t num_peaks,
                                           const TString &info_label)
    : TGMainFrame(parent, 1400, 900) {

  hist_ = hist;
  fit_func_ = fit_func;
  info_label_text_ = info_label;
  range_low_ = range_low;
  range_high_ = range_high;
  original_range_low_ = range_low;
  original_range_high_ = range_high;
  hist_x_min_ = hist->GetXaxis()->GetXmin();
  hist_x_max_ = hist->GetXaxis()->GetXmax();
  num_peaks_ = num_peaks;
  num_params_ = num_peaks * 10 + 2;

  accepted_ = kFALSE;
  done_ = kFALSE;
  needs_redraw_ = kFALSE;
  syncing_ = kFALSE;

  if (fit_func_->GetNpar() != num_params_) {
    std::cerr << "InteractiveFitEditor: expected " << num_params_
              << " parameters but TF1 has " << fit_func_->GetNpar()
              << std::endl;
  }

  // Save original state
  original_params_ = new Double_t[num_params_];
  original_bounds_low_ = new Double_t[num_params_];
  original_bounds_high_ = new Double_t[num_params_];
  original_fixed_ = new Bool_t[num_params_];
  current_bounds_low_ = new Double_t[num_params_];
  current_bounds_high_ = new Double_t[num_params_];

  for (Int_t i = 0; i < num_params_; i++) {
    original_params_[i] = fit_func_->GetParameter(i);

    Double_t lo = 0, hi = 0;
    fit_func_->GetParLimits(i, lo, hi);

    // lo >= hi means the parameter is fixed (FixParameter sets lo == hi)
    original_fixed_[i] = (lo >= hi);

    if (original_fixed_[i]) {
      // FixParameter overwrites limits; infer reasonable defaults
      GetDefaultBounds(i, lo, hi);
    }

    original_bounds_low_[i] = lo;
    original_bounds_high_[i] = hi;
    current_bounds_low_[i] = lo;
    current_bounds_high_[i] = hi;
  }

  // Allocate widget arrays
  sliders_ = new TGHSlider *[num_params_];
  value_entries_ = new TGNumberEntry *[num_params_];
  fix_checks_ = new TGCheckButton *[num_params_];
  lo_bound_entries_ = new TGNumberEntry *[num_params_];
  hi_bound_entries_ = new TGNumberEntry *[num_params_];

  hist_draw_ = nullptr;
  bkg_draw_ = nullptr;
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

  SetWindowName("Interactive Fit Editor");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  UpdateCanvas();
}

InteractiveFitEditor::~InteractiveFitEditor() {
  if (redraw_timer_) {
    redraw_timer_->TurnOff();
    delete redraw_timer_;
  }

  delete[] original_params_;
  delete[] original_bounds_low_;
  delete[] original_bounds_high_;
  delete[] original_fixed_;
  delete[] current_bounds_low_;
  delete[] current_bounds_high_;

  delete[] sliders_;
  delete[] value_entries_;
  delete[] fix_checks_;
  delete[] lo_bound_entries_;
  delete[] hi_bound_entries_;

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
  delete bkg_draw_;
  delete hist_draw_;
  delete residual_pad_;
  delete main_pad_;
}

// GUI Construction

void InteractiveFitEditor::BuildGUI() {
  // Main horizontal split: canvas (left) | controls (right)
  TGHorizontalFrame *main_frame = new TGHorizontalFrame(this, 1400, 850);
  AddFrame(main_frame,
           new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

  embedded_canvas_ =
      new TRootEmbeddedCanvas("FitEditorCanvas", main_frame, 850, 800);
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

  // Tab widget for peaks + background
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

void InteractiveFitEditor::BuildPeakTab(TGCompositeFrame *parent,
                                        Int_t peak_idx) {
  Int_t offset = peak_idx * 10;

  // Gaussian group
  TGGroupFrame *gaus_grp = new TGGroupFrame(parent, "Gaussian", kVerticalFrame);
  parent->AddFrame(gaus_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 1));
  AddParamRow(gaus_grp, offset + 0, "Mu");
  AddParamRow(gaus_grp, offset + 1, "Sigma");
  AddParamRow(gaus_grp, offset + 2, "Amplitude");

  // Step group
  TGGroupFrame *step_grp = new TGGroupFrame(parent, "Step", kVerticalFrame);
  parent->AddFrame(step_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(step_grp, offset + 3, "Step Amp");

  // Low tail group
  TGGroupFrame *ltail_grp =
      new TGGroupFrame(parent, "Low-Side Tails", kVerticalFrame);
  parent->AddFrame(ltail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 1));
  AddParamRow(ltail_grp, offset + 4, "Exp Amp");
  AddParamRow(ltail_grp, offset + 5, "Exp Decay/Sigma");
  AddParamRow(ltail_grp, offset + 6, "Lin Amp");
  AddParamRow(ltail_grp, offset + 7, "Lin Slope");

  // High tail group
  TGGroupFrame *htail_grp =
      new TGGroupFrame(parent, "High-Side Tail", kVerticalFrame);
  parent->AddFrame(htail_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 1, 3));
  AddParamRow(htail_grp, offset + 8, "Exp Amp");
  AddParamRow(htail_grp, offset + 9, "Exp Decay/Sigma");
}

void InteractiveFitEditor::BuildBackgroundTab(TGCompositeFrame *parent) {
  TGGroupFrame *bkg_grp =
      new TGGroupFrame(parent, "Background", kVerticalFrame);
  parent->AddFrame(bkg_grp, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
  AddParamRow(bkg_grp, BkgConstIdx(), "Constant");
  AddParamRow(bkg_grp, BkgSlopeIdx(), "Slope");
}

void InteractiveFitEditor::AddParamRow(TGCompositeFrame *parent,
                                       Int_t param_idx, const char *name) {
  TGHorizontalFrame *row = new TGHorizontalFrame(parent, 490, 28);
  parent->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

  Double_t val = fit_func_->GetParameter(param_idx);
  Bool_t fixed = original_fixed_[param_idx];

  // Label
  TGLabel *label = new TGLabel(row, name);
  row->AddFrame(label, new TGLayoutHints(kLHintsCenterY, 2, 4, 2, 2));
  label->Resize(90, 20);

  // Slider
  TGHSlider *slider =
      new TGHSlider(row, 140, kSlider1, kSliderBase + param_idx);
  slider->SetRange(0, kSliderRes);
  slider->SetPosition(ValToSlider(param_idx, val));
  slider->Associate(this);
  row->AddFrame(slider,
                new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 2, 2, 2));
  sliders_[param_idx] = slider;

  // Value entry
  TGNumberEntry *entry = new TGNumberEntry(
      row, val, 8, kEntryBase + param_idx, TGNumberFormat::kNESReal,
      TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
  entry->GetNumberEntry()->Associate(this);
  row->AddFrame(entry, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  value_entries_[param_idx] = entry;

  // Fix checkbox
  TGCheckButton *fix_cb = new TGCheckButton(row, "Fix", kFixBase + param_idx);
  fix_cb->Associate(this);
  if (fixed) {
    fix_cb->SetState(kButtonDown);
  }
  row->AddFrame(fix_cb, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
  fix_checks_[param_idx] = fix_cb;

  // Low bound entry
  TGNumberEntry *lo_entry = new TGNumberEntry(
      row, current_bounds_low_[param_idx], 6, kLoBoundBase + param_idx,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  lo_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(lo_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  lo_bound_entries_[param_idx] = lo_entry;

  // High bound entry
  TGNumberEntry *hi_entry = new TGNumberEntry(
      row, current_bounds_high_[param_idx], 6, kHiBoundBase + param_idx,
      TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
      TGNumberFormat::kNELNoLimits);
  hi_entry->GetNumberEntry()->Associate(this);
  row->AddFrame(hi_entry, new TGLayoutHints(kLHintsCenterY, 1, 1, 2, 2));
  hi_bound_entries_[param_idx] = hi_entry;

  // Disable slider and entry if parameter is fixed
  if (fixed) {
    slider->SetEnabled(kFALSE);
    entry->GetNumberEntry()->SetEnabled(kFALSE);
  }
}

// Drawing

void InteractiveFitEditor::InitDrawing() {
  TCanvas *canvas = embedded_canvas_->GetCanvas();
  canvas->cd();

  main_pad_ = new TPad("editor_main", "editor_main", 0, 0.3, 1, 1.0);
  main_pad_->SetBottomMargin(0.04);
  main_pad_->SetGridx(1);
  main_pad_->SetGridy(1);
  main_pad_->SetTopMargin(0.08);
  main_pad_->Draw();

  residual_pad_ = new TPad("editor_res", "editor_res", 0, 0, 1, 0.3);
  residual_pad_->SetTopMargin(0.04);
  residual_pad_->SetBottomMargin(0.35);
  residual_pad_->SetGridx(1);
  residual_pad_->SetGridy(1);
  residual_pad_->Draw();

  // Main pad
  main_pad_->cd();

  hist_draw_ = static_cast<TH1 *>(hist_->Clone("hist_editor_draw"));
  hist_draw_->SetDirectory(0);
  Double_t min_disp = 0.9 * range_low_;
  Double_t max_disp = 1.1 * range_high_;
  hist_draw_->GetXaxis()->SetRangeUser(min_disp, max_disp);
  hist_draw_->GetXaxis()->SetLabelSize(0);
  hist_draw_->GetXaxis()->SetTitleSize(0);
  hist_draw_->SetLineColor(kViolet);
  hist_draw_->SetLineWidth(2);
  hist_draw_->Draw();

  fit_func_->SetLineColor(kAzure);
  fit_func_->SetNpx(1000);
  fit_func_->Draw("same");

  bkg_draw_ = new TF1("bkg_editor", FittingFunctions::LinearBackground,
                      range_low_, range_high_, 2);
  bkg_draw_->SetParameter(0, fit_func_->GetParameter(BkgConstIdx()));
  bkg_draw_->SetParameter(1, fit_func_->GetParameter(BkgSlopeIdx()));
  bkg_draw_->SetLineColor(kGreen);
  bkg_draw_->SetNpx(1000);
  bkg_draw_->Draw("same");

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
  chi2_label_->SetTextAlign(31); // right-aligned
  chi2_label_->Draw();

  main_pad_->SetLogy(kTRUE);

  // Residual pad
  residual_pad_->cd();

  n_res_points_ = 0;
  Int_t nbins = hist_->GetNbinsX();
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = hist_->GetBinCenter(i);
    if (x < range_low_ || x > range_high_)
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

  zero_line_ = new TF1("zero_editor", "0", ax_min, ax_max);
  zero_line_->SetLineColor(kBlack);
  zero_line_->SetLineStyle(2);
  zero_line_->SetLineWidth(2);
  zero_line_->Draw("same");

  plus3_line_ = new TF1("plus3_editor", "3", ax_min, ax_max);
  plus3_line_->SetLineColor(kGray + 2);
  plus3_line_->SetLineStyle(3);
  plus3_line_->SetLineWidth(2);
  plus3_line_->Draw("same");

  minus3_line_ = new TF1("minus3_editor", "-3", ax_min, ax_max);
  minus3_line_->SetLineColor(kGray + 2);
  minus3_line_->SetLineStyle(3);
  minus3_line_->SetLineWidth(2);
  minus3_line_->Draw("same");

  UpdateCompPoints();
  UpdateResPoints();
}

void InteractiveFitEditor::UpdateCanvas() {
  UpdateCompPoints();
  UpdateResPoints();

  // Update background TF1 parameters
  bkg_draw_->SetParameter(0, fit_func_->GetParameter(BkgConstIdx()));
  bkg_draw_->SetParameter(1, fit_func_->GetParameter(BkgSlopeIdx()));

  // Compute live chi2/ndf
  Double_t chi2_sum = 0;
  Int_t ndf_bins = 0;
  Int_t nbins_chi2 = hist_->GetNbinsX();
  for (Int_t i = 1; i <= nbins_chi2; i++) {
    Double_t x = hist_->GetBinCenter(i);
    if (x < range_low_ || x > range_high_)
      continue;
    Double_t data = hist_->GetBinContent(i);
    Double_t error = hist_->GetBinError(i);
    if (error > 0 && data > 0) {
      Double_t fit_val = fit_func_->Eval(x);
      Double_t residual = (data - fit_val) / error;
      chi2_sum += residual * residual;
      ndf_bins++;
    }
  }
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

  // Update zero line to match current range
  zero_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);
  plus3_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);
  minus3_line_->SetRange(0.9 * range_low_, 1.1 * range_high_);

  // Lock residual Y range
  res_graph_->GetYaxis()->SetRangeUser(-5.5, 5.5);

  TCanvas *canvas = embedded_canvas_->GetCanvas();
  main_pad_->Modified();
  residual_pad_->Modified();
  canvas->Modified();
  canvas->Update();

  needs_redraw_ = kFALSE;
}

void InteractiveFitEditor::UpdateCompPoints() {
  Double_t x_step = (range_high_ - range_low_) / (kNDrawPts - 1);

  Double_t bkg_const = fit_func_->GetParameter(BkgConstIdx());
  Double_t bkg_slope = fit_func_->GetParameter(BkgSlopeIdx());

  for (Int_t p = 0; p < num_peaks_; p++) {
    Int_t off = p * 10;

    Double_t mu = fit_func_->GetParameter(off + 0);
    Double_t sigma = fit_func_->GetParameter(off + 1);
    Double_t gaus_amp = fit_func_->GetParameter(off + 2);
    // Amplitude params are ratios; convert to absolute for component drawing
    Double_t step_amp = fit_func_->GetParameter(off + 3) * gaus_amp;
    Double_t lexp_amp = fit_func_->GetParameter(off + 4) * gaus_amp;
    Double_t lexp_dec = fit_func_->GetParameter(off + 5);
    Double_t llin_amp = fit_func_->GetParameter(off + 6) * gaus_amp;
    Double_t llin_slp = fit_func_->GetParameter(off + 7);
    Double_t hexp_amp = fit_func_->GetParameter(off + 8) * gaus_amp;
    Double_t hexp_dec = fit_func_->GetParameter(off + 9);

    for (Int_t i = 0; i < kNDrawPts; i++) {
      Double_t x = range_low_ + i * x_step;
      Double_t x_arr[1] = {x};

      Double_t bkg_par[2] = {bkg_const, bkg_slope};
      Double_t bkg_val = FittingFunctions::LinearBackground(x_arr, bkg_par);

      Double_t gaus_par[3] = {mu, sigma, gaus_amp};
      comp_graphs_[p][0]->SetPoint(
          i, x, FittingFunctions::Gaussian(x_arr, gaus_par) + bkg_val);

      Double_t step_par[3] = {mu, sigma, step_amp};
      comp_graphs_[p][1]->SetPoint(
          i, x, FittingFunctions::Step(x_arr, step_par) + bkg_val);

      Double_t lt_par[6] = {mu, sigma, lexp_amp, lexp_dec, llin_amp, llin_slp};
      comp_graphs_[p][2]->SetPoint(
          i, x, FittingFunctions::LowTail(x_arr, lt_par) + bkg_val);

      Double_t ht_par[4] = {mu, sigma, hexp_amp, hexp_dec};
      comp_graphs_[p][3]->SetPoint(
          i, x, FittingFunctions::HighTail(x_arr, ht_par) + bkg_val);
    }
  }
}

void InteractiveFitEditor::UpdateResPoints() {
  Int_t nbins = hist_->GetNbinsX();
  Int_t pt = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = hist_->GetBinCenter(i);
    if (x < range_low_ || x > range_high_)
      continue;
    Double_t data = hist_->GetBinContent(i);
    Double_t error = hist_->GetBinError(i);
    if (error > 0 && data > 0) {
      Double_t fit_val = fit_func_->Eval(x);
      Double_t pull = (data - fit_val) / error;
      res_graph_->SetPoint(pt, x, pull);
      pt++;
    }
  }
  res_graph_->Set(pt);
  n_res_points_ = pt;

  // Update residual x-axis to match current range
  Double_t disp_lo = 0.9 * range_low_;
  Double_t disp_hi = 1.1 * range_high_;
  res_graph_->GetXaxis()->SetLimits(disp_lo, disp_hi);
}

// Widget Synchronization

void InteractiveFitEditor::SyncAllWidgets() {
  syncing_ = kTRUE;
  for (Int_t i = 0; i < num_params_; i++) {
    SyncWidget(i);
  }
  syncing_ = kFALSE;
}

void InteractiveFitEditor::SyncWidget(Int_t param_idx) {
  Double_t val = fit_func_->GetParameter(param_idx);

  // Update fixed state from current TF1
  Double_t lo = 0, hi = 0;
  fit_func_->GetParLimits(param_idx, lo, hi);
  Bool_t fixed = (lo >= hi);

  // Update bounds if parameter is not fixed
  if (!fixed && lo < hi) {
    current_bounds_low_[param_idx] = lo;
    current_bounds_high_[param_idx] = hi;
  }

  // Clamp value to bounds for slider
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

// Event Handlers

void InteractiveFitEditor::OnSliderMoved(Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(param_idx))
    return;

  Int_t pos = sliders_[param_idx]->GetPosition();
  Double_t val = SliderToVal(param_idx, pos);
  fit_func_->SetParameter(param_idx, val);

  syncing_ = kTRUE;
  value_entries_[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveFitEditor::OnEntryChanged(Int_t param_idx) {
  if (syncing_)
    return;
  if (IsFixed(param_idx))
    return;

  Double_t val = value_entries_[param_idx]->GetNumber();

  // Expand bounds if the entered value falls outside
  Bool_t bounds_changed = kFALSE;
  if (val < current_bounds_low_[param_idx]) {
    current_bounds_low_[param_idx] = val;
    bounds_changed = kTRUE;
  }
  if (val > current_bounds_high_[param_idx]) {
    current_bounds_high_[param_idx] = val;
    bounds_changed = kTRUE;
  }

  fit_func_->SetParameter(param_idx, val);
  fit_func_->SetParLimits(param_idx, current_bounds_low_[param_idx],
                          current_bounds_high_[param_idx]);

  syncing_ = kTRUE;
  sliders_[param_idx]->SetPosition(ValToSlider(param_idx, val));
  if (bounds_changed) {
    lo_bound_entries_[param_idx]->SetNumber(current_bounds_low_[param_idx]);
    hi_bound_entries_[param_idx]->SetNumber(current_bounds_high_[param_idx]);
  }
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveFitEditor::OnBoundsChanged(Int_t param_idx) {
  if (syncing_)
    return;

  Double_t new_lo = lo_bound_entries_[param_idx]->GetNumber();
  Double_t new_hi = hi_bound_entries_[param_idx]->GetNumber();

  if (new_lo >= new_hi)
    return; // invalid bounds, ignore

  current_bounds_low_[param_idx] = new_lo;
  current_bounds_high_[param_idx] = new_hi;

  // Clamp current value to new bounds
  Double_t val = fit_func_->GetParameter(param_idx);
  if (val < new_lo)
    val = new_lo;
  if (val > new_hi)
    val = new_hi;
  fit_func_->SetParameter(param_idx, val);

  if (!IsFixed(param_idx)) {
    fit_func_->SetParLimits(param_idx, new_lo, new_hi);
  }

  syncing_ = kTRUE;
  sliders_[param_idx]->SetPosition(ValToSlider(param_idx, val));
  value_entries_[param_idx]->SetNumber(val);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

void InteractiveFitEditor::OnFixToggled(Int_t param_idx) {
  Bool_t now_fixed = (fix_checks_[param_idx]->GetState() == kButtonDown);

  if (now_fixed) {
    Double_t val = fit_func_->GetParameter(param_idx);
    fit_func_->FixParameter(param_idx, val);
    sliders_[param_idx]->SetEnabled(kFALSE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kFALSE);
  } else {
    // Restore bounds and release
    fit_func_->SetParLimits(param_idx, current_bounds_low_[param_idx],
                            current_bounds_high_[param_idx]);
    sliders_[param_idx]->SetEnabled(kTRUE);
    value_entries_[param_idx]->GetNumberEntry()->SetEnabled(kTRUE);
  }

  needs_redraw_ = kTRUE;
}

void InteractiveFitEditor::OnRangeChanged() {
  if (syncing_)
    return;

  Double_t new_lo = range_slider_->GetMinPosition();
  Double_t new_hi = range_slider_->GetMaxPosition();

  if (new_lo >= new_hi)
    return;

  range_low_ = new_lo;
  range_high_ = new_hi;

  fit_func_->SetRange(range_low_, range_high_);
  bkg_draw_->SetRange(range_low_, range_high_);

  // Update histogram display range
  hist_draw_->GetXaxis()->SetRangeUser(0.9 * range_low_, 1.1 * range_high_);

  // Sync number entries
  syncing_ = kTRUE;
  range_lo_entry_->SetNumber(range_low_);
  range_hi_entry_->SetNumber(range_high_);
  syncing_ = kFALSE;

  needs_redraw_ = kTRUE;
}

// Actions

void InteractiveFitEditor::DoRefit() {
  // Release all non-fixed parameters and use current values as initial guesses
  for (Int_t i = 0; i < num_params_; i++) {
    if (!IsFixed(i)) {
      fit_func_->SetParLimits(i, current_bounds_low_[i],
                              current_bounds_high_[i]);
    }
  }

  hist_->Fit(fit_func_, "LSMRBEN");

  SyncAllWidgets();
  needs_redraw_ = kTRUE;
  UpdateCanvas();
}

void InteractiveFitEditor::DoAccept() {
  accepted_ = kTRUE;
  done_ = kTRUE;
}

void InteractiveFitEditor::DoCancel() {
  // Restore original parameters and fixed states
  for (Int_t i = 0; i < num_params_; i++) {
    fit_func_->SetParameter(i, original_params_[i]);
    if (original_fixed_[i]) {
      fit_func_->FixParameter(i, original_params_[i]);
    } else {
      fit_func_->SetParLimits(i, original_bounds_low_[i],
                              original_bounds_high_[i]);
    }
  }
  // Restore original range
  fit_func_->SetRange(original_range_low_, original_range_high_);
  accepted_ = kFALSE;
  done_ = kTRUE;
}

void InteractiveFitEditor::DoReset() {
  // Restore original parameters but don't close
  for (Int_t i = 0; i < num_params_; i++) {
    fit_func_->SetParameter(i, original_params_[i]);
    current_bounds_low_[i] = original_bounds_low_[i];
    current_bounds_high_[i] = original_bounds_high_[i];
    if (original_fixed_[i]) {
      fit_func_->FixParameter(i, original_params_[i]);
    } else {
      fit_func_->SetParLimits(i, original_bounds_low_[i],
                              original_bounds_high_[i]);
    }
  }

  // Restore range
  range_low_ = original_range_low_;
  range_high_ = original_range_high_;
  fit_func_->SetRange(range_low_, range_high_);
  bkg_draw_->SetRange(range_low_, range_high_);
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

// Event Loop

Bool_t InteractiveFitEditor::ProcessMessage(Long_t msg, Long_t parm1,
                                            Long_t parm2) {
  if (syncing_)
    return kTRUE;

  (void)parm2; // unused directly; slider position read from widget

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
        // Sync slider from entries
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

Bool_t InteractiveFitEditor::HandleTimer(TTimer *timer) {
  if (timer == redraw_timer_ && needs_redraw_) {
    UpdateCanvas();
  }
  return kTRUE;
}

void InteractiveFitEditor::CloseWindow() { DoCancel(); }

// Helpers

Int_t InteractiveFitEditor::ValToSlider(Int_t param_idx, Double_t val) {
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

Double_t InteractiveFitEditor::SliderToVal(Int_t param_idx, Int_t pos) {
  Double_t lo = current_bounds_low_[param_idx];
  Double_t hi = current_bounds_high_[param_idx];
  Double_t frac = static_cast<Double_t>(pos) / kSliderRes;
  return lo + frac * (hi - lo);
}

Bool_t InteractiveFitEditor::IsFixed(Int_t param_idx) {
  return (fix_checks_[param_idx]->GetState() == kButtonDown);
}

void InteractiveFitEditor::GetDefaultBounds(Int_t param_idx, Double_t &lo,
                                            Double_t &hi) {
  Double_t peak_height = hist_->GetMaximum();
  Double_t range_width = range_high_ - range_low_;

  if (param_idx == BkgConstIdx()) {
    lo = 0;
    hi = peak_height;
    return;
  }
  if (param_idx == BkgSlopeIdx()) {
    lo = -peak_height / range_width;
    hi = peak_height / range_width;
    return;
  }

  Int_t local = param_idx % 10;
  switch (local) {
  case 0: // Mu
    lo = range_low_;
    hi = range_high_;
    break;
  case 1: // Sigma
    lo = range_width * 0.001;
    hi = range_width * 0.5;
    break;
  case 2: // GausAmplitude
    lo = 0;
    hi = peak_height * 2;
    break;
  case 3: // StepAmplitude ratio
    lo = 0;
    hi = 0.5;
    break;
  case 4: // LowExpTailAmplitude ratio
    lo = 0;
    hi = 0.5;
    break;
  case 5: // LowExpTailRatio
    lo = 1.0;
    hi = 100;
    break;
  case 6: // LowLinTailAmplitude ratio
    lo = 0;
    hi = 0.5;
    break;
  case 7: // LowLinTailSlope
    lo = -1;
    hi = 1;
    break;
  case 8: // HighExpTailAmplitude ratio
    lo = 0;
    hi = 0.5;
    break;
  case 9: // HighExpTailRatio
    lo = 1.0;
    hi = 100;
    break;
  default:
    lo = -1000;
    hi = 1000;
    break;
  }
}

Int_t InteractiveFitEditor::PeakStyle(Int_t peak_idx) {
  if (peak_idx == 0)
    return 1; // solid
  if (peak_idx == 1)
    return 3; // dashed
  return 4;   // dash-dot
}

// Launcher

Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func, Double_t range_low,
                                  Double_t range_high, Int_t num_peaks,
                                  const TString &info_label) {
  if (!gClient) {
    std::cerr << "InteractiveFitEditor: GUI not available (gClient is null). "
              << "Make sure you are not in batch mode." << std::endl;
    return kFALSE;
  }

  // Swallow recoverable X protocol errors for the editor's lifetime so a
  // transient bad redraw doesn't trip ROOT's crashing default handler.
  AUXErrorHandlerSave xerr_save = AUInstallTolerantXErrorHandler();

  InteractiveFitEditor *editor =
      new InteractiveFitEditor(gClient->GetRoot(), hist, fit_func, range_low,
                               range_high, num_peaks, info_label);

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

  // DontCallClose prevents TGMainFrame from triggering gApplication->Terminate
  // when the window is destroyed, which would kill the whole macro
  editor->DontCallClose();
  editor->UnmapWindow();
  gSystem->ProcessEvents();
  delete editor;

  AURestoreXErrorHandler(xerr_save);
  return result;
}
