#include "PlottingUtils.hpp"

PlotSaveFormat PlottingUtils::save_format_ = PlotSaveFormat::kPNG;
Bool_t PlottingUtils::preferences_set_ = kFALSE;
Width_t PlottingUtils::line_width_ = 2;
TString PlottingUtils::plots_base_dir_ = "plots";

void PlottingUtils::SetPlotsBaseDir(const TString &dir) {
  TString d = dir;
  while (d.Length() > 0 && d[d.Length() - 1] == '/')
    d.Chop();
  plots_base_dir_ = d;
}

TString PlottingUtils::GetPlotsBaseDir() { return plots_base_dir_; }

void PlottingUtils::WarnIfNotConfigured(const TString method_name) {
  if (!preferences_set_)
    std::cout << "WARNING: " << method_name
              << " called before SetStylePreferences()." << std::endl;
}

void PlottingUtils::SetStylePreferences(PlotSaveFormat save_format) {
  save_format_ = save_format;
  preferences_set_ = kTRUE;
  line_width_ = save_format_ == PlotSaveFormat::kPNG ? 2 : 1;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.06, "XY");
  gStyle->SetLabelSize(0.06, "XY");
  gStyle->SetLegendFont(132);
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTextFont(42);
  gStyle->SetHistLineWidth(line_width_);
  gStyle->SetLineWidth(line_width_);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(line_width_);
  gStyle->SetGridColor(kGray);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void PlottingUtils::ConfigureGraph(TGraph *graph, Int_t color,
                                   const TString title) {
  WarnIfNotConfigured("ConfigureGraph");
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetLineWidth(line_width_);
}

void PlottingUtils::ConfigureGraph(TGraphErrors *graph, Int_t color,
                                   const TString title) {
  WarnIfNotConfigured("ConfigureGraph");
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(line_width_);
}

void PlottingUtils::ConfigureHistogram(TH1 *hist, Int_t color,
                                       const TString title) {
  if (!hist)
    return;
  WarnIfNotConfigured("ConfigureHistogram");

  hist->SetLineColor(color);
  hist->SetTitle(title);
  hist->SetLineWidth(line_width_);
  hist->SetFillStyle(0);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetNoExponent(kTRUE);
  hist->GetYaxis()->SetNdivisions(50109);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.2);
}

void PlottingUtils::Configure2DHistogram(TH2 *hist, TCanvas *canvas,
                                         const TString title) {
  if (!hist)
    return;
  if (!canvas)
    return;
  WarnIfNotConfigured("Configure2DHistogram");

  hist->SetTitle(title);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1);
  hist->GetXaxis()->SetNdivisions(506);
  hist->GetYaxis()->SetNdivisions(506);

  canvas->SetLogz(kTRUE);
  canvas->SetRightMargin(0.15);
}

void PlottingUtils::ConfigureAndDrawGraph(TGraph *graph, Int_t color,
                                          const TString title) {
  if (!graph)
    return;

  ConfigureGraph(graph, color, title);
  graph->Draw();
}

void PlottingUtils::ConfigureAndDrawHistogram(TH1 *hist, Int_t color,
                                              const TString title) {
  if (!hist)
    return;

  ConfigureHistogram(hist, color, title);
  hist->Draw("HIST");
}

void PlottingUtils::ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                                const TString title) {
  if (!hist)
    return;
  if (!canvas)
    return;

  Configure2DHistogram(hist, canvas, title);
  hist->Draw("COLZ");
}

TCanvas *PlottingUtils::GetConfiguredCanvas(Bool_t logy) {
  WarnIfNotConfigured("GetConfiguredCanvas");
  TCanvas *canvas = new TCanvas(GetRandomName(), "", 1200, 800);

  canvas->SetGridx(1);
  canvas->SetGridy(1);
  canvas->SetLogy(logy);

  canvas->SetTicks(1, 1);
  gPad->SetTicks(1, 1);

  return canvas;
}

void PlottingUtils::SaveFigure(TCanvas *canvas, TString output_name,
                               TString output_subdirectory,
                               PlotSaveOptions save_options) {
  WarnIfNotConfigured("SaveFigure");
  canvas->SetLogy(kFALSE);
  canvas->Modified();
  canvas->Update();

  TString extension = (save_format_ == PlotSaveFormat::kPNG) ? ".png" : ".pdf";
  TString output_filename = output_name + extension;

  TString base = GetPlotsBaseDir();
  TString full_dir =
      output_subdirectory == "" ? base : base + "/" + output_subdirectory;

  if (gSystem->AccessPathName(full_dir)) {
    gSystem->mkdir(full_dir, kTRUE);
  }

  TString prefix = full_dir + "/";

  if (save_options != PlotSaveOptions::kLOG)
    canvas->Print(prefix + output_filename);

  if (save_options != PlotSaveOptions::kLINEAR) {
    TList *primitives = canvas->GetListOfPrimitives();
    Int_t size = primitives->GetSize();

    for (Int_t i = 0; i < size; i++) {
      TObject *object = primitives->At(i);
      if (object->InheritsFrom(TH2::Class())) {
        std::cout << std::endl;
        std::cerr << "ERROR: Used PlotSaveOptions::kLOG for 2D histogram."
                  << std::endl;
        std::cout
            << "This option is exclusive to 1D histograms/graphs because it refers to the y axis."
            << std::endl;
        std::cout << "Use PlotSaveOptions::kLINEAR." << std::endl;
        std::cout
            << "The z axis is already log by default if you used PlottingUtils::Configure2DHistogram()."
            << std::endl;
        std::cout << "Plot " << prefix + "log_" + output_filename
                  << " was not saved." << std::endl;
        std::exit(1);
      };
    };

    canvas->SetLogy(kTRUE);
    canvas->Modified();
    canvas->Update();
    canvas->Print(prefix + "log_" + output_filename);

    canvas->SetLogy(kFALSE);
    canvas->Modified();
    canvas->Update();
  }
}

std::vector<Int_t> PlottingUtils::GetDefaultColors() {
  return {kRed + 1,   kBlue + 1,   kGreen + 2,  kOrange + 1,  kMagenta + 1,
          kCyan + 2,  kViolet + 1, kSpring - 1, kPink + 1,    kTeal + 2,
          kAzure + 2, kYellow + 1, kOrange - 3, kMagenta - 3, kCyan - 6,
          kRed - 4,   kBlue - 4,   kGreen - 6,  kViolet - 4,  kSpring + 5,
          kPink - 3,  kTeal - 5,   kAzure - 3,  kOrange + 7};
}

TLegend *PlottingUtils::AddLegend(Double_t x1, Double_t x2, Double_t y1,
                                  Double_t y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(30);
  leg->SetTextFont(43);
  leg->Draw();

  return leg;
}

TLatex *PlottingUtils::AddText(const TString label, Double_t x, Double_t y,
                               Double_t angle) {
  TLatex *text = new TLatex(x, y, label);
  text->SetNDC();
  text->SetTextSize(30);
  text->SetTextAlign(33);
  text->SetTextFont(43);
  text->SetTextAngle(angle);
  text->Draw();

  return text;
}

TString PlottingUtils::GetRandomName() {
  static TRandom3 generator(0);
  Double_t number = generator.Rndm();
  TString name = Form("name%.7f", number);
  return name;
}

void PlottingUtils::PlotFitWithResiduals(
    TH1 *hist, TGraph *total_graph,
    const std::vector<TGraph *> &component_graphs, Float_t fit_range_low,
    Float_t fit_range_high, const TString &output_name,
    const TString &output_subdirectory, const TString &label, Bool_t logy) {
  TCanvas *canvas = GetConfiguredCanvas(kFALSE);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad1->SetBottomMargin(0.04);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low;
  Float_t max_hist_value = 1.1 * fit_range_high;

  hist->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  hist->GetXaxis()->SetLabelSize(0);
  hist->GetXaxis()->SetTitleSize(0);
  hist->SetLineColor(kViolet);
  hist->GetYaxis()->SetTitleOffset(1);
  hist->SetLineWidth(line_width_);
  hist->Draw();

  pad1->SetTickx(0);
  if (total_graph)
    total_graph->Draw("L same");
  for (std::size_t i = 0; i < component_graphs.size(); i++) {
    if (component_graphs[i])
      component_graphs[i]->Draw("L same");
  }

  pad2->cd();

  Int_t nbins = hist->GetNbinsX();
  TGraph *residuals = new TGraph();
  Int_t point_counter = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = hist->GetBinCenter(i);
    if (x < fit_range_low || x > fit_range_high)
      continue;
    Double_t data = hist->GetBinContent(i);
    Double_t fit_val = total_graph ? total_graph->Eval(x) : 0.0;
    Double_t error = hist->GetBinError(i);

    if (error > 0 && data > 0) {
      Double_t pull = (data - fit_val) / error;
      residuals->SetPoint(point_counter, x, pull);
      point_counter++;
    }
  }

  residuals->SetMarkerStyle(20);
  residuals->SetMarkerSize(0.8);
  residuals->SetMarkerColor(kAzure);
  residuals->SetLineColor(kAzure);
  residuals->SetTitle("");
  Double_t actual_min =
      hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetFirst());
  Double_t actual_max =
      hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());
  residuals->GetXaxis()->SetLimits(actual_min, actual_max);
  residuals->GetYaxis()->SetTitle("#delta/#sigma");
  residuals->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  residuals->GetXaxis()->SetTitleSize(0.13);
  residuals->GetYaxis()->SetTitleSize(0.13);
  residuals->GetXaxis()->SetLabelSize(0.12);
  residuals->GetYaxis()->SetLabelSize(0.12);
  residuals->GetXaxis()->SetTitleOffset(1.0);
  residuals->GetYaxis()->SetTitleOffset(0.3);
  residuals->GetYaxis()->SetNdivisions(505);
  residuals->GetXaxis()->SetNdivisions(510);
  residuals->GetYaxis()->CenterTitle(kTRUE);
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(line_width_);
  zero_line->Draw("same");

  TF1 *plus3_line = new TF1("plus3_line", "3", actual_min, actual_max);
  plus3_line->SetLineColor(kGray + 2);
  plus3_line->SetLineStyle(3);
  plus3_line->SetLineWidth(line_width_);
  plus3_line->Draw("same");

  TF1 *minus3_line = new TF1("minus3_line", "-3", actual_min, actual_max);
  minus3_line->SetLineColor(kGray + 2);
  minus3_line->SetLineStyle(3);
  minus3_line->SetLineWidth(line_width_);
  minus3_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(logy);
  if (label.Length() > 0) {
    AddText(label, 0.85, 0.85);
  }

  PlotSaveOptions save_opts =
      logy ? PlotSaveOptions::kLOG : PlotSaveOptions::kLINEAR;
  SaveFigure(canvas, output_name, output_subdirectory, save_opts);

  PlotPullHistogram(residuals, output_name, output_subdirectory);

  gROOT->GetListOfCanvases()->Remove(canvas);
  canvas->SetBatch(kTRUE);
  delete canvas;
  delete residuals;
  delete zero_line;
  delete plus3_line;
  delete minus3_line;
}

void PlottingUtils::PlotPullHistogram(TGraph *residuals,
                                      const TString &output_name,
                                      const TString &output_subdirectory) {
  Int_t npoints = residuals->GetN();
  if (npoints == 0)
    return;

  TH1D *pull_hist =
      new TH1D("pull_hist", ";#delta/#sigma;Counts", 82, -5.5, 5.5);

  Double_t *y = residuals->GetY();
  for (Int_t i = 0; i < npoints; i++) {
    pull_hist->Fill(y[i]);
  }

  TH1D *gauss_ref = new TH1D("gauss_ref", "", 82, -5.5, 5.5);
  for (Int_t i = 1; i <= gauss_ref->GetNbinsX(); i++) {
    Double_t lo = gauss_ref->GetBinLowEdge(i);
    Double_t hi = lo + gauss_ref->GetBinWidth(i);
    gauss_ref->SetBinContent(i, npoints * (TMath::Freq(hi) - TMath::Freq(lo)));
  }

  Double_t ks_pvalue = pull_hist->KolmogorovTest(gauss_ref);

  TCanvas *hist_canvas = GetConfiguredCanvas(kFALSE);
  ConfigureAndDrawHistogram(pull_hist, kAzure);
  AddText(TString::Format("KS p = %.3f", ks_pvalue), 0.85, 0.85);

  SaveFigure(hist_canvas, "residuals_" + output_name,
             output_subdirectory + "/residual_hists", PlotSaveOptions::kLINEAR);

  delete hist_canvas;
  delete gauss_ref;
  delete pull_hist;
}
