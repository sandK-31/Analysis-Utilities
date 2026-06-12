#ifndef PLOTTINGUTILS_H
#define PLOTTINGUTILS_H

#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

enum class PlotSaveOptions { kLINEAR, kLOG, kBOTH };
enum class PlotSaveFormat { kPNG, kPDF };

class PlottingUtils {
public:
  static void
  SetStylePreferences(PlotSaveFormat save_format = PlotSaveFormat::kPNG);
  static void ConfigureAndDrawGraph(TGraph *graph, Int_t color,
                                    const TString title = "");
  static void ConfigureAndDrawHistogram(TH1 *hist, Int_t color,
                                        const TString title = "");
  static void ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                          const TString title = "");

  static void ConfigureGraph(TGraph *graph, Int_t color,
                             const TString title = "");
  static void ConfigureGraph(TGraphErrors *graph, Int_t color,
                             const TString title = "");

  static void ConfigureHistogram(TH1 *hist, Int_t color,
                                 const TString title = "");
  static void Configure2DHistogram(TH2 *hist, TCanvas *canvas,
                                   const TString title = "");

  static TCanvas *GetConfiguredCanvas(Bool_t logy = kFALSE);
  static void SaveFigure(TCanvas *canvas, TString output_filename,
                         TString output_subdirectory = "",
                         PlotSaveOptions save_options = PlotSaveOptions::kBOTH);

  static void SetPlotsBaseDir(const TString &dir);
  static TString GetPlotsBaseDir();

  static void
  PlotFitWithResiduals(TH1 *hist, TGraph *total_graph,
                       const std::vector<TGraph *> &component_graphs,
                       Float_t fit_range_low, Float_t fit_range_high,
                       const TString &output_name,
                       const TString &output_subdirectory = "fits",
                       const TString &label = "", Bool_t logy = kTRUE);

  static TLegend *AddLegend(Double_t x1 = 0.7, Double_t x2 = 0.9,
                            Double_t y1 = 0.7, Double_t y2 = 0.9);
  static TLatex *AddText(const TString label, Double_t x = 0.9,
                         Double_t y = 0.85, Double_t angle = 0);

  static std::vector<Int_t> GetDefaultColors();

  static TString GetRandomName();
  static Width_t GetLineWidth() { return line_width_; };

private:
  static PlotSaveFormat save_format_;
  static Bool_t preferences_set_;
  static Width_t line_width_;
  static TString plots_base_dir_;
  static void WarnIfNotConfigured(const TString method_name);
  static void PlotPullHistogram(TGraph *residuals, const TString &output_name,
                                const TString &output_subdirectory);
};

#endif
