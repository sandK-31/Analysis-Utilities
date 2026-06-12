#include "FittingUtils.hpp"

void FittingUtils::SaveInteractiveParams(const TString &input_name,
                                         const TString &peak_name) {
  TString fits_dir = PlottingUtils::GetPlotsBaseDir() + "/fits";
  gSystem->mkdir(fits_dir, kTRUE);
  TString filename = fits_dir + "/" + peak_name + "_" + input_name + ".fits";
  std::ofstream out(filename.Data());
  if (!out.is_open()) {
    std::cerr << "WARNING: Could not save interactive params to " << filename
              << std::endl;
    return;
  }
  Int_t npar = fit_function_->GetNpar();
  out << std::setprecision(15);
  out << "RANGE " << fit_range_low_ << " " << fit_range_high_ << "\n";
  for (Int_t i = 0; i < npar; i++) {
    Double_t lo = 0, hi = 0;
    fit_function_->GetParLimits(i, lo, hi);
    Bool_t fixed = (lo >= hi);
    out << fit_function_->GetParName(i) << " " << fit_function_->GetParameter(i)
        << " " << (fixed ? 1 : 0) << "\n";
  }
  out.close();
  std::cout << "Saved interactive params to " << filename << std::endl;
}

Bool_t FittingUtils::LoadInteractiveParams(const TString &input_name,
                                           const TString &peak_name) {
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + peak_name +
                     "_" + input_name + ".fits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return kFALSE;

  Int_t npar = fit_function_->GetNpar();
  std::string token;
  Int_t idx = 0;

  // First line: RANGE low high
  in >> token;
  if (token == "RANGE") {
    Double_t rlo, rhi;
    in >> rlo >> rhi;
    fit_range_low_ = rlo;
    fit_range_high_ = rhi;
    fit_function_->SetRange(rlo, rhi);
  }

  // Remaining lines: name value fixed
  Double_t value;
  Int_t fixed;
  while (in >> token >> value >> fixed && idx < npar) {
    fit_function_->SetParameter(idx, value);
    if (fixed)
      fit_function_->FixParameter(idx, value);
    idx++;
  }
  in.close();

  if (idx != npar) {
    std::cerr << "WARNING: Parameter count mismatch in " << filename
              << " (expected " << npar << ", got " << idx << ")" << std::endl;
    return kFALSE;
  }

  std::cout << "Loaded interactive params from " << filename << std::endl;
  return kTRUE;
}

Double_t FittingFunctions::Gaussian(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x[0] - mu) / sigma;
  Double_t gaus_amplitude = par[2];
  return gaus_amplitude * TMath::Exp(-0.5 * z * z);
}

Double_t FittingFunctions::LinearBackground(Double_t *x, Double_t *par) {
  Double_t bkg_constant = par[0];
  Double_t lin_bkg_slope = par[1];
  return lin_bkg_slope * x[0] + bkg_constant;
}

Double_t FittingFunctions::LowTail(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t exp_tail_amplitude = par[2];
  Double_t exp_tail_ratio = par[3];
  Double_t lin_tail_amplitude = par[4];
  Double_t lin_tail_slope = par[5];

  if (sigma <= 0)
    return 0;

  if (exp_tail_amplitude == 0 && lin_tail_amplitude == 0)
    return 0;

  Double_t y = x[0] - mu;
  Double_t tau = exp_tail_ratio * sigma;

  Double_t exp_term = (exp_tail_amplitude == 0 || tau <= 0)
                          ? 0
                          : exp_tail_amplitude * TMath::Exp(y / tau);
  Double_t lin_term = lin_tail_amplitude == 0
                          ? 0
                          : lin_tail_amplitude * (1.0 + lin_tail_slope * y);
  Double_t erfc_term = 1.0 - TMath::Erf(y / (TMath::Sqrt(2) * sigma));

  return (exp_term + lin_term) * erfc_term;
}

Double_t FittingFunctions::HighTail(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t exp_tail_amplitude = par[2];
  Double_t exp_tail_ratio = par[3];

  if (sigma <= 0)
    return 0;

  if (exp_tail_amplitude == 0)
    return 0;

  Double_t y = mu - x[0];
  Double_t tau = exp_tail_ratio * sigma;
  if (tau <= 0)
    return 0;

  Double_t exp_term = exp_tail_amplitude * TMath::Exp(y / tau);
  Double_t erfc_term = 1.0 - TMath::Erf(y / (TMath::Sqrt(2) * sigma));

  return exp_term * erfc_term;
}

Double_t FittingFunctions::Step(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  if (sigma <= 0)
    return 0;

  Double_t z = (x[0] - mu) / sigma;
  Double_t step_amplitude = par[2];

  Double_t denominator = TMath::Power(1 + TMath::Exp(z), 2);
  if (denominator < 1e-100)
    return 0;

  return step_amplitude / denominator;
}

Double_t FittingFunctions::PeakFunction(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t gaus_amplitude = par[2];
  // par[3,4,6,8] are ratios of gaus_amplitude
  Double_t step_amplitude = par[3] * gaus_amplitude;
  Double_t low_exp_tail_amplitude = par[4] * gaus_amplitude;
  Double_t low_exp_tail_ratio = par[5];
  Double_t low_lin_tail_amplitude = par[6] * gaus_amplitude;
  Double_t low_lin_tail_slope = par[7];
  Double_t high_exp_tail_amplitude = par[8] * gaus_amplitude;
  Double_t high_exp_tail_ratio = par[9];
  Double_t bkg_constant = par[10];
  Double_t lin_bkg_slope = par[11];

  Double_t gaus_par[3] = {mu, sigma, gaus_amplitude};
  Double_t bkg_par[2] = {bkg_constant, lin_bkg_slope};
  Double_t step_par[3] = {mu, sigma, step_amplitude};
  Double_t low_tail_par[6] = {mu,
                              sigma,
                              low_exp_tail_amplitude,
                              low_exp_tail_ratio,
                              low_lin_tail_amplitude,
                              low_lin_tail_slope};
  Double_t high_tail_par[4] = {mu, sigma, high_exp_tail_amplitude,
                               high_exp_tail_ratio};

  return Gaussian(x, gaus_par) + LinearBackground(x, bkg_par) +
         Step(x, step_par) + LowTail(x, low_tail_par) +
         HighTail(x, high_tail_par);
}

Double_t FittingFunctions::DoublePeakFunction(Double_t *x, Double_t *par) {
  // Peak 1: params 0-9
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t step_amplitude1 = par[3] * gaus_amplitude1;
  Double_t low_exp_tail_amplitude1 = par[4] * gaus_amplitude1;
  Double_t low_exp_tail_ratio1 = par[5];
  Double_t low_lin_tail_amplitude1 = par[6] * gaus_amplitude1;
  Double_t low_lin_tail_slope1 = par[7];
  Double_t high_exp_tail_amplitude1 = par[8] * gaus_amplitude1;
  Double_t high_exp_tail_ratio1 = par[9];

  // Peak 2: params 10-19
  Double_t mu2 = par[10];
  Double_t sigma2 = par[11];
  Double_t gaus_amplitude2 = par[12];
  Double_t step_amplitude2 = par[13] * gaus_amplitude2;
  Double_t low_exp_tail_amplitude2 = par[14] * gaus_amplitude2;
  Double_t low_exp_tail_ratio2 = par[15];
  Double_t low_lin_tail_amplitude2 = par[16] * gaus_amplitude2;
  Double_t low_lin_tail_slope2 = par[17];
  Double_t high_exp_tail_amplitude2 = par[18] * gaus_amplitude2;
  Double_t high_exp_tail_ratio2 = par[19];

  // Background: params 20-21
  Double_t bkg_const = par[20];
  Double_t bkg_slope = par[21];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t step1_par[3] = {mu1, sigma1, step_amplitude1};
  Double_t low_tail1_par[6] = {mu1,
                               sigma1,
                               low_exp_tail_amplitude1,
                               low_exp_tail_ratio1,
                               low_lin_tail_amplitude1,
                               low_lin_tail_slope1};
  Double_t high_tail1_par[4] = {mu1, sigma1, high_exp_tail_amplitude1,
                                high_exp_tail_ratio1};

  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t step2_par[3] = {mu2, sigma2, step_amplitude2};
  Double_t low_tail2_par[6] = {mu2,
                               sigma2,
                               low_exp_tail_amplitude2,
                               low_exp_tail_ratio2,
                               low_lin_tail_amplitude2,
                               low_lin_tail_slope2};
  Double_t high_tail2_par[4] = {mu2, sigma2, high_exp_tail_amplitude2,
                                high_exp_tail_ratio2};

  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Step(x, step1_par) +
         LowTail(x, low_tail1_par) + HighTail(x, high_tail1_par) +
         Gaussian(x, gaus2_par) + Step(x, step2_par) +
         LowTail(x, low_tail2_par) + HighTail(x, high_tail2_par) +
         LinearBackground(x, bkg_par);
}

Double_t FittingFunctions::TriplePeakFunction(Double_t *x, Double_t *par) {
  // Peak 1: params 0-9
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t step_amplitude1 = par[3] * gaus_amplitude1;
  Double_t low_exp_tail_amplitude1 = par[4] * gaus_amplitude1;
  Double_t low_exp_tail_ratio1 = par[5];
  Double_t low_lin_tail_amplitude1 = par[6] * gaus_amplitude1;
  Double_t low_lin_tail_slope1 = par[7];
  Double_t high_exp_tail_amplitude1 = par[8] * gaus_amplitude1;
  Double_t high_exp_tail_ratio1 = par[9];

  // Peak 2: params 10-19
  Double_t mu2 = par[10];
  Double_t sigma2 = par[11];
  Double_t gaus_amplitude2 = par[12];
  Double_t step_amplitude2 = par[13] * gaus_amplitude2;
  Double_t low_exp_tail_amplitude2 = par[14] * gaus_amplitude2;
  Double_t low_exp_tail_ratio2 = par[15];
  Double_t low_lin_tail_amplitude2 = par[16] * gaus_amplitude2;
  Double_t low_lin_tail_slope2 = par[17];
  Double_t high_exp_tail_amplitude2 = par[18] * gaus_amplitude2;
  Double_t high_exp_tail_ratio2 = par[19];

  // Peak 3: params 20-29
  Double_t mu3 = par[20];
  Double_t sigma3 = par[21];
  Double_t gaus_amplitude3 = par[22];
  Double_t step_amplitude3 = par[23] * gaus_amplitude3;
  Double_t low_exp_tail_amplitude3 = par[24] * gaus_amplitude3;
  Double_t low_exp_tail_ratio3 = par[25];
  Double_t low_lin_tail_amplitude3 = par[26] * gaus_amplitude3;
  Double_t low_lin_tail_slope3 = par[27];
  Double_t high_exp_tail_amplitude3 = par[28] * gaus_amplitude3;
  Double_t high_exp_tail_ratio3 = par[29];

  // Background: params 30-31
  Double_t bkg_const = par[30];
  Double_t bkg_slope = par[31];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t step1_par[3] = {mu1, sigma1, step_amplitude1};
  Double_t low_tail1_par[6] = {mu1,
                               sigma1,
                               low_exp_tail_amplitude1,
                               low_exp_tail_ratio1,
                               low_lin_tail_amplitude1,
                               low_lin_tail_slope1};
  Double_t high_tail1_par[4] = {mu1, sigma1, high_exp_tail_amplitude1,
                                high_exp_tail_ratio1};

  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t step2_par[3] = {mu2, sigma2, step_amplitude2};
  Double_t low_tail2_par[6] = {mu2,
                               sigma2,
                               low_exp_tail_amplitude2,
                               low_exp_tail_ratio2,
                               low_lin_tail_amplitude2,
                               low_lin_tail_slope2};
  Double_t high_tail2_par[4] = {mu2, sigma2, high_exp_tail_amplitude2,
                                high_exp_tail_ratio2};

  Double_t gaus3_par[3] = {mu3, sigma3, gaus_amplitude3};
  Double_t step3_par[3] = {mu3, sigma3, step_amplitude3};
  Double_t low_tail3_par[6] = {mu3,
                               sigma3,
                               low_exp_tail_amplitude3,
                               low_exp_tail_ratio3,
                               low_lin_tail_amplitude3,
                               low_lin_tail_slope3};
  Double_t high_tail3_par[4] = {mu3, sigma3, high_exp_tail_amplitude3,
                                high_exp_tail_ratio3};

  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Step(x, step1_par) +
         LowTail(x, low_tail1_par) + HighTail(x, high_tail1_par) +
         Gaussian(x, gaus2_par) + Step(x, step2_par) +
         LowTail(x, low_tail2_par) + HighTail(x, high_tail2_par) +
         Gaussian(x, gaus3_par) + Step(x, step3_par) +
         LowTail(x, low_tail3_par) + HighTail(x, high_tail3_par) +
         LinearBackground(x, bkg_par);
}

FittingUtils::FittingUtils(TH1 *working_hist, Float_t fit_range_low,
                           Float_t fit_range_high, Bool_t use_flat_background,
                           Bool_t use_step, Bool_t use_low_exp_tail,
                           Bool_t use_low_lin_tail, Bool_t use_high_exp_tail) {

  working_hist_ = static_cast<TH1 *>(working_hist->Clone());
  // Detach from gDirectory so its lifetime is bound to this object, not to
  // whatever TFile happens to be current at construction time.
  working_hist_->SetDirectory(nullptr);
  fit_range_low_ = fit_range_low;
  fit_range_high_ = fit_range_high;
  use_flat_background_ = use_flat_background;
  use_step_ = use_step;
  use_low_exp_tail_ = use_low_exp_tail;
  use_low_lin_tail_ = use_low_lin_tail;
  use_high_exp_tail_ = use_high_exp_tail;
  use_manual_init_ = kFALSE;
  interactive_ = kFALSE;

  fit_function_ = new TF1("PeakFunction", &FittingFunctions::PeakFunction,
                          fit_range_low_, fit_range_high_, 12);

  fit_function_->SetParName(0, "Mu");
  fit_function_->SetParName(1, "Sigma");
  fit_function_->SetParName(2, "GausAmplitude");
  fit_function_->SetParName(3, "StepAmplitude");
  fit_function_->SetParName(4, "LowExpTailAmplitude");
  fit_function_->SetParName(5, "LowExpTailRatio");
  fit_function_->SetParName(6, "LowLinTailAmplitude");
  fit_function_->SetParName(7, "LowLinTailSlope");
  fit_function_->SetParName(8, "HighExpTailAmplitude");
  fit_function_->SetParName(9, "HighExpTailRatio");
  fit_function_->SetParName(10, "BkgConstant");
  fit_function_->SetParName(11, "LinBkgSlope");

  Double_t mu_init = (fit_range_low_ + fit_range_high_) / 2;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;

  Double_t max_bin = working_hist_->GetMaximumBin();
  Double_t peak_height = working_hist_->GetBinContent(max_bin);
  Double_t bkg_estimate = EstimateBackground();

  // Gaussian parameters
  fit_function_->SetParLimits(0, fit_range_low_, fit_range_high_);
  fit_function_->SetParameter(0, mu_init);
  fit_function_->SetParLimits(1, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParameter(1, sigma_init);
  fit_function_->SetParLimits(2, 0, peak_height * 0.999);
  fit_function_->SetParameter(2, peak_height * 0.999);

  // Step ratio (fraction of gaussian amplitude)
  fit_function_->SetParLimits(3, 0, 0.5);
  if (use_step_)
    fit_function_->SetParameter(3, 0);
  else
    fit_function_->FixParameter(3, 0);

  // Low tail parameters (amplitudes as ratios of gaussian)
  fit_function_->SetParLimits(4, 0, 0.5);
  if (use_low_exp_tail_)
    fit_function_->SetParameter(4, 0.1);
  else
    fit_function_->FixParameter(4, 0);

  fit_function_->SetParLimits(5, 1.0, 100);
  if (use_low_exp_tail_)
    fit_function_->SetParameter(5, 1.5);
  else
    fit_function_->FixParameter(5, 1);

  fit_function_->SetParLimits(6, 0, 0.5);
  if (use_low_lin_tail_)
    fit_function_->SetParameter(6, 0.1);
  else
    fit_function_->FixParameter(6, 0);

  fit_function_->SetParLimits(7, -0.1, 0.1);
  if (use_low_lin_tail_)
    fit_function_->SetParameter(7, 0);
  else
    fit_function_->FixParameter(7, 0);

  // High tail parameters (amplitude as ratio of gaussian)
  fit_function_->SetParLimits(8, 0, 0.5);
  if (use_high_exp_tail_)
    fit_function_->SetParameter(8, 0.1);
  else
    fit_function_->FixParameter(8, 0);

  fit_function_->SetParLimits(9, 1.0, 100);
  if (use_high_exp_tail_)
    fit_function_->SetParameter(9, 1.5);
  else
    fit_function_->FixParameter(9, 1);

  // Background parameters
  fit_function_->SetParLimits(10, 0, peak_height * 0.999);
  fit_function_->SetParameter(10, bkg_estimate);

  if (!use_flat_background_)
    fit_function_->SetParLimits(11, -1000, 1000);
  else
    fit_function_->FixParameter(11, 0);

  std::cout << "Fit configuration:" << std::endl;
  std::cout << std::endl;
  if (use_flat_background_) {
    std::cout << "Background: FLAT" << std::endl;
  } else {
    std::cout << "Background: LINEAR" << std::endl;
  }
  std::cout << "Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "Low exponential tail: "
            << (use_low_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "Low linear tail: "
            << (use_low_lin_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "High exponential tail: "
            << (use_high_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
}

FittingUtils::~FittingUtils() {
  delete fit_function_;
  fit_function_ = nullptr;
  delete working_hist_;
  working_hist_ = nullptr;
}

void FittingUtils::SortPeaksByMu(Int_t num_peaks) {
  // Bubble sort peak blocks (10 params each) by ascending mu
  for (Int_t i = 0; i < num_peaks - 1; i++) {
    for (Int_t j = 0; j < num_peaks - i - 1; j++) {
      Double_t mu_j = fit_function_->GetParameter(j * 10);
      Double_t mu_next = fit_function_->GetParameter((j + 1) * 10);
      if (mu_j > mu_next) {
        std::cout << "Sorting peaks: swapping peak " << j + 1 << " (mu=" << mu_j
                  << ") and peak " << j + 2 << " (mu=" << mu_next << ")"
                  << std::endl;
        for (Int_t k = 0; k < 10; k++) {
          Int_t idx_a = j * 10 + k;
          Int_t idx_b = (j + 1) * 10 + k;

          Double_t tmp_val = fit_function_->GetParameter(idx_a);
          Double_t tmp_err = fit_function_->GetParError(idx_a);
          Double_t lo_a, hi_a, lo_b, hi_b;
          fit_function_->GetParLimits(idx_a, lo_a, hi_a);
          fit_function_->GetParLimits(idx_b, lo_b, hi_b);

          fit_function_->SetParameter(idx_a,
                                      fit_function_->GetParameter(idx_b));
          fit_function_->SetParError(idx_a, fit_function_->GetParError(idx_b));
          fit_function_->SetParLimits(idx_a, lo_b, hi_b);

          fit_function_->SetParameter(idx_b, tmp_val);
          fit_function_->SetParError(idx_b, tmp_err);
          fit_function_->SetParLimits(idx_b, lo_a, hi_a);
        }
      }
    }
  }
}

void FittingUtils::SetManualParameters(const std::vector<Double_t> &params) {
  if (params.size() != (size_t)fit_function_->GetNpar()) {
    std::cerr << "ERROR: Manual parameters size (" << params.size()
              << ") doesn't match number of fit parameters ("
              << fit_function_->GetNpar() << ")" << std::endl;
    return;
  }

  manual_params_ = params;
  use_manual_init_ = kTRUE;

  // Apply the parameters immediately
  for (size_t i = 0; i < params.size(); i++) {
    fit_function_->FixParameter(i, params[i]);
  }

  std::cout << "Manual parameters set:" << std::endl;
  for (size_t i = 0; i < params.size(); i++) {
    std::cout << "  Par[" << i << "] " << fit_function_->GetParName(i) << " = "
              << params[i] << std::endl;
  }
}

void FittingUtils::SetManualParameter(Int_t index, Double_t value) {
  if (index < 0 || index >= fit_function_->GetNpar()) {
    std::cerr << "ERROR: Parameter index " << index << " out of range [0, "
              << fit_function_->GetNpar() - 1 << "]" << std::endl;
    return;
  }

  if (!use_manual_init_) {
    manual_params_.resize(fit_function_->GetNpar(), 0.0);
    use_manual_init_ = kTRUE;
  }

  manual_params_[index] = value;
  fit_function_->FixParameter(index, value);

  std::cout << "Set Par[" << index << "] " << fit_function_->GetParName(index)
            << " = " << value << std::endl;
}

void FittingUtils::AppendPeakGraphs(std::vector<TGraph *> &components,
                                    Int_t param_offset, Style_t line_style,
                                    TF1 *background, Int_t npts,
                                    Double_t x_step) {
  Double_t ga = fit_function_->GetParameter(param_offset + 2);
  Width_t line_width = PlottingUtils::GetLineWidth();

  TF1 *peak =
      new TF1(PlottingUtils::GetRandomName(), FittingFunctions::Gaussian,
              fit_range_low_, fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(param_offset + 0));
  peak->SetParameter(1, fit_function_->GetParameter(param_offset + 1));
  peak->SetParameter(2, ga);
  TGraph *peak_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    peak_graph->SetPoint(i, x, peak->Eval(x) + background->Eval(x));
  }
  peak_graph->SetLineColor(kBlack);
  peak_graph->SetLineStyle(line_style);
  peak_graph->SetLineWidth(line_width);
  components.push_back(peak_graph);
  delete peak;

  TF1 *step = new TF1(PlottingUtils::GetRandomName(), FittingFunctions::Step,
                      fit_range_low_, fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(param_offset + 0));
  step->SetParameter(1, fit_function_->GetParameter(param_offset + 1));
  step->SetParameter(2, fit_function_->GetParameter(param_offset + 3) * ga);
  TGraph *step_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    step_graph->SetPoint(i, x, step->Eval(x) + background->Eval(x));
  }
  step_graph->SetLineColor(kGray);
  step_graph->SetLineStyle(line_style);
  step_graph->SetLineWidth(line_width);
  components.push_back(step_graph);
  delete step;

  TF1 *low_tail =
      new TF1(PlottingUtils::GetRandomName(), FittingFunctions::LowTail,
              fit_range_low_, fit_range_high_, 6);
  low_tail->SetParameter(0, fit_function_->GetParameter(param_offset + 0));
  low_tail->SetParameter(1, fit_function_->GetParameter(param_offset + 1));
  low_tail->SetParameter(2, fit_function_->GetParameter(param_offset + 4) * ga);
  low_tail->SetParameter(3, fit_function_->GetParameter(param_offset + 5));
  low_tail->SetParameter(4, fit_function_->GetParameter(param_offset + 6) * ga);
  low_tail->SetParameter(5, fit_function_->GetParameter(param_offset + 7));
  TGraph *low_tail_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    low_tail_graph->SetPoint(i, x, low_tail->Eval(x) + background->Eval(x));
  }
  low_tail_graph->SetLineColor(kRed);
  low_tail_graph->SetLineStyle(line_style);
  low_tail_graph->SetLineWidth(line_width);
  components.push_back(low_tail_graph);
  delete low_tail;

  TF1 *high_tail =
      new TF1(PlottingUtils::GetRandomName(), FittingFunctions::HighTail,
              fit_range_low_, fit_range_high_, 4);
  high_tail->SetParameter(0, fit_function_->GetParameter(param_offset + 0));
  high_tail->SetParameter(1, fit_function_->GetParameter(param_offset + 1));
  high_tail->SetParameter(2,
                          fit_function_->GetParameter(param_offset + 8) * ga);
  high_tail->SetParameter(3, fit_function_->GetParameter(param_offset + 9));
  TGraph *high_tail_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    high_tail_graph->SetPoint(i, x, high_tail->Eval(x) + background->Eval(x));
  }
  high_tail_graph->SetLineColor(kOrange);
  high_tail_graph->SetLineStyle(line_style);
  high_tail_graph->SetLineWidth(line_width);
  components.push_back(high_tail_graph);
  delete high_tail;
}

void FittingUtils::PlotFitSinglePeak(const TString input_name,
                                     const TString peak_name,
                                     const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();

  TGraph *total_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    total_graph->SetPoint(i, x, fit_function_->Eval(x));
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(10));
  background->SetParameter(1, fit_function_->GetParameter(11));

  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    background_graph->SetPoint(i, x, background->Eval(x));
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);

  Double_t ga = fit_function_->GetParameter(2);

  TF1 *peak = new TF1("gaussian", FittingFunctions::Gaussian, fit_range_low_,
                      fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  TGraph *peak_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    peak_graph->SetPoint(i, x, peak->Eval(x) + background->Eval(x));
  }
  peak_graph->SetLineColor(kBlack);
  peak_graph->SetLineWidth(line_width);
  components.push_back(peak_graph);
  delete peak;

  if (TMath::Abs(fit_function_->GetParameter(3)) > 1e-6) {
    TF1 *step = new TF1("step", FittingFunctions::Step, fit_range_low_,
                        fit_range_high_, 3);
    step->SetParameter(0, fit_function_->GetParameter(0));
    step->SetParameter(1, fit_function_->GetParameter(1));
    step->SetParameter(2, fit_function_->GetParameter(3) * ga);
    TGraph *step_graph = new TGraph(npts);
    for (Int_t i = 0; i < npts; i++) {
      Double_t x = fit_range_low_ + i * x_step;
      step_graph->SetPoint(i, x, step->Eval(x) + background->Eval(x));
    }
    step_graph->SetLineColor(kGray);
    step_graph->SetLineWidth(line_width);
    components.push_back(step_graph);
    delete step;
  }

  if (TMath::Abs(fit_function_->GetParameter(4)) > 1e-6 ||
      TMath::Abs(fit_function_->GetParameter(6)) > 1e-6) {
    TF1 *low_tail = new TF1("lowtail", FittingFunctions::LowTail,
                            fit_range_low_, fit_range_high_, 6);
    low_tail->SetParameter(0, fit_function_->GetParameter(0));
    low_tail->SetParameter(1, fit_function_->GetParameter(1));
    low_tail->SetParameter(2, fit_function_->GetParameter(4) * ga);
    low_tail->SetParameter(3, fit_function_->GetParameter(5));
    low_tail->SetParameter(4, fit_function_->GetParameter(6) * ga);
    low_tail->SetParameter(5, fit_function_->GetParameter(7));
    TGraph *low_tail_graph = new TGraph(npts);
    for (Int_t i = 0; i < npts; i++) {
      Double_t x = fit_range_low_ + i * x_step;
      low_tail_graph->SetPoint(i, x, low_tail->Eval(x) + background->Eval(x));
    }
    low_tail_graph->SetLineColor(kRed);
    low_tail_graph->SetLineWidth(line_width);
    components.push_back(low_tail_graph);
    delete low_tail;
  }

  if (TMath::Abs(fit_function_->GetParameter(8)) > 1e-6) {
    TF1 *high_tail = new TF1("hightail", FittingFunctions::HighTail,
                             fit_range_low_, fit_range_high_, 4);
    high_tail->SetParameter(0, fit_function_->GetParameter(0));
    high_tail->SetParameter(1, fit_function_->GetParameter(1));
    high_tail->SetParameter(2, fit_function_->GetParameter(8) * ga);
    high_tail->SetParameter(3, fit_function_->GetParameter(9));
    TGraph *high_tail_graph = new TGraph(npts);
    for (Int_t i = 0; i < npts; i++) {
      Double_t x = fit_range_low_ + i * x_step;
      high_tail_graph->SetPoint(i, x, high_tail->Eval(x) + background->Eval(x));
    }
    high_tail_graph->SetLineColor(kOrange);
    high_tail_graph->SetLineWidth(line_width);
    components.push_back(high_tail_graph);
    delete high_tail;
  }

  delete background;

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

void FittingUtils::PlotFitDoublePeak(const TString input_name,
                                     const TString peak_name,
                                     const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();

  TGraph *total_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    total_graph->SetPoint(i, x, fit_function_->Eval(x));
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(20));
  background->SetParameter(1, fit_function_->GetParameter(21));

  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    background_graph->SetPoint(i, x, background->Eval(x));
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);

  AppendPeakGraphs(components, 0, 1, background, npts, x_step);
  AppendPeakGraphs(components, 10, 3, background, npts, x_step);

  delete background;

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

void FittingUtils::PlotFitTriplePeak(const TString input_name,
                                     const TString peak_name,
                                     const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();

  TGraph *total_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    total_graph->SetPoint(i, x, fit_function_->Eval(x));
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(30));
  background->SetParameter(1, fit_function_->GetParameter(31));

  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t x = fit_range_low_ + i * x_step;
    background_graph->SetPoint(i, x, background->Eval(x));
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);

  AppendPeakGraphs(components, 0, 1, background, npts, x_step);
  AppendPeakGraphs(components, 10, 3, background, npts, x_step);
  AppendPeakGraphs(components, 20, 4, background, npts, x_step);

  delete background;

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

Double_t FittingUtils::EstimateBackground() {
  Int_t left_bin = working_hist_->FindBin(fit_range_low_);
  Int_t right_bin = working_hist_->FindBin(fit_range_high_);

  Int_t n_sideband = (right_bin - left_bin) / 10;
  Double_t left_avg = 0, right_avg = 0;

  for (Int_t i = 0; i < n_sideband; i++) {
    left_avg += working_hist_->GetBinContent(left_bin + i);
    right_avg += working_hist_->GetBinContent(right_bin - i);
  }

  return (left_avg + right_avg) / (2.0 * n_sideband);
}

Double_t FittingUtils::ClampToBounds(Int_t param_index, Double_t value) {
  Double_t low, high;
  fit_function_->GetParLimits(param_index, low, high);

  if (low < high) {
    return TMath::Max(low, TMath::Min(value, high));
  }
  return value;
}

FitResult FittingUtils::FitSinglePeak(const TString input_name,
                                      const TString peak_name) {
  FitResult results;
  results.peaks.emplace_back(); // 1 peak, default -1

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      TFitResultPtr refit = working_hist_->Fit(fit_function_, "LSMRBENR+");
      if (refit.Get() && refit->IsValid())
        final_chi2 = refit->Chi2() / refit->Ndf();
      else if (fit_function_->GetNDF() > 0)
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(working_hist_, fit_function_,
                                     fit_range_low_, fit_range_high_, 1,
                                     peak_name + " / " + input_name)) {
        Double_t rlo_tmp, rhi_tmp;
        fit_function_->GetRange(rlo_tmp, rhi_tmp);
        fit_range_low_ = rlo_tmp;
        fit_range_high_ = rhi_tmp;
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    // Automated fitting: initial fit, component testing, final fit

    // Fix all optional components to disabled state for baseline fit
    fit_function_->FixParameter(3, 0); // StepAmplitude
    fit_function_->FixParameter(4, 0); // LowExpTailAmplitude
    fit_function_->FixParameter(5, 1); // LowExpTailRatio
    fit_function_->FixParameter(6, 0); // LowLinTailAmplitude
    fit_function_->FixParameter(7, 0); // LowLinTailSlope
    fit_function_->FixParameter(8, 0); // HighExpTailAmplitude
    fit_function_->FixParameter(9, 1); // HighExpTailRatio

    if (use_flat_background_) {
      fit_function_->FixParameter(11, 0);
    }

    if (use_manual_init_) {
      std::cout << "Using manually initialized parameters" << std::endl;
      for (size_t i = 0; i < manual_params_.size(); i++) {
        fit_function_->SetParameter(i, manual_params_[i]);
      }
      std::cout << "Skipping auto-initialization, using provided values"
                << std::endl;
    } else {
      if (!use_flat_background_) {
        TF1 *bkg_only = new TF1("bkg_temp", FittingFunctions::LinearBackground,
                                fit_range_low_, fit_range_high_, 2);
        Double_t exclude_low =
            fit_function_->GetParameter(0) - 3 * fit_function_->GetParameter(1);
        working_hist_->Fit(bkg_only, "QN0R", "", fit_range_low_, exclude_low);
        fit_function_->SetParameter(
            10, ClampToBounds(10, bkg_only->GetParameter(0)));
        fit_function_->SetParameter(
            11, ClampToBounds(11, bkg_only->GetParameter(1)));
        delete bkg_only;
      }
    }

    // Initial fit with just Gaussian + Background
    TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (!initial_fit.Get() || !initial_fit->IsValid()) {
      std::cout << "ERROR: Initial fit failed" << std::endl;
      return results;
    }

    Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

    Double_t gaus_amp = TMath::Abs(fit_function_->GetParameter(2));
    Double_t peak_height =
        working_hist_->GetBinContent(working_hist_->GetMaximumBin());
    Double_t range_width = fit_range_high_ - fit_range_low_;
    Double_t bkg_estimate = EstimateBackground();

    std::cout << "Background mode: "
              << (use_flat_background_ ? "FLAT" : "LINEAR") << std::endl;

    Int_t npar = fit_function_->GetNpar();
    std::vector<Double_t> best_params(npar);
    std::vector<Double_t> best_errors(npar);
    for (Int_t i = 0; i < npar; i++) {
      best_params[i] = fit_function_->GetParameter(i);
      best_errors[i] = fit_function_->GetParError(i);
    }

    //  Low-side group testing (step + low tails together)
    Bool_t any_low_side = use_step_ || use_low_exp_tail_ || use_low_lin_tail_;

    if (any_low_side) {
      std::cout << "Testing low-side component group..." << std::endl;

      if (use_step_) {
        fit_function_->ReleaseParameter(3);
        fit_function_->SetParLimits(3, 0, 0.5);
        if (!use_manual_init_)
          fit_function_->SetParameter(3, 0.15);
      }
      if (use_low_exp_tail_) {
        fit_function_->ReleaseParameter(4);
        fit_function_->ReleaseParameter(5);
        fit_function_->SetParLimits(4, 0, 0.5);
        fit_function_->SetParLimits(5, 1.0, 100);
        if (!use_manual_init_) {
          fit_function_->SetParameter(4, 0.15);
          fit_function_->SetParameter(5, 1.5);
        }
      }
      if (use_low_lin_tail_) {
        fit_function_->ReleaseParameter(6);
        fit_function_->ReleaseParameter(7);
        fit_function_->SetParLimits(6, 0, 0.5);
        fit_function_->SetParLimits(7, -0.1, 0.1);
        if (!use_manual_init_) {
          fit_function_->SetParameter(6, 0.15);
          fit_function_->SetParameter(7, 0);
        }
      }

      TFitResultPtr group_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (group_fit.Get() && group_fit->IsValid()) {
        Double_t chi2_group = group_fit->Chi2() / group_fit->Ndf();
        std::cout << "Low-side group chi2/ndf: " << chi2_group << " vs "
                  << best_chi2 << std::endl;

        if (chi2_group < best_chi2) {
          std::cout << "Low-side group ACCEPTED, pruning..." << std::endl;
          best_chi2 = chi2_group;
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }

          // Prune step
          if (use_step_) {
            fit_function_->FixParameter(3, 0);
            TFitResultPtr prune_fit =
                working_hist_->Fit(fit_function_, "LSMBNQ0R");
            if (prune_fit.Get() && prune_fit->IsValid() &&
                prune_fit->Chi2() / prune_fit->Ndf() <= best_chi2) {
              std::cout << "  Step pruned (not needed)" << std::endl;
              best_chi2 = prune_fit->Chi2() / prune_fit->Ndf();
              for (Int_t i = 0; i < npar; i++) {
                best_params[i] = fit_function_->GetParameter(i);
                best_errors[i] = fit_function_->GetParError(i);
              }
            } else {
              std::cout << "  Step retained" << std::endl;
              fit_function_->ReleaseParameter(3);
              fit_function_->SetParLimits(3, 0, 0.5);
              for (Int_t i = 0; i < npar; i++) {
                fit_function_->SetParameter(i, best_params[i]);
                fit_function_->SetParError(i, best_errors[i]);
              }
            }
          }

          // Prune low exp tail
          if (use_low_exp_tail_) {
            fit_function_->FixParameter(4, 0);
            fit_function_->FixParameter(5, 1);
            TFitResultPtr prune_fit =
                working_hist_->Fit(fit_function_, "LSMBNQ0R");
            if (prune_fit.Get() && prune_fit->IsValid() &&
                prune_fit->Chi2() / prune_fit->Ndf() <= best_chi2) {
              std::cout << "  Low exp tail pruned (not needed)" << std::endl;
              best_chi2 = prune_fit->Chi2() / prune_fit->Ndf();
              for (Int_t i = 0; i < npar; i++) {
                best_params[i] = fit_function_->GetParameter(i);
                best_errors[i] = fit_function_->GetParError(i);
              }
            } else {
              std::cout << "  Low exp tail retained" << std::endl;
              fit_function_->ReleaseParameter(4);
              fit_function_->ReleaseParameter(5);
              fit_function_->SetParLimits(4, 0, 0.5);
              fit_function_->SetParLimits(5, 1.0, 100);
              for (Int_t i = 0; i < npar; i++) {
                fit_function_->SetParameter(i, best_params[i]);
                fit_function_->SetParError(i, best_errors[i]);
              }
            }
          }

          // Prune low lin tail
          if (use_low_lin_tail_) {
            fit_function_->FixParameter(6, 0);
            fit_function_->FixParameter(7, 0);
            TFitResultPtr prune_fit =
                working_hist_->Fit(fit_function_, "LSMBNQ0R");
            if (prune_fit.Get() && prune_fit->IsValid() &&
                prune_fit->Chi2() / prune_fit->Ndf() <= best_chi2) {
              std::cout << "  Low lin tail pruned (not needed)" << std::endl;
              best_chi2 = prune_fit->Chi2() / prune_fit->Ndf();
              for (Int_t i = 0; i < npar; i++) {
                best_params[i] = fit_function_->GetParameter(i);
                best_errors[i] = fit_function_->GetParError(i);
              }
            } else {
              std::cout << "  Low lin tail retained" << std::endl;
              fit_function_->ReleaseParameter(6);
              fit_function_->ReleaseParameter(7);
              fit_function_->SetParLimits(6, 0, 0.5);
              fit_function_->SetParLimits(7, -0.1, 0.1);
              for (Int_t i = 0; i < npar; i++) {
                fit_function_->SetParameter(i, best_params[i]);
                fit_function_->SetParError(i, best_errors[i]);
              }
            }
          }
        } else {
          std::cout << "Low-side group REJECTED" << std::endl;
          if (use_step_)
            fit_function_->FixParameter(3, 0);
          if (use_low_exp_tail_) {
            fit_function_->FixParameter(4, 0);
            fit_function_->FixParameter(5, 1);
          }
          if (use_low_lin_tail_) {
            fit_function_->FixParameter(6, 0);
            fit_function_->FixParameter(7, 0);
          }
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "Low-side group fit FAILED" << std::endl;
        if (use_step_)
          fit_function_->FixParameter(3, 0);
        if (use_low_exp_tail_) {
          fit_function_->FixParameter(4, 0);
          fit_function_->FixParameter(5, 1);
        }
        if (use_low_lin_tail_) {
          fit_function_->FixParameter(6, 0);
          fit_function_->FixParameter(7, 0);
        }
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  High tail testing (independent)
    if (use_high_exp_tail_) {
      std::cout << "Testing high exponential tail..." << std::endl;

      fit_function_->ReleaseParameter(8);
      fit_function_->ReleaseParameter(9);
      fit_function_->SetParLimits(8, 0, 0.5);
      fit_function_->SetParLimits(9, 1.0, 100);

      if (!use_manual_init_) {
        fit_function_->SetParameter(8, 0.15);
        fit_function_->SetParameter(9, 1.5);
      }

      TFitResultPtr htail_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (htail_fit.Get() && htail_fit->IsValid()) {
        Double_t chi2_with_htail = htail_fit->Chi2() / htail_fit->Ndf();
        std::cout << "Chi2/ndf: " << chi2_with_htail << " vs " << best_chi2
                  << std::endl;

        if (chi2_with_htail < best_chi2) {
          std::cout << "High exp tail ACCEPTED" << std::endl;
          best_chi2 = chi2_with_htail;
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          std::cout << "High exp tail REJECTED" << std::endl;
          fit_function_->FixParameter(8, 0);
          fit_function_->FixParameter(9, 1);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "High exp tail fit FAILED" << std::endl;
        fit_function_->FixParameter(8, 0);
        fit_function_->FixParameter(9, 1);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  Final fit
    std::cout << "Final fit with selected components..." << std::endl;
    for (Int_t i = 0; i < npar; i++) {
      fit_function_->SetParameter(i, best_params[i]);
      fit_function_->SetParError(i, best_errors[i]);
    }

    if (use_flat_background_) {
      fit_function_->FixParameter(11, 0);
    }

    TFitResultPtr final_fit = working_hist_->Fit(fit_function_, "LSMRBENR+");

    if (final_fit.Get() && final_fit->IsValid()) {
      final_chi2 = final_fit->Chi2() / final_fit->Ndf();
      fit_valid = kTRUE;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    }
  }

  if (fit_valid) {
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitSinglePeak(input_name, peak_name, chi2label);

    PeakFitResult peak;
    peak.mu = fit_function_->GetParameter(0);
    peak.mu_error = fit_function_->GetParError(0);
    peak.sigma = fit_function_->GetParameter(1);
    peak.sigma_error = fit_function_->GetParError(1);
    peak.gaus_amplitude = fit_function_->GetParameter(2);
    peak.gaus_amplitude_error = fit_function_->GetParError(2);
    // Convert ratios back to absolute amplitudes
    Double_t ga = peak.gaus_amplitude;
    peak.step_amplitude = fit_function_->GetParameter(3) * ga;
    peak.step_amplitude_error = fit_function_->GetParError(3) * ga;
    peak.low_exp_tail_amplitude = fit_function_->GetParameter(4) * ga;
    peak.low_exp_tail_amplitude_error = fit_function_->GetParError(4) * ga;
    peak.low_exp_tail_ratio = fit_function_->GetParameter(5);
    peak.low_exp_tail_ratio_error = fit_function_->GetParError(5);
    peak.low_lin_tail_amplitude = fit_function_->GetParameter(6) * ga;
    peak.low_lin_tail_amplitude_error = fit_function_->GetParError(6) * ga;
    peak.low_lin_tail_slope = fit_function_->GetParameter(7);
    peak.low_lin_tail_slope_error = fit_function_->GetParError(7);
    peak.high_exp_tail_amplitude = fit_function_->GetParameter(8) * ga;
    peak.high_exp_tail_amplitude_error = fit_function_->GetParError(8) * ga;
    peak.high_exp_tail_ratio = fit_function_->GetParameter(9);
    peak.high_exp_tail_ratio_error = fit_function_->GetParError(9);

    results.peaks[0] = peak;
    results.bkg_constant = fit_function_->GetParameter(10);
    results.bkg_constant_error = fit_function_->GetParError(10);
    results.lin_bkg_slope = fit_function_->GetParameter(11);
    results.lin_bkg_slope_error = fit_function_->GetParError(11);
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Fit did not converge" << std::endl;
  }

  return results;
}

FitResult FittingUtils::FitDoublePeak(const TString input_name,
                                      const TString peak_name,
                                      Double_t mu1_init, Double_t mu2_init) {
  FitResult results;
  results.peaks.emplace_back(); // peak 1, default -1
  results.peaks.emplace_back(); // peak 2, default -1

  if (mu1_init > mu2_init) {
    std::cout << "WARNING: mu1_init > mu2_init, swapping initial values"
              << std::endl;
    Double_t temp = mu1_init;
    mu1_init = mu2_init;
    mu2_init = temp;
  }

  delete fit_function_;
  fit_function_ = new TF1("DoublePeak", &FittingFunctions::DoublePeakFunction,
                          fit_range_low_, fit_range_high_, 22);

  // Peak 1 param names (offset 0)
  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowExpTailAmplitude1");
  fit_function_->SetParName(5, "LowExpTailRatio1");
  fit_function_->SetParName(6, "LowLinTailAmplitude1");
  fit_function_->SetParName(7, "LowLinTailSlope1");
  fit_function_->SetParName(8, "HighExpTailAmplitude1");
  fit_function_->SetParName(9, "HighExpTailRatio1");

  // Peak 2 param names (offset 10)
  fit_function_->SetParName(10, "Mu2");
  fit_function_->SetParName(11, "Sigma2");
  fit_function_->SetParName(12, "GausAmplitude2");
  fit_function_->SetParName(13, "StepAmplitude2");
  fit_function_->SetParName(14, "LowExpTailAmplitude2");
  fit_function_->SetParName(15, "LowExpTailRatio2");
  fit_function_->SetParName(16, "LowLinTailAmplitude2");
  fit_function_->SetParName(17, "LowLinTailSlope2");
  fit_function_->SetParName(18, "HighExpTailAmplitude2");
  fit_function_->SetParName(19, "HighExpTailRatio2");

  // Background param names
  fit_function_->SetParName(20, "BkgConst");
  fit_function_->SetParName(21, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  // Peak 1 and 2 limits (amplitude params are ratios)
  for (Int_t p = 0; p < 2; p++) {
    Int_t o = p * 10;
    fit_function_->SetParLimits(o + 0, fit_range_low_, fit_range_high_);
    fit_function_->SetParLimits(o + 1, range_width * 0.001, range_width * 0.5);
    fit_function_->SetParLimits(o + 2, 0, peak_height * 0.999);
    fit_function_->SetParLimits(o + 3, 0, 0.5);
    fit_function_->SetParLimits(o + 4, 0, 0.5);
    fit_function_->SetParLimits(o + 5, 1.0, 100);
    fit_function_->SetParLimits(o + 6, 0, 0.5);
    fit_function_->SetParLimits(o + 7, -0.1, 0.1);
    fit_function_->SetParLimits(o + 8, 0, 0.5);
    fit_function_->SetParLimits(o + 9, 1.0, 100);
  }

  fit_function_->SetParameter(0, mu1_init);
  fit_function_->SetParameter(1, sigma_init);
  fit_function_->SetParameter(2, peak_height * 0.999);

  fit_function_->SetParameter(10, mu2_init);
  fit_function_->SetParameter(11, sigma_init);
  fit_function_->SetParameter(12, peak_height * 0.999);

  // Fix disabled optional components on both peaks
  for (Int_t p = 0; p < 2; p++) {
    Int_t o = p * 10;
    if (use_step_)
      fit_function_->SetParameter(o + 3, 0);
    else
      fit_function_->FixParameter(o + 3, 0);

    if (use_low_exp_tail_) {
      fit_function_->SetParameter(o + 4, 0);
      fit_function_->SetParameter(o + 5, 1.5);
    } else {
      fit_function_->FixParameter(o + 4, 0);
      fit_function_->FixParameter(o + 5, 1);
    }

    if (use_low_lin_tail_) {
      fit_function_->SetParameter(o + 6, 0);
      fit_function_->SetParameter(o + 7, 0);
    } else {
      fit_function_->FixParameter(o + 6, 0);
      fit_function_->FixParameter(o + 7, 0);
    }

    if (use_high_exp_tail_) {
      fit_function_->SetParameter(o + 8, 0);
      fit_function_->SetParameter(o + 9, 1.5);
    } else {
      fit_function_->FixParameter(o + 8, 0);
      fit_function_->FixParameter(o + 9, 1);
    }
  }

  // Background
  fit_function_->SetParLimits(20, 0, peak_height * 0.999);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(21, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  }
  fit_function_->SetParameter(20, bkg_estimate);
  fit_function_->SetParameter(21, 0);
  if (use_flat_background_) {
    fit_function_->FixParameter(21, 0);
  }

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      TFitResultPtr refit = working_hist_->Fit(fit_function_, "LSMRBENR+");
      if (refit.Get() && refit->IsValid())
        final_chi2 = refit->Chi2() / refit->Ndf();
      else if (fit_function_->GetNDF() > 0)
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(working_hist_, fit_function_,
                                     fit_range_low_, fit_range_high_, 2,
                                     peak_name + " / " + input_name)) {
        Double_t rlo_tmp, rhi_tmp;
        fit_function_->GetRange(rlo_tmp, rhi_tmp);
        fit_range_low_ = rlo_tmp;
        fit_range_high_ = rhi_tmp;
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    // Fix all optional components for baseline fit
    fit_function_->FixParameter(3, 0);
    fit_function_->FixParameter(4, 0);
    fit_function_->FixParameter(5, 1);
    fit_function_->FixParameter(6, 0);
    fit_function_->FixParameter(7, 0);
    fit_function_->FixParameter(8, 0);
    fit_function_->FixParameter(9, 1);
    fit_function_->FixParameter(13, 0);
    fit_function_->FixParameter(14, 0);
    fit_function_->FixParameter(15, 1);
    fit_function_->FixParameter(16, 0);
    fit_function_->FixParameter(17, 0);
    fit_function_->FixParameter(18, 0);
    fit_function_->FixParameter(19, 1);

    // Initial fit
    TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (!initial_fit.Get() || !initial_fit->IsValid()) {
      std::cout << "ERROR: Initial double peak fit failed" << std::endl;
      return results;
    }

    Double_t gaus_amp1 = TMath::Abs(fit_function_->GetParameter(2));
    Double_t gaus_amp2 = TMath::Abs(fit_function_->GetParameter(12));

    Int_t npar = fit_function_->GetNpar();
    std::vector<Double_t> best_params(npar);
    std::vector<Double_t> best_errors(npar);
    for (Int_t i = 0; i < npar; i++) {
      best_params[i] = fit_function_->GetParameter(i);
      best_errors[i] = fit_function_->GetParError(i);
    }

    Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

    //  Low-side group for peak 1 (offset 0)
    {
      std::cout << "Testing low-side group for peak1..." << std::endl;

      fit_function_->ReleaseParameter(3);
      fit_function_->SetParLimits(3, 0, 0.5);
      fit_function_->SetParameter(3, 0.15);

      fit_function_->ReleaseParameter(4);
      fit_function_->ReleaseParameter(5);
      fit_function_->SetParLimits(4, 0, 0.5);
      fit_function_->SetParLimits(5, 1.0, 100);
      fit_function_->SetParameter(4, 0.15);
      fit_function_->SetParameter(5, 1.5);

      fit_function_->ReleaseParameter(6);
      fit_function_->ReleaseParameter(7);
      fit_function_->SetParLimits(6, 0, 0.5);
      fit_function_->SetParLimits(7, -0.1, 0.1);
      fit_function_->SetParameter(6, 0.15);
      fit_function_->SetParameter(7, 0);

      TFitResultPtr group_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (group_fit.Get() && group_fit->IsValid() &&
          group_fit->Chi2() / group_fit->Ndf() < best_chi2) {
        std::cout << "Low-side group peak1 ACCEPTED, pruning..." << std::endl;
        best_chi2 = group_fit->Chi2() / group_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }

        // Prune step1
        fit_function_->FixParameter(3, 0);
        TFitResultPtr p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  Step1 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(3);
          fit_function_->SetParLimits(3, 0, 0.5);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low exp tail1
        fit_function_->FixParameter(4, 0);
        fit_function_->FixParameter(5, 1);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowExpTail1 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(4);
          fit_function_->ReleaseParameter(5);
          fit_function_->SetParLimits(4, 0, 0.5);
          fit_function_->SetParLimits(5, 1.0, 100);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low lin tail1
        fit_function_->FixParameter(6, 0);
        fit_function_->FixParameter(7, 0);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowLinTail1 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(6);
          fit_function_->ReleaseParameter(7);
          fit_function_->SetParLimits(6, 0, 0.5);
          fit_function_->SetParLimits(7, -0.1, 0.1);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "Low-side group peak1 REJECTED" << std::endl;
        fit_function_->FixParameter(3, 0);
        fit_function_->FixParameter(4, 0);
        fit_function_->FixParameter(5, 1);
        fit_function_->FixParameter(6, 0);
        fit_function_->FixParameter(7, 0);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  High tail for peak 2 (outer component, no inter-peak overlap)
    {
      std::cout << "Testing high tail for peak2..." << std::endl;
      fit_function_->ReleaseParameter(18);
      fit_function_->ReleaseParameter(19);
      fit_function_->SetParLimits(18, 0, 0.5);
      fit_function_->SetParLimits(19, 1.0, 100);
      fit_function_->SetParameter(18, 0.15);
      fit_function_->SetParameter(19, 1.5);

      TFitResultPtr ht_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");
      if (ht_fit.Get() && ht_fit->IsValid() &&
          ht_fit->Chi2() / ht_fit->Ndf() < best_chi2) {
        std::cout << "HighTail2 ACCEPTED" << std::endl;
        best_chi2 = ht_fit->Chi2() / ht_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "HighTail2 REJECTED" << std::endl;
        fit_function_->FixParameter(18, 0);
        fit_function_->FixParameter(19, 1);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  Inter-peak group: peak1 high tail + peak2 low-side (both affect the
    //  region between the two peaks, so they must be tested jointly)
    {
      std::cout
          << "Testing inter-peak group (peak1 high tail + peak2 low-side)..."
          << std::endl;

      // Release peak1 high tail
      fit_function_->ReleaseParameter(8);
      fit_function_->ReleaseParameter(9);
      fit_function_->SetParLimits(8, 0, 0.5);
      fit_function_->SetParLimits(9, 1.0, 100);
      fit_function_->SetParameter(8, 0.15);
      fit_function_->SetParameter(9, 1.5);

      // Release peak2 low-side group
      fit_function_->ReleaseParameter(13);
      fit_function_->SetParLimits(13, 0, 0.5);
      fit_function_->SetParameter(13, 0.15);

      fit_function_->ReleaseParameter(14);
      fit_function_->ReleaseParameter(15);
      fit_function_->SetParLimits(14, 0, 0.5);
      fit_function_->SetParLimits(15, 1.0, 100);
      fit_function_->SetParameter(14, 0.15);
      fit_function_->SetParameter(15, 1.5);

      fit_function_->ReleaseParameter(16);
      fit_function_->ReleaseParameter(17);
      fit_function_->SetParLimits(16, 0, 0.5);
      fit_function_->SetParLimits(17, -0.1, 0.1);
      fit_function_->SetParameter(16, 0.15);
      fit_function_->SetParameter(17, 0);

      TFitResultPtr group_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (group_fit.Get() && group_fit->IsValid() &&
          group_fit->Chi2() / group_fit->Ndf() < best_chi2) {
        std::cout << "Inter-peak group ACCEPTED, pruning..." << std::endl;
        best_chi2 = group_fit->Chi2() / group_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }

        // Prune peak1 high tail
        fit_function_->FixParameter(8, 0);
        fit_function_->FixParameter(9, 1);
        TFitResultPtr p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  HighTail1 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(8);
          fit_function_->ReleaseParameter(9);
          fit_function_->SetParLimits(8, 0, 0.5);
          fit_function_->SetParLimits(9, 1.0, 100);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune step2
        fit_function_->FixParameter(13, 0);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  Step2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(13);
          fit_function_->SetParLimits(13, 0, 0.5);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low exp tail2
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowExpTail2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(14);
          fit_function_->ReleaseParameter(15);
          fit_function_->SetParLimits(14, 0, 0.5);
          fit_function_->SetParLimits(15, 1.0, 100);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low lin tail2
        fit_function_->FixParameter(16, 0);
        fit_function_->FixParameter(17, 0);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowLinTail2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(16);
          fit_function_->ReleaseParameter(17);
          fit_function_->SetParLimits(16, 0, 0.5);
          fit_function_->SetParLimits(17, -0.1, 0.1);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "Inter-peak group REJECTED" << std::endl;
        fit_function_->FixParameter(8, 0);
        fit_function_->FixParameter(9, 1);
        fit_function_->FixParameter(13, 0);
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        fit_function_->FixParameter(16, 0);
        fit_function_->FixParameter(17, 0);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  Final fit
    std::cout << "Final fit with selected components..." << std::endl;
    for (Int_t i = 0; i < npar; i++) {
      fit_function_->SetParameter(i, best_params[i]);
      fit_function_->SetParError(i, best_errors[i]);
    }
    if (use_flat_background_) {
      fit_function_->FixParameter(21, 0);
    }

    TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

    if (fit_result.Get() && fit_result->IsValid()) {
      final_chi2 = fit_result->Chi2() / fit_result->Ndf();
      std::cout << "Double peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
      fit_valid = kTRUE;
    } else {
      std::cout << "ERROR: Double peak fit failed to converge" << std::endl;
    }
  }

  if (fit_valid) {
    SortPeaksByMu(2);

    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    for (Int_t pk = 0; pk < 2; pk++) {
      Int_t o = pk * 10;
      PeakFitResult p;
      p.mu = fit_function_->GetParameter(o + 0);
      p.mu_error = fit_function_->GetParError(o + 0);
      p.sigma = fit_function_->GetParameter(o + 1);
      p.sigma_error = fit_function_->GetParError(o + 1);
      p.gaus_amplitude = fit_function_->GetParameter(o + 2);
      p.gaus_amplitude_error = fit_function_->GetParError(o + 2);
      // Convert ratios back to absolute amplitudes
      Double_t ga = p.gaus_amplitude;
      p.step_amplitude = fit_function_->GetParameter(o + 3) * ga;
      p.step_amplitude_error = fit_function_->GetParError(o + 3) * ga;
      p.low_exp_tail_amplitude = fit_function_->GetParameter(o + 4) * ga;
      p.low_exp_tail_amplitude_error = fit_function_->GetParError(o + 4) * ga;
      p.low_exp_tail_ratio = fit_function_->GetParameter(o + 5);
      p.low_exp_tail_ratio_error = fit_function_->GetParError(o + 5);
      p.low_lin_tail_amplitude = fit_function_->GetParameter(o + 6) * ga;
      p.low_lin_tail_amplitude_error = fit_function_->GetParError(o + 6) * ga;
      p.low_lin_tail_slope = fit_function_->GetParameter(o + 7);
      p.low_lin_tail_slope_error = fit_function_->GetParError(o + 7);
      p.high_exp_tail_amplitude = fit_function_->GetParameter(o + 8) * ga;
      p.high_exp_tail_amplitude_error = fit_function_->GetParError(o + 8) * ga;
      p.high_exp_tail_ratio = fit_function_->GetParameter(o + 9);
      p.high_exp_tail_ratio_error = fit_function_->GetParError(o + 9);
      results.peaks[pk] = p;
    }

    results.bkg_constant = fit_function_->GetParameter(20);
    results.bkg_constant_error = fit_function_->GetParError(20);
    results.lin_bkg_slope = fit_function_->GetParameter(21);
    results.lin_bkg_slope_error = fit_function_->GetParError(21);
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Double peak fit failed" << std::endl;
  }

  return results;
}

FitResult FittingUtils::FitDoublePeak(const TString input_name,
                                      const TString peak_name,
                                      const PeakFitResult &constrained_peak,
                                      Double_t mu2_init) {
  FitResult results;
  results.peaks.emplace_back(); // peak 1, default -1
  results.peaks.emplace_back(); // peak 2, default -1

  delete fit_function_;
  fit_function_ = new TF1("DoublePeak", &FittingFunctions::DoublePeakFunction,
                          fit_range_low_, fit_range_high_, 22);

  // Peak 1 param names (offset 0)
  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowExpTailAmplitude1");
  fit_function_->SetParName(5, "LowExpTailRatio1");
  fit_function_->SetParName(6, "LowLinTailAmplitude1");
  fit_function_->SetParName(7, "LowLinTailSlope1");
  fit_function_->SetParName(8, "HighExpTailAmplitude1");
  fit_function_->SetParName(9, "HighExpTailRatio1");

  // Peak 2 param names (offset 10)
  fit_function_->SetParName(10, "Mu2");
  fit_function_->SetParName(11, "Sigma2");
  fit_function_->SetParName(12, "GausAmplitude2");
  fit_function_->SetParName(13, "StepAmplitude2");
  fit_function_->SetParName(14, "LowExpTailAmplitude2");
  fit_function_->SetParName(15, "LowExpTailRatio2");
  fit_function_->SetParName(16, "LowLinTailAmplitude2");
  fit_function_->SetParName(17, "LowLinTailSlope2");
  fit_function_->SetParName(18, "HighExpTailAmplitude2");
  fit_function_->SetParName(19, "HighExpTailRatio2");

  fit_function_->SetParName(20, "BkgConst");
  fit_function_->SetParName(21, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  // Constrained peak 1: fix all shape params (mu, sigma, amplitude ratios,
  // decay, slope). Only gaussian amplitude is free — tail/step amplitudes
  // scale automatically via their fixed ratios.
  {
    const PeakFitResult &cp = constrained_peak;

    fit_function_->FixParameter(0, cp.mu);
    fit_function_->FixParameter(1, cp.sigma);

    fit_function_->SetParLimits(2, 0, peak_height * 2);
    fit_function_->SetParameter(2, cp.gaus_amplitude);

    fit_function_->FixParameter(3, cp.step_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(4,
                                cp.low_exp_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(
        5, cp.low_exp_tail_ratio > 0 ? cp.low_exp_tail_ratio : 1);
    fit_function_->FixParameter(6,
                                cp.low_lin_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(7, cp.low_lin_tail_slope);

    fit_function_->FixParameter(8,
                                cp.high_exp_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(
        9, cp.high_exp_tail_ratio > 0 ? cp.high_exp_tail_ratio : 1);
  }

  // Free peak 2 (all optional components fixed to 0)
  fit_function_->SetParLimits(10, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(11, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(12, 0, peak_height * 0.999);
  fit_function_->SetParLimits(13, 0, 0.5);
  fit_function_->SetParLimits(14, 0, 0.5);
  fit_function_->SetParLimits(15, 1.0, 100);
  fit_function_->SetParLimits(16, 0, 0.5);
  fit_function_->SetParLimits(17, -0.1, 0.1);
  fit_function_->SetParLimits(18, 0, 0.5);
  fit_function_->SetParLimits(19, 1.0, 100);

  fit_function_->SetParameter(10, mu2_init);
  fit_function_->SetParameter(11, sigma_init);
  fit_function_->SetParameter(12, peak_height * 0.999);

  // Fix disabled optional components on free peak 2
  if (use_step_)
    fit_function_->SetParameter(13, 0);
  else
    fit_function_->FixParameter(13, 0);

  if (use_low_exp_tail_) {
    fit_function_->SetParameter(14, 0);
    fit_function_->SetParameter(15, 1.5);
  } else {
    fit_function_->FixParameter(14, 0);
    fit_function_->FixParameter(15, 1);
  }

  if (use_low_lin_tail_) {
    fit_function_->SetParameter(16, 0);
    fit_function_->SetParameter(17, 0);
  } else {
    fit_function_->FixParameter(16, 0);
    fit_function_->FixParameter(17, 0);
  }

  if (use_high_exp_tail_) {
    fit_function_->SetParameter(18, 0);
    fit_function_->SetParameter(19, 1.5);
  } else {
    fit_function_->FixParameter(18, 0);
    fit_function_->FixParameter(19, 1);
  }

  // Background
  fit_function_->SetParLimits(20, 0, peak_height * 0.999);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(21, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  }
  fit_function_->SetParameter(20, bkg_estimate);
  fit_function_->SetParameter(21, 0);
  if (use_flat_background_) {
    fit_function_->FixParameter(21, 0);
  }

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      TFitResultPtr refit = working_hist_->Fit(fit_function_, "LSMRBENR+");
      if (refit.Get() && refit->IsValid())
        final_chi2 = refit->Chi2() / refit->Ndf();
      else if (fit_function_->GetNDF() > 0)
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(working_hist_, fit_function_,
                                     fit_range_low_, fit_range_high_, 2,
                                     peak_name + " / " + input_name)) {
        Double_t rlo_tmp, rhi_tmp;
        fit_function_->GetRange(rlo_tmp, rhi_tmp);
        fit_range_low_ = rlo_tmp;
        fit_range_high_ = rhi_tmp;
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    // Fix all optional components on peak 2 for baseline fit
    fit_function_->FixParameter(13, 0);
    fit_function_->FixParameter(14, 0);
    fit_function_->FixParameter(15, 1);
    fit_function_->FixParameter(16, 0);
    fit_function_->FixParameter(17, 0);
    fit_function_->FixParameter(18, 0);
    fit_function_->FixParameter(19, 1);

    // Initial fit
    TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (!initial_fit.Get() || !initial_fit->IsValid()) {
      std::cout << "ERROR: Initial double peak fit (constrained) failed"
                << std::endl;
      return results;
    }

    Double_t gaus_amp2 = TMath::Abs(fit_function_->GetParameter(12));

    Int_t npar = fit_function_->GetNpar();
    std::vector<Double_t> best_params(npar);
    std::vector<Double_t> best_errors(npar);
    for (Int_t i = 0; i < npar; i++) {
      best_params[i] = fit_function_->GetParameter(i);
      best_errors[i] = fit_function_->GetParError(i);
    }

    Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

    //  Low-side group for peak 2 only (offset 10)
    {
      std::cout << "Testing low-side group for peak2..." << std::endl;

      fit_function_->ReleaseParameter(13);
      fit_function_->SetParLimits(13, 0, 0.5);
      fit_function_->SetParameter(13, 0.15);

      fit_function_->ReleaseParameter(14);
      fit_function_->ReleaseParameter(15);
      fit_function_->SetParLimits(14, 0, 0.5);
      fit_function_->SetParLimits(15, 1.0, 100);
      fit_function_->SetParameter(14, 0.15);
      fit_function_->SetParameter(15, 1.5);

      fit_function_->ReleaseParameter(16);
      fit_function_->ReleaseParameter(17);
      fit_function_->SetParLimits(16, 0, 0.5);
      fit_function_->SetParLimits(17, -0.1, 0.1);
      fit_function_->SetParameter(16, 0.15);
      fit_function_->SetParameter(17, 0);

      TFitResultPtr group_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (group_fit.Get() && group_fit->IsValid() &&
          group_fit->Chi2() / group_fit->Ndf() < best_chi2) {
        std::cout << "Low-side group peak2 ACCEPTED, pruning..." << std::endl;
        best_chi2 = group_fit->Chi2() / group_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }

        // Prune step2
        fit_function_->FixParameter(13, 0);
        TFitResultPtr p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  Step2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(13);
          fit_function_->SetParLimits(13, 0, 0.5);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low exp tail2
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowExpTail2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(14);
          fit_function_->ReleaseParameter(15);
          fit_function_->SetParLimits(14, 0, 0.5);
          fit_function_->SetParLimits(15, 1.0, 100);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low lin tail2
        fit_function_->FixParameter(16, 0);
        fit_function_->FixParameter(17, 0);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowLinTail2 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(16);
          fit_function_->ReleaseParameter(17);
          fit_function_->SetParLimits(16, 0, 0.5);
          fit_function_->SetParLimits(17, -0.1, 0.1);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "Low-side group peak2 REJECTED" << std::endl;
        fit_function_->FixParameter(13, 0);
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        fit_function_->FixParameter(16, 0);
        fit_function_->FixParameter(17, 0);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  High tail for peak 2
    {
      std::cout << "Testing high tail for peak2..." << std::endl;
      fit_function_->ReleaseParameter(18);
      fit_function_->ReleaseParameter(19);
      fit_function_->SetParLimits(18, 0, 0.5);
      fit_function_->SetParLimits(19, 1.0, 100);
      fit_function_->SetParameter(18, 0.15);
      fit_function_->SetParameter(19, 1.5);

      TFitResultPtr ht_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");
      if (ht_fit.Get() && ht_fit->IsValid() &&
          ht_fit->Chi2() / ht_fit->Ndf() < best_chi2) {
        std::cout << "HighTail2 ACCEPTED" << std::endl;
        best_chi2 = ht_fit->Chi2() / ht_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "HighTail2 REJECTED" << std::endl;
        fit_function_->FixParameter(18, 0);
        fit_function_->FixParameter(19, 1);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  Final fit
    std::cout << "Final fit with selected components..." << std::endl;
    for (Int_t i = 0; i < npar; i++) {
      fit_function_->SetParameter(i, best_params[i]);
      fit_function_->SetParError(i, best_errors[i]);
    }
    if (use_flat_background_) {
      fit_function_->FixParameter(21, 0);
    }

    TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

    if (fit_result.Get() && fit_result->IsValid()) {
      final_chi2 = fit_result->Chi2() / fit_result->Ndf();
      std::cout << "Double peak fit (constrained) converged successfully"
                << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
      fit_valid = kTRUE;
    } else {
      std::cout << "ERROR: Double peak fit (constrained) failed to converge"
                << std::endl;
    }
  }

  if (fit_valid) {
    SortPeaksByMu(2);

    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    for (Int_t pk = 0; pk < 2; pk++) {
      Int_t o = pk * 10;
      PeakFitResult p;
      p.mu = fit_function_->GetParameter(o + 0);
      p.mu_error = fit_function_->GetParError(o + 0);
      p.sigma = fit_function_->GetParameter(o + 1);
      p.sigma_error = fit_function_->GetParError(o + 1);
      p.gaus_amplitude = fit_function_->GetParameter(o + 2);
      p.gaus_amplitude_error = fit_function_->GetParError(o + 2);
      // Convert ratios back to absolute amplitudes
      Double_t ga = p.gaus_amplitude;
      p.step_amplitude = fit_function_->GetParameter(o + 3) * ga;
      p.step_amplitude_error = fit_function_->GetParError(o + 3) * ga;
      p.low_exp_tail_amplitude = fit_function_->GetParameter(o + 4) * ga;
      p.low_exp_tail_amplitude_error = fit_function_->GetParError(o + 4) * ga;
      p.low_exp_tail_ratio = fit_function_->GetParameter(o + 5);
      p.low_exp_tail_ratio_error = fit_function_->GetParError(o + 5);
      p.low_lin_tail_amplitude = fit_function_->GetParameter(o + 6) * ga;
      p.low_lin_tail_amplitude_error = fit_function_->GetParError(o + 6) * ga;
      p.low_lin_tail_slope = fit_function_->GetParameter(o + 7);
      p.low_lin_tail_slope_error = fit_function_->GetParError(o + 7);
      p.high_exp_tail_amplitude = fit_function_->GetParameter(o + 8) * ga;
      p.high_exp_tail_amplitude_error = fit_function_->GetParError(o + 8) * ga;
      p.high_exp_tail_ratio = fit_function_->GetParameter(o + 9);
      p.high_exp_tail_ratio_error = fit_function_->GetParError(o + 9);
      results.peaks[pk] = p;
    }

    results.bkg_constant = fit_function_->GetParameter(20);
    results.bkg_constant_error = fit_function_->GetParError(20);
    results.lin_bkg_slope = fit_function_->GetParameter(21);
    results.lin_bkg_slope_error = fit_function_->GetParError(21);
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Double peak fit (constrained) failed" << std::endl;
  }

  return results;
}

FitResult FittingUtils::FitTriplePeak(const TString input_name,
                                      const TString peak_name,
                                      const FitResult &constrained_peaks,
                                      Double_t mu3_init) {
  FitResult results;
  results.peaks.emplace_back(); // peak 1, default -1
  results.peaks.emplace_back(); // peak 2, default -1
  results.peaks.emplace_back(); // peak 3, default -1

  delete fit_function_;
  fit_function_ = new TF1("TriplePeak", &FittingFunctions::TriplePeakFunction,
                          fit_range_low_, fit_range_high_, 32);

  // Peak 1 param names (offset 0)
  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowExpTailAmplitude1");
  fit_function_->SetParName(5, "LowExpTailRatio1");
  fit_function_->SetParName(6, "LowLinTailAmplitude1");
  fit_function_->SetParName(7, "LowLinTailSlope1");
  fit_function_->SetParName(8, "HighExpTailAmplitude1");
  fit_function_->SetParName(9, "HighExpTailRatio1");

  // Peak 2 param names (offset 10)
  fit_function_->SetParName(10, "Mu2");
  fit_function_->SetParName(11, "Sigma2");
  fit_function_->SetParName(12, "GausAmplitude2");
  fit_function_->SetParName(13, "StepAmplitude2");
  fit_function_->SetParName(14, "LowExpTailAmplitude2");
  fit_function_->SetParName(15, "LowExpTailRatio2");
  fit_function_->SetParName(16, "LowLinTailAmplitude2");
  fit_function_->SetParName(17, "LowLinTailSlope2");
  fit_function_->SetParName(18, "HighExpTailAmplitude2");
  fit_function_->SetParName(19, "HighExpTailRatio2");

  // Peak 3 param names (offset 20)
  fit_function_->SetParName(20, "Mu3");
  fit_function_->SetParName(21, "Sigma3");
  fit_function_->SetParName(22, "GausAmplitude3");
  fit_function_->SetParName(23, "StepAmplitude3");
  fit_function_->SetParName(24, "LowExpTailAmplitude3");
  fit_function_->SetParName(25, "LowExpTailRatio3");
  fit_function_->SetParName(26, "LowLinTailAmplitude3");
  fit_function_->SetParName(27, "LowLinTailSlope3");
  fit_function_->SetParName(28, "HighExpTailAmplitude3");
  fit_function_->SetParName(29, "HighExpTailRatio3");

  fit_function_->SetParName(30, "BkgConst");
  fit_function_->SetParName(31, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  // Constrained peaks: fix shape (mu, sigma, decay, slope), allow amplitudes
  // to float with [0, ...] bounds. Disable components that were off in the
  // background fit.
  const PeakFitResult *cpeaks[2] = {&constrained_peaks.peaks[0],
                                    &constrained_peaks.peaks[1]};
  for (Int_t pk = 0; pk < 2; pk++) {
    Int_t o = pk * 10;
    const PeakFitResult &cp = *cpeaks[pk];

    // Mu and sigma: fixed (shape)
    fit_function_->FixParameter(o + 0, cp.mu);
    fit_function_->FixParameter(o + 1, cp.sigma);

    // Gaussian amplitude: free, positive
    fit_function_->SetParLimits(o + 2, 0, peak_height * 2);
    fit_function_->SetParameter(o + 2, cp.gaus_amplitude);

    // Fix all amplitude ratios and shape params — they scale with
    // gaus_amplitude
    fit_function_->FixParameter(o + 3, cp.step_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(o + 4,
                                cp.low_exp_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(
        o + 5, cp.low_exp_tail_ratio > 0 ? cp.low_exp_tail_ratio : 1);
    fit_function_->FixParameter(o + 6,
                                cp.low_lin_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(o + 7, cp.low_lin_tail_slope);
    fit_function_->FixParameter(o + 8,
                                cp.high_exp_tail_amplitude / cp.gaus_amplitude);
    fit_function_->FixParameter(
        o + 9, cp.high_exp_tail_ratio > 0 ? cp.high_exp_tail_ratio : 1);
  }

  // Free peak 3 (offset 20)
  fit_function_->SetParLimits(20, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(21, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(22, 0, peak_height * 0.999);
  fit_function_->SetParLimits(23, 0, 0.5);
  fit_function_->SetParLimits(24, 0, 0.5);
  fit_function_->SetParLimits(25, 1.0, 100);
  fit_function_->SetParLimits(26, 0, 0.5);
  fit_function_->SetParLimits(27, -0.1, 0.1);
  fit_function_->SetParLimits(28, 0, 0.5);
  fit_function_->SetParLimits(29, 1.0, 100);

  fit_function_->SetParameter(20, mu3_init);
  fit_function_->SetParameter(21, sigma_init);
  fit_function_->SetParameter(22, peak_height * 0.999);

  // Fix disabled optional components on free peak 3
  if (use_step_)
    fit_function_->SetParameter(23, 0);
  else
    fit_function_->FixParameter(23, 0);

  if (use_low_exp_tail_) {
    fit_function_->SetParameter(24, 0);
    fit_function_->SetParameter(25, 1.5);
  } else {
    fit_function_->FixParameter(24, 0);
    fit_function_->FixParameter(25, 1);
  }

  if (use_low_lin_tail_) {
    fit_function_->SetParameter(26, 0);
    fit_function_->SetParameter(27, 0);
  } else {
    fit_function_->FixParameter(26, 0);
    fit_function_->FixParameter(27, 0);
  }

  if (use_high_exp_tail_) {
    fit_function_->SetParameter(28, 0);
    fit_function_->SetParameter(29, 1.5);
  } else {
    fit_function_->FixParameter(28, 0);
    fit_function_->FixParameter(29, 1);
  }

  // Background
  fit_function_->SetParLimits(30, 0, peak_height * 0.999);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(31, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  }
  fit_function_->SetParameter(30, bkg_estimate);
  fit_function_->SetParameter(31, 0);
  if (use_flat_background_) {
    fit_function_->FixParameter(31, 0);
  }

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      TFitResultPtr refit = working_hist_->Fit(fit_function_, "LSMRBENR+");
      if (refit.Get() && refit->IsValid())
        final_chi2 = refit->Chi2() / refit->Ndf();
      else if (fit_function_->GetNDF() > 0)
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveFitEditor(working_hist_, fit_function_,
                                     fit_range_low_, fit_range_high_, 3,
                                     peak_name + " / " + input_name)) {
        Double_t rlo_tmp, rhi_tmp;
        fit_function_->GetRange(rlo_tmp, rhi_tmp);
        fit_range_low_ = rlo_tmp;
        fit_range_high_ = rhi_tmp;
        final_chi2 = fit_function_->GetChisquare() / fit_function_->GetNDF();
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    // Fix all optional components on peak 3 for baseline fit
    fit_function_->FixParameter(23, 0);
    fit_function_->FixParameter(24, 0);
    fit_function_->FixParameter(25, 1);
    fit_function_->FixParameter(26, 0);
    fit_function_->FixParameter(27, 0);
    fit_function_->FixParameter(28, 0);
    fit_function_->FixParameter(29, 1);

    // Initial fit
    TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (!initial_fit.Get() || !initial_fit->IsValid()) {
      std::cout << "ERROR: Initial triple peak fit failed" << std::endl;
      return results;
    }

    Double_t gaus_amp3 = TMath::Abs(fit_function_->GetParameter(22));

    Int_t npar = fit_function_->GetNpar();
    std::vector<Double_t> best_params(npar);
    std::vector<Double_t> best_errors(npar);
    for (Int_t i = 0; i < npar; i++) {
      best_params[i] = fit_function_->GetParameter(i);
      best_errors[i] = fit_function_->GetParError(i);
    }

    Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

    //  Low-side group for peak 3 (offset 20)
    {
      std::cout << "Testing low-side group for peak3..." << std::endl;

      fit_function_->ReleaseParameter(23);
      fit_function_->SetParLimits(23, 0, 0.5);
      fit_function_->SetParameter(23, 0.15);

      fit_function_->ReleaseParameter(24);
      fit_function_->ReleaseParameter(25);
      fit_function_->SetParLimits(24, 0, 0.5);
      fit_function_->SetParLimits(25, 1.0, 100);
      fit_function_->SetParameter(24, 0.15);
      fit_function_->SetParameter(25, 1.5);

      fit_function_->ReleaseParameter(26);
      fit_function_->ReleaseParameter(27);
      fit_function_->SetParLimits(26, 0, 0.5);
      fit_function_->SetParLimits(27, -0.1, 0.1);
      fit_function_->SetParameter(26, 0.15);
      fit_function_->SetParameter(27, 0);

      TFitResultPtr group_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

      if (group_fit.Get() && group_fit->IsValid() &&
          group_fit->Chi2() / group_fit->Ndf() < best_chi2) {
        std::cout << "Low-side group peak3 ACCEPTED, pruning..." << std::endl;
        best_chi2 = group_fit->Chi2() / group_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }

        // Prune step3
        fit_function_->FixParameter(23, 0);
        TFitResultPtr p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  Step3 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(23);
          fit_function_->SetParLimits(23, 0, 0.5);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low exp tail3
        fit_function_->FixParameter(24, 0);
        fit_function_->FixParameter(25, 1);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowExpTail3 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(24);
          fit_function_->ReleaseParameter(25);
          fit_function_->SetParLimits(24, 0, 0.5);
          fit_function_->SetParLimits(25, 1.0, 100);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }

        // Prune low lin tail3
        fit_function_->FixParameter(26, 0);
        fit_function_->FixParameter(27, 0);
        p = working_hist_->Fit(fit_function_, "LSMBNQ0R");
        if (p.Get() && p->IsValid() && p->Chi2() / p->Ndf() <= best_chi2) {
          std::cout << "  LowLinTail3 pruned" << std::endl;
          best_chi2 = p->Chi2() / p->Ndf();
          for (Int_t i = 0; i < npar; i++) {
            best_params[i] = fit_function_->GetParameter(i);
            best_errors[i] = fit_function_->GetParError(i);
          }
        } else {
          fit_function_->ReleaseParameter(26);
          fit_function_->ReleaseParameter(27);
          fit_function_->SetParLimits(26, 0, 0.5);
          fit_function_->SetParLimits(27, -0.1, 0.1);
          for (Int_t i = 0; i < npar; i++) {
            fit_function_->SetParameter(i, best_params[i]);
            fit_function_->SetParError(i, best_errors[i]);
          }
        }
      } else {
        std::cout << "Low-side group peak3 REJECTED" << std::endl;
        fit_function_->FixParameter(23, 0);
        fit_function_->FixParameter(24, 0);
        fit_function_->FixParameter(25, 1);
        fit_function_->FixParameter(26, 0);
        fit_function_->FixParameter(27, 0);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  High tail for peak 3
    {
      std::cout << "Testing high tail for peak3..." << std::endl;
      fit_function_->ReleaseParameter(28);
      fit_function_->ReleaseParameter(29);
      fit_function_->SetParLimits(28, 0, 0.5);
      fit_function_->SetParLimits(29, 1.0, 100);
      fit_function_->SetParameter(28, 0.15);
      fit_function_->SetParameter(29, 1.5);

      TFitResultPtr ht_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");
      if (ht_fit.Get() && ht_fit->IsValid() &&
          ht_fit->Chi2() / ht_fit->Ndf() < best_chi2) {
        std::cout << "HighTail3 ACCEPTED" << std::endl;
        best_chi2 = ht_fit->Chi2() / ht_fit->Ndf();
        for (Int_t i = 0; i < npar; i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "HighTail3 REJECTED" << std::endl;
        fit_function_->FixParameter(28, 0);
        fit_function_->FixParameter(29, 1);
        for (Int_t i = 0; i < npar; i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    }

    //  Final fit
    std::cout << "Final fit with selected components..." << std::endl;
    for (Int_t i = 0; i < npar; i++) {
      fit_function_->SetParameter(i, best_params[i]);
      fit_function_->SetParError(i, best_errors[i]);
    }
    if (use_flat_background_) {
      fit_function_->FixParameter(31, 0);
    }

    TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

    if (fit_result.Get() && fit_result->IsValid()) {
      final_chi2 = fit_result->Chi2() / fit_result->Ndf();
      std::cout << "Triple peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
      fit_valid = kTRUE;
    } else {
      std::cout << "ERROR: Triple peak fit failed to converge" << std::endl;
    }
  }

  if (fit_valid) {
    SortPeaksByMu(3);

    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitTriplePeak(input_name, peak_name, chi2label);

    for (Int_t pk = 0; pk < 3; pk++) {
      Int_t o = pk * 10;
      PeakFitResult p;
      p.mu = fit_function_->GetParameter(o + 0);
      p.mu_error = fit_function_->GetParError(o + 0);
      p.sigma = fit_function_->GetParameter(o + 1);
      p.sigma_error = fit_function_->GetParError(o + 1);
      p.gaus_amplitude = fit_function_->GetParameter(o + 2);
      p.gaus_amplitude_error = fit_function_->GetParError(o + 2);
      // Convert ratios back to absolute amplitudes
      Double_t ga = p.gaus_amplitude;
      p.step_amplitude = fit_function_->GetParameter(o + 3) * ga;
      p.step_amplitude_error = fit_function_->GetParError(o + 3) * ga;
      p.low_exp_tail_amplitude = fit_function_->GetParameter(o + 4) * ga;
      p.low_exp_tail_amplitude_error = fit_function_->GetParError(o + 4) * ga;
      p.low_exp_tail_ratio = fit_function_->GetParameter(o + 5);
      p.low_exp_tail_ratio_error = fit_function_->GetParError(o + 5);
      p.low_lin_tail_amplitude = fit_function_->GetParameter(o + 6) * ga;
      p.low_lin_tail_amplitude_error = fit_function_->GetParError(o + 6) * ga;
      p.low_lin_tail_slope = fit_function_->GetParameter(o + 7);
      p.low_lin_tail_slope_error = fit_function_->GetParError(o + 7);
      p.high_exp_tail_amplitude = fit_function_->GetParameter(o + 8) * ga;
      p.high_exp_tail_amplitude_error = fit_function_->GetParError(o + 8) * ga;
      p.high_exp_tail_ratio = fit_function_->GetParameter(o + 9);
      p.high_exp_tail_ratio_error = fit_function_->GetParError(o + 9);
      results.peaks[pk] = p;
    }

    results.bkg_constant = fit_function_->GetParameter(30);
    results.bkg_constant_error = fit_function_->GetParError(30);
    results.lin_bkg_slope = fit_function_->GetParameter(31);
    results.lin_bkg_slope_error = fit_function_->GetParError(31);
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Triple peak fit failed" << std::endl;
  }

  return results;
}
