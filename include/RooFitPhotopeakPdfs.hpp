#ifndef ROOFITPHOTOPEAKPDFS_H
#define ROOFITPHOTOPEAKPDFS_H

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooFit/EvalContext.h>
#include <RooRealProxy.h>

class RooStepShelf : public RooAbsPdf {
public:
  RooStepShelf() {}
  RooStepShelf(const char *name, const char *title, RooAbsReal &x,
               RooAbsReal &mu, RooAbsReal &sigma);
  RooStepShelf(const RooStepShelf &other, const char *name = nullptr);
  TObject *clone(const char *newname = nullptr) const override {
    return new RooStepShelf(*this, newname);
  }

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  void doEval(RooFit::EvalContext &ctx) const override;
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  Double_t evaluate() const override;

  ClassDefOverride(RooStepShelf, 1)
};

class RooLowExpTail : public RooAbsPdf {
public:
  RooLowExpTail() {}
  RooLowExpTail(const char *name, const char *title, RooAbsReal &x,
                RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &tau);
  RooLowExpTail(const RooLowExpTail &other, const char *name = nullptr);
  TObject *clone(const char *newname = nullptr) const override {
    return new RooLowExpTail(*this, newname);
  }

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  void doEval(RooFit::EvalContext &ctx) const override;
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy tau_;
  Double_t evaluate() const override;

  ClassDefOverride(RooLowExpTail, 1)
};

class RooLowLinTail : public RooAbsPdf {
public:
  RooLowLinTail() {}
  RooLowLinTail(const char *name, const char *title, RooAbsReal &x,
                RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &slope);
  RooLowLinTail(const RooLowLinTail &other, const char *name = nullptr);
  TObject *clone(const char *newname = nullptr) const override {
    return new RooLowLinTail(*this, newname);
  }

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  void doEval(RooFit::EvalContext &ctx) const override;
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy slope_;
  Double_t evaluate() const override;

  ClassDefOverride(RooLowLinTail, 1)
};

class RooHighExpTail : public RooAbsPdf {
public:
  RooHighExpTail() {}
  RooHighExpTail(const char *name, const char *title, RooAbsReal &x,
                 RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &tau);
  RooHighExpTail(const RooHighExpTail &other, const char *name = nullptr);
  TObject *clone(const char *newname = nullptr) const override {
    return new RooHighExpTail(*this, newname);
  }

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  void doEval(RooFit::EvalContext &ctx) const override;
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy tau_;
  Double_t evaluate() const override;

  ClassDefOverride(RooHighExpTail, 1)
};

#endif
