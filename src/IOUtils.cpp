#include "IOUtils.hpp"
#include <TDirectory.h>
#include <TROOT.h>
#include <mutex>

namespace {
TString g_root_files_base = "root_files";
Bool_t g_thread_safe = kFALSE;
std::recursive_mutex g_file_open_mutex;

TString JoinPath(const TString &subpath) {
  if (gSystem->IsAbsoluteFileName(subpath))
    return subpath;
  char *tmp = gSystem->ConcatFileName(g_root_files_base.Data(), subpath.Data());
  TString full(tmp);
  delete[] tmp;
  return full;
}

void ApplyCompressionDefaults(TFile *file) {
  // ZSTD compression for every TTree/TFile produced by this package; set on the
  // file before branches are created so the setting propagates to the baskets.
  if (file && !file->IsZombie()) {
    file->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZSTD);
    file->SetCompressionLevel(5);
  }
}
} // namespace

void IO::SetRootFilesBaseDir(const TString &dir) {
  TString d = dir;
  while (d.Length() > 0 && d[d.Length() - 1] == '/')
    d.Chop();
  g_root_files_base = d;
}

TString IO::GetRootFilesBaseDir() { return g_root_files_base; }

void IO::SetThreadSafe(Bool_t enabled) {
  g_thread_safe = enabled;
  if (enabled) {
    ROOT::EnableThreadSafety();
  }
}

Bool_t IO::IsThreadSafe() { return g_thread_safe; }

IO::ScopedRootLock::ScopedRootLock() : engaged_(g_thread_safe) {
  if (engaged_)
    g_file_open_mutex.lock();
}

IO::ScopedRootLock::~ScopedRootLock() {
  if (engaged_)
    g_file_open_mutex.unlock();
}

TFile *IO::OpenForReading(const TString &subpath) {
  TString full = JoinPath(subpath);
  if (g_thread_safe) {
    std::lock_guard<std::recursive_mutex> lock(g_file_open_mutex);
    TDirectory::TContext ctxGuard;
    return new TFile(full, "READ");
  }
  return new TFile(full, "READ");
}

TFile *IO::OpenForWriting(const TString &subpath, const TString mode) {
  TString full = JoinPath(subpath);
  TString parent = gSystem->DirName(full);
  gSystem->mkdir(parent, kTRUE);
  if (g_thread_safe) {
    std::lock_guard<std::recursive_mutex> lock(g_file_open_mutex);
    TFile *file = new TFile(full, mode);
    ApplyCompressionDefaults(file);
    return file;
  }
  TFile *file = new TFile(full, mode);
  ApplyCompressionDefaults(file);
  return file;
}
