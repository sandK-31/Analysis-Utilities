#ifndef IOUTILS_H
#define IOUTILS_H

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

namespace IO {
void SetRootFilesBaseDir(const TString &dir);
TString GetRootFilesBaseDir();
void SetThreadSafe(Bool_t enabled = kTRUE);
Bool_t IsThreadSafe();
TFile *OpenForReading(const TString &subpath);
TFile *OpenForWriting(const TString &subpath, const TString mode = "RECREATE");

class ScopedRootLock {
public:
  ScopedRootLock();
  ~ScopedRootLock();

private:
  Bool_t engaged_;
  ScopedRootLock(const ScopedRootLock &);
  ScopedRootLock &operator=(const ScopedRootLock &);
};
} // namespace IO

#endif
