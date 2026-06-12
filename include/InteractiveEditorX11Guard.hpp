#ifndef INTERACTIVEEDITORX11GUARD_H
#define INTERACTIVEEDITORX11GUARD_H

#include <Rtypes.h>
#include <dlfcn.h>

// Guard against ROOT's default X11 error handler turning a recoverable X
// protocol error into a SIGSEGV while an interactive editor pumps its own
// event loop.
//
// When a TGTextEntry (e.g. the number entries) is redrawn through the TTF
// backend, TGX11TTF::DrawString issues an XGetGeometry round-trip on the
// widget's drawable. If that drawable is momentarily stale, the X server
// answers with BadWindow/BadDrawable. Xlib then calls ROOT's
// RootX11ErrorHandler, which looks the offending window up by id and calls a
// virtual method on it to format the message; if that window object is gone
// the handler itself crashes. The recoverable protocol error thus becomes a
// fatal SIGSEGV that takes down the whole ROOT session.
//
// While the editor loop runs we install a no-op handler that swallows such
// errors (Xlib ignores the return value for non-fatal errors and carries on),
// then restore the previous handler. XSetErrorHandler is resolved at runtime
// via dlsym so this library carries no direct link-time dependency on libX11;
// the symbol is provided by libX11, loaded transitively by ROOT's GUI
// libraries whenever a gClient exists.

typedef int (*AUXErrorHandler)(void *display, void *event);
typedef AUXErrorHandler (*AUXSetErrorHandlerFn)(AUXErrorHandler);

inline int AUIgnoreXError(void *, void *) { return 0; }

struct AUXErrorHandlerSave {
  AUXSetErrorHandlerFn set_fn;
  AUXErrorHandler prev;
};

// Install the tolerant handler. Returns what is needed to restore the prior
// state; if XSetErrorHandler cannot be resolved the call is a no-op and
// restoring is harmless.
inline AUXErrorHandlerSave AUInstallTolerantXErrorHandler() {
  AUXErrorHandlerSave s;
  s.set_fn = (AUXSetErrorHandlerFn)dlsym(RTLD_DEFAULT, "XSetErrorHandler");
  if (!s.set_fn) {
    // libX11 may have been brought in with local scope; promote and retry.
    void *h = dlopen("libX11.so.6", RTLD_NOW | RTLD_GLOBAL | RTLD_NOLOAD);
    if (h) {
      s.set_fn = (AUXSetErrorHandlerFn)dlsym(h, "XSetErrorHandler");
    }
  }
  s.prev = s.set_fn ? s.set_fn(AUIgnoreXError) : (AUXErrorHandler)0;
  return s;
}

inline void AURestoreXErrorHandler(AUXErrorHandlerSave s) {
  if (s.set_fn) {
    s.set_fn(s.prev);
  }
}

#endif
