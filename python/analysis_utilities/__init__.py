"""Analysis utilities for nuclear measurement data."""

from analysis_utilities.init_utils import set_root_preferences
from analysis_utilities.io import (get_root_files_base_dir, open_for_reading,
                                   open_for_writing, set_root_files_base_dir)

__version__ = "@VERSION@"

__all__ = [
    "load_cpp_library",
    "set_root_preferences",
    "set_root_files_base_dir",
    "get_root_files_base_dir",
    "open_for_reading",
    "open_for_writing",
]

_cpp_loaded = False


def load_cpp_library():
    """Load the C++ Analysis-Utilities library into ROOT.

    After calling this, ROOT.PlottingUtils, ROOT.InitUtils,
    ROOT.PlotSaveFormat, and ROOT.PlotSaveOptions are available.

    Returns:
        ROOT module (for convenience)
    """
    global _cpp_loaded
    import ROOT

    if not _cpp_loaded:
        if ROOT.gSystem.Load("libanalysis-utils") < 0:
            raise RuntimeError(
                "Could not load libanalysis-utils.so. "
                "Make sure LD_LIBRARY_PATH includes the library directory.")
        ROOT.gInterpreter.Declare('#include "PlottingUtils.hpp"')
        ROOT.gInterpreter.Declare('#include "InitUtils.hpp"')
        ROOT.gInterpreter.Declare('#include "IOUtils.hpp"')
        _cpp_loaded = True

    return ROOT
