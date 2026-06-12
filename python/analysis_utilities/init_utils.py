"""Python port of :cpp:class:`InitUtils` (the parts applicable to Python).

Binary file converters stay on the C++ side; only the ROOT environment
setup (styling, batch mode, standard output directories) is mirrored here.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Any, Optional, Union

PathLike = Union[str, os.PathLike]


def set_root_preferences(
    save_format: Optional[Any] = None,
    plots_dir: Optional[PathLike] = None,
    root_files_dir: Optional[PathLike] = None,
) -> Any:
    """Configure ROOT environment and plotting defaults.

    Mirrors :cpp:func:`InitUtils::SetROOTPreferences`: loads the C++
    library, calls :cpp:func:`PlottingUtils::SetStylePreferences`, forces
    the ROOT style, enables batch mode, and ensures the plots and root_files
    directories exist.

    Parameters
    ----------
    save_format : ROOT.PlotSaveFormat, optional
        Output format for plots. Defaults to
        :cpp:enumerator:`PlotSaveFormat::kPNG` if not provided.
    plots_dir : str | os.PathLike, optional
        Base directory for plot output. Resolved to an absolute path and
        forwarded to :cpp:func:`PlottingUtils::SetPlotsBaseDir` so that
        :cpp:func:`PlottingUtils::SaveFigure` writes
        ``<plots_dir>/<subdir>/<name>.<ext>``. If omitted, a warning is
        emitted and the CWD-relative default ``"plots"`` is used.
    root_files_dir : str | os.PathLike, optional
        Directory for ROOT files. Resolved to an absolute path. If omitted,
        a warning is emitted and the CWD-relative default ``"root_files"``
        is used. The C++ IO layer (binary converters) still writes to
        ``"root_files"`` relative to CWD; this controls only the directory
        that is created here.

    Returns
    -------
    ROOT
        The ROOT module, for convenience.
    """
    from analysis_utilities import load_cpp_library

    ROOT = load_cpp_library()
    if save_format is None:
        save_format = ROOT.PlotSaveFormat.kPNG
    ROOT.PlottingUtils.SetStylePreferences(save_format)
    ROOT.gROOT.ForceStyle(True)
    ROOT.gROOT.SetBatch(True)

    if plots_dir is None:
        warnings.warn(
            "set_root_preferences called without plots_dir; defaulting to "
            "CWD-relative 'plots'. Pass an absolute path to drive output "
            "into a project root.",
            stacklevel=2,
        )
        resolved_plots_dir = Path("plots").resolve()
    else:
        resolved_plots_dir = Path(plots_dir).expanduser().resolve()
    ROOT.PlottingUtils.SetPlotsBaseDir(str(resolved_plots_dir))

    if root_files_dir is None:
        warnings.warn(
            "set_root_preferences called without root_files_dir; "
            "defaulting to CWD-relative 'root_files'.",
            stacklevel=2,
        )
        resolved_root_files_dir = Path("root_files").resolve()
    else:
        from analysis_utilities.io import (get_root_files_base_dir,
                                           set_root_files_base_dir)

        set_root_files_base_dir(root_files_dir)
        resolved_root_files_dir = Path(get_root_files_base_dir())

    resolved_plots_dir.mkdir(parents=True, exist_ok=True)
    resolved_root_files_dir.mkdir(parents=True, exist_ok=True)

    return ROOT
