"""Tests for analysis_utilities.init_utils."""

from __future__ import annotations

import warnings

import pytest

ROOT = pytest.importorskip("ROOT")

try:
    from analysis_utilities import load_cpp_library, set_root_preferences

    load_cpp_library()
except (RuntimeError, OSError) as exc:
    pytest.skip(f"C++ library not loadable: {exc}", allow_module_level=True)


def test_set_root_preferences_routes_savefigure_into_plots_dir(tmp_path) -> None:
    plots_dir = tmp_path / "plots"
    root_files_dir = tmp_path / "root_files"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        set_root_preferences(plots_dir=plots_dir, root_files_dir=root_files_dir)

    assert str(ROOT.PlottingUtils.GetPlotsBaseDir()) == str(plots_dir.resolve())

    canvas = ROOT.TCanvas("test_canvas", "", 100, 100)
    ROOT.PlottingUtils.SaveFigure(
        canvas, "smoke", "calibration", ROOT.PlotSaveOptions.kLINEAR
    )

    expected = plots_dir / "calibration" / "smoke.png"
    assert expected.exists(), f"Expected {expected} to be written"


def test_set_plots_base_dir_strips_trailing_slash(tmp_path) -> None:
    target = tmp_path / "with_slash"
    ROOT.PlottingUtils.SetPlotsBaseDir(str(target) + "/")
    assert str(ROOT.PlottingUtils.GetPlotsBaseDir()) == str(target)


def test_set_root_preferences_warns_without_dirs() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        set_root_preferences()
    messages = [str(w.message) for w in caught]
    assert any("plots_dir" in m for m in messages)
    assert any("root_files_dir" in m for m in messages)
