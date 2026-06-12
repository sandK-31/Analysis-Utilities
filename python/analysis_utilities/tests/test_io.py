"""Tests for analysis_utilities.io read/write helpers."""

from __future__ import annotations

import warnings

import pytest

ROOT = pytest.importorskip("ROOT")

try:
    from analysis_utilities import (get_root_files_base_dir, load_cpp_library,
                                      open_for_reading, open_for_writing,
                                      set_root_files_base_dir,
                                      set_root_preferences)

    load_cpp_library()
except (RuntimeError, OSError) as exc:
    pytest.skip(f"C++ library not loadable: {exc}", allow_module_level=True)


def test_open_for_writing_creates_nested_parents(tmp_path) -> None:
    set_root_files_base_dir(tmp_path)
    fout = open_for_writing("deep/nested/dir/file.root")
    try:
        assert fout.IsOpen()
    finally:
        fout.Close()
    expected = tmp_path / "deep" / "nested" / "dir" / "file.root"
    assert expected.exists()


def test_round_trip_read_write(tmp_path) -> None:
    set_root_files_base_dir(tmp_path)

    fout = open_for_writing("round/trip.root")
    h = ROOT.TH1F("h", "h", 10, 0.0, 1.0)
    h.SetBinContent(5, 42.0)
    h.Write()
    fout.Close()

    fin = open_for_reading("round/trip.root")
    try:
        assert fin.IsOpen()
        h_read = fin.Get("h")
        assert h_read is not None
        assert h_read.GetBinContent(5) == 42.0
    finally:
        fin.Close()


def test_absolute_path_passthrough(tmp_path) -> None:
    other_base = tmp_path / "ignored_base"
    set_root_files_base_dir(other_base)

    abs_target = tmp_path / "abs_target.root"
    fout = open_for_writing(str(abs_target))
    try:
        assert fout.IsOpen()
    finally:
        fout.Close()
    assert abs_target.exists()
    # Absolute path should bypass the base dir entirely; base dir was
    # never touched, so it shouldn't have been created.
    assert not other_base.exists()


def test_set_root_preferences_propagates_root_files_base(tmp_path) -> None:
    plots_dir = tmp_path / "plots"
    root_files_dir = tmp_path / "root_files"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        set_root_preferences(plots_dir=plots_dir, root_files_dir=root_files_dir)
    assert get_root_files_base_dir() == str(root_files_dir.resolve())
    assert (
        str(ROOT.IO.GetRootFilesBaseDir()) == str(root_files_dir.resolve())
    )


def test_set_root_files_base_dir_strips_trailing_slash(tmp_path) -> None:
    target = tmp_path / "no_slash"
    set_root_files_base_dir(str(target) + "/")
    assert get_root_files_base_dir() == str(target.resolve())
    assert str(ROOT.IO.GetRootFilesBaseDir()) == str(target.resolve())
