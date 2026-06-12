import array as _array
import os
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import ROOT

_DEFAULT_CACHE_DIR = "df_cache"

_root_files_base = "root_files"

PathLike = Union[str, os.PathLike]


def set_root_files_base_dir(dir: PathLike) -> None:
    """Set the base directory used by :func:`open_for_reading` and
    :func:`open_for_writing` for relative subpaths.

    The argument is resolved to an absolute path via
    ``Path(dir).expanduser().resolve()`` and stored as a string. When the
    C++ library is loadable, the resolved value is also forwarded to
    :cpp:func:`IO::SetRootFilesBaseDir` so Python and C++ share state.
    """
    global _root_files_base
    resolved = str(Path(dir).expanduser().resolve())
    _root_files_base = resolved
    try:
        from analysis_utilities import load_cpp_library

        load_cpp_library().IO.SetRootFilesBaseDir(resolved)
    except (ImportError, RuntimeError, OSError):
        pass


def get_root_files_base_dir() -> str:
    """Return the current root-files base directory (absolute path string)."""
    return _root_files_base


def _join_path(subpath: str) -> str:
    if ROOT.gSystem.IsAbsoluteFileName(subpath):
        return subpath
    return str(ROOT.gSystem.ConcatFileName(_root_files_base, subpath))


def open_for_reading(subpath: str) -> "ROOT.TFile":
    """Open ``<root_files_base>/<subpath>`` for reading.

    Absolute ``subpath`` values are passed through untouched. Returns a
    ``ROOT.TFile`` opened in ``"READ"`` mode.
    """
    full = _join_path(subpath)
    return ROOT.TFile(full, "READ")


def open_for_writing(subpath: str, mode: str = "RECREATE") -> "ROOT.TFile":
    """Open ``<root_files_base>/<subpath>`` for writing, creating parents.

    The parent directory of the resolved path is created via
    ``ROOT.gSystem.mkdir(..., True)``. Absolute ``subpath`` values are
    passed through untouched. Returns a ``ROOT.TFile`` opened in
    ``mode`` (default ``"RECREATE"``).
    """
    full = _join_path(subpath)
    parent = ROOT.gSystem.DirName(full)
    ROOT.gSystem.mkdir(parent, True)
    return ROOT.TFile(full, mode)


_TYPE_MAP = {
    "Float_t": ("f", np.float32),
    "float": ("f", np.float32),
    "Double_t": ("d", np.float64),
    "double": ("d", np.float64),
    "Long64_t": ("q", np.int64),
    "long long": ("q", np.int64),
    "ULong64_t": ("Q", np.uint64),
    "unsigned long long": ("Q", np.uint64),
    "Int_t": ("i", np.int32),
    "int": ("i", np.int32),
    "UInt_t": ("I", np.uint32),
    "unsigned int": ("I", np.uint32),
    "Short_t": ("h", np.int16),
    "short": ("h", np.int16),
    "UShort_t": ("H", np.uint16),
    "unsigned short": ("H", np.uint16),
    "UChar_t": ("B", np.uint8),
    "unsigned char": ("B", np.uint8),
}


def _cache_key(root_files, tree_name):
    """Build a cache filename from the ROOT file paths and tree name."""
    base = "_".join(
        os.path.splitext(os.path.basename(p))[0] for p in root_files)
    return f"{base}_{tree_name}"


def _newest_mtime(root_files):
    """Return the newest modification time among root_files."""
    return max(os.path.getmtime(p) for p in root_files)


def load_tree_data(
    root_files,
    tree_name="features",
    scalar_branches=None,
    array_branch=None,
    max_events=None,
    cache_dir=_DEFAULT_CACHE_DIR,
):
    """Load TTree data into numpy arrays and pandas DataFrame.

    Results are cached as pickle (DataFrame) and .npy (waveform array)
    files inside cache_dir.  The cache is invalidated when
    any source ROOT file is newer than the cached file.

    Parameters
    ----------
    root_files : str or list of str
        Path(s) to ROOT file(s).
    tree_name : str
        Name of the TTree to read.
    scalar_branches : list of str or None
        Scalar branch names to load. If None, auto-detects all
        scalar (non-array) branches.  Caching requires all branches
        (i.e. scalar_branches=None).
    array_branch : str or None
        Name of a TArrayF/TArrayS branch to load as a 2-D numpy array.
    max_events : int or None
        Maximum number of events to load.
    cache_dir : str or None
        Directory for cached DataFrames.  Set to None to disable
        caching.  Defaults to "df_cache" (relative to cwd).

    Returns
    -------
    features_df : pandas.DataFrame
        Scalar branch data.
    waveforms : numpy.ndarray or None
        2-D array (n_events, n_samples) if array_branch is given, else None.
        Only returned if array_branch is not None.
    """
    if isinstance(root_files, str):
        root_files = [root_files]

    use_cache = cache_dir is not None
    if use_cache and scalar_branches is not None:
        raise ValueError(
            "Caching requires loading all branches (scalar_branches=None). "
            "Either omit scalar_branches or set cache_dir=None.")

    if use_cache:
        os.makedirs(cache_dir, exist_ok=True)
        key = _cache_key(root_files, tree_name)
        pkl_path = os.path.join(cache_dir, f"{key}.pkl")
        npy_path = os.path.join(cache_dir, f"{key}.npy")

        if os.path.exists(pkl_path):
            src_mtime = _newest_mtime(root_files)
            if src_mtime <= os.path.getmtime(pkl_path):
                print(f"Loading cached DataFrame: {pkl_path}")
                df = pd.read_pickle(pkl_path)
                if array_branch:
                    wf = np.load(npy_path) if os.path.exists(
                        npy_path) else None
                    return df, wf
                return df

    chain = ROOT.TChain(tree_name)
    for path in root_files:
        if chain.Add(path) == 0:
            raise FileNotFoundError(f"Could not add {path} to TChain")

    n_total = chain.GetEntries()
    if n_total == 0:
        raise ValueError(f"TChain is empty (files: {root_files})")

    if scalar_branches is None:
        scalar_branches = []
        branch_list = chain.GetListOfBranches()
        for i in range(branch_list.GetEntries()):
            br = branch_list.At(i)
            name = br.GetName()
            if name != array_branch:
                scalar_branches.append(name)

    if max_events is not None:
        n_to_read = min(n_total, max_events)
    else:
        n_to_read = n_total

    chain.SetBranchStatus("*", 0)
    for name in scalar_branches:
        br = chain.GetBranch(name)
        if br is None:
            raise ValueError(
                f"Branch '{name}' not found in tree '{tree_name}'")
        chain.SetBranchStatus(name, 1)
    if array_branch:
        if chain.GetBranch(array_branch) is None:
            raise ValueError(
                f"Branch '{array_branch}' not found in tree '{tree_name}'")
        chain.SetBranchStatus(array_branch, 1)

    buffers = {}
    np_dtypes = {}
    for name in scalar_branches:
        leaf = chain.GetBranch(name).GetLeaf(name)
        type_name = leaf.GetTypeName()

        entry = _TYPE_MAP.get(type_name)
        if entry is None:
            print(
                f"Warning: skipping branch '{name}' (unsupported type '{type_name}')"
            )
            continue
        typecode, dt = entry
        buf = _array.array(typecode, [0])
        np_dtypes[name] = dt

        chain.SetBranchAddress(name, buf)
        buffers[name] = buf

    if array_branch:
        br_class = chain.GetBranch(array_branch).GetClassName()
        if "TArrayS" in br_class:
            arr_obj = ROOT.TArrayS()
            arr_dtype = np.int16
        else:
            arr_obj = ROOT.TArrayF()
            arr_dtype = np.float32
        chain.SetBranchAddress(array_branch, arr_obj)
        chain.GetEntry(0)
        wf_size = arr_obj.GetSize()
        waveforms = np.empty((n_to_read, wf_size), dtype=arr_dtype)
    else:
        waveforms = None

    scalar_data = {
        name: np.empty(n_to_read, dtype=np_dtypes[name])
        for name in buffers
    }

    read_count = 0
    entry_idx = 0

    while read_count < n_to_read:
        if entry_idx >= n_total:
            break

        nb = chain.GetEntry(entry_idx)
        if nb <= 0:
            entry_idx += 1
            continue

        for name, buf in buffers.items():
            scalar_data[name][read_count] = buf[0]

        if array_branch:
            waveforms[read_count] = np.frombuffer(arr_obj.GetArray(),
                                                  dtype=arr_dtype,
                                                  count=wf_size)

        read_count += 1
        entry_idx += 1

    if read_count < n_to_read:
        for name in scalar_data:
            scalar_data[name] = scalar_data[name][:read_count]
        if waveforms is not None:
            waveforms = waveforms[:read_count]

    features_df = pd.DataFrame(scalar_data)

    if use_cache:
        features_df.to_pickle(pkl_path)
        print(f"Cached DataFrame to {pkl_path}")
        if waveforms is not None:
            np.save(npy_path, waveforms)
            print(f"Cached waveforms to {npy_path}")

    if array_branch:
        return features_df, waveforms
    else:
        return features_df
