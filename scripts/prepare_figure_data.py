"""Unpack the archived exports in ../results into the folder layout read by the figure scripts.

Usage: python -I prepare_figure_data.py [data folder]      (default: ./figure_data)

Creates <data folder>/E14, E15, E16, E17, E19, r1, r2, r3, each holding the contents of one
export, as expected by make_figures_JPR.py and make_figures_JPR_SI.py. The archives are not modified.
"""
import os
import shutil
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")
EXPORTS = {
    "E14": "Experiment_14_Revision_EXPORT.zip",
    "E15": "Experiment_15_Model_Revision_EXPORT.zip",
    "E16": "Experiment_16_Revision_Figures_EXPORT.zip",
    "E17": "Experiment_17_Resolution_Community_Memory_EXPORT.zip",
    "E19": "Experiment_19_Debridement_From_Treated_State_EXPORT.zip",
    "r1": "Robustness_R1_alpha_epsilon_memory_map_EXPORT.zip",
    "r2": "Robustness_R2_multivariate_structural_sample_EXPORT.zip",
    "r3": "Robustness_R3_functional_form_EXPORT.zip",
}


def unpack(archive, target):
    """Extract `archive` so that its files sit directly in `target`, whether or not the archive
    wraps them in a single top-level folder."""
    tmp = target + "_tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    shutil.rmtree(target, ignore_errors=True)
    with zipfile.ZipFile(archive) as z:
        z.extractall(tmp)
    entries = os.listdir(tmp)
    inner = os.path.join(tmp, entries[0]) if len(entries) == 1 and os.path.isdir(os.path.join(tmp, entries[0])) else tmp
    shutil.move(inner, target)
    shutil.rmtree(tmp, ignore_errors=True)


def main():
    out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "figure_data")
    os.makedirs(out, exist_ok=True)
    for key, name in EXPORTS.items():
        unpack(os.path.join(RESULTS, name), os.path.join(out, key))
        print(f"{key:4s} <- {name}")
    print("figure data in", out)


if __name__ == "__main__":
    main()
