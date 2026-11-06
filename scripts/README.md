# Scripts

| Script | What it does |
|---|---|
| `prepare_figure_data.py` | Unpacks the archived exports in `../results/` into the folder layout read by the figure scripts (`E14`, `E15`, `E16`, `E17`, `E19`, `r1`, `r2`, `r3`) |
| `make_figures_JPR.py` | Main-text Figures 1–5 of the JPR paper, with the source data of every panel |
| `make_figures_JPR_SI.py` | Supporting-information Figures S1–S8, with the source data of every panel |
| `make_figure_S22.py` | Supporting-information Figure S9 (named S22 in an earlier draft), from the Experiment 17 export |
| `style.py` | Shared plotting style, imported by the two `make_figures_JPR` scripts |
| `numbers_ledger_JPR.csv` | Every value quoted in the main text, with the export and file it comes from; rows corrected on 6 Nov 2026 carry a note |
| `export_reference_trajectories.py` | Writes `../reference/reference_trajectories.json` for checking other implementations |

## Redrawing the figures of the paper

From the repository root, with the package installed (`pip install -e ".[notebooks]"`):

```bash
python -I scripts/prepare_figure_data.py figure_data
python -I scripts/make_figures_JPR.py    figure_data figures_main
python -I scripts/make_figures_JPR_SI.py figure_data figures_si
python -I scripts/make_figure_S22.py     figure_data/E17 figures_s9
```

Each figure is written as PDF, 600-dpi PNG and LZW-compressed TIFF, with one `_source.csv` per panel. Figures 1C, S1 and S2 are recomputed from the model with the `periomemory` package; every other panel is read from the exports. The whole run takes about a minute.

## Where each figure comes from

| Figure | Export |
|---|---|
| 1C, S1, S2 | recomputed with `periomemory` 1.0.0 |
| 1D, 2A–D, 2F, 5A; S1F, S4, S5, S6 | Experiment_14_Revision |
| 2E, 2F, 4C, 5B; S4, S6, S8B | Experiment_15_Model_Revision |
| 3, 4C, 5B; S4–S6, S7A | Experiment_16_Revision_Figures |
| 4A, 4B, 4D–F, 5D; S9 | Experiment_17_Resolution_Community_Memory |
| 5C; S7B, S8A | Experiment_19_Debridement_From_Treated_State |
| S3 | Robustness_R1–R3 |
