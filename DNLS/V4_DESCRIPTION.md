# V4 deposit note for fig8/fig9 relabeling

The V3 deposit named two long-time PDFs `fig8_ratio_collapse.pdf` and
`fig9_T1e6_saturation.pdf`, but no dedicated producing script for the original
`fig9_T1e6_saturation.pdf` survived in the repository history or on the local
workstation.

For the V4 documentation patch, the cheapest science-equivalent fix is to reuse
and relabel the existing local long-time analysis figures instead of claiming a
fresh dedicated T=10^6 figure-generation pipeline:

- `fig8_ratio_collapse.pdf` is the relabeled export of `fig_long_alpha_fit.png`.
- `fig9_T1e6_saturation.pdf` is the relabeled export of `fig_long_ipr_vs_t.png`.

Both relabeled PDFs were generated from the current repository copies of the
local PNG outputs and are included only to preserve the V3 figure-slot naming.
They should be described transparently in the Zenodo V4 record as relabeled
stand-ins, not as newly regenerated dedicated figure scripts.

Current provenance in this repository:

- source data: `ipr_vs_time.csv`
- analysis script: `analyze_long_time.py`
- local figure exports: `fig_long_alpha_fit.png`, `fig_long_ipr_vs_t.png`
- relabeled V4 slot files: `DNLS/fig8_ratio_collapse.pdf`,
  `DNLS/fig9_T1e6_saturation.pdf`

If a true T=10^6-specific replacement is later required, rerun
`dnls_long_time.py` with a matching `ipr_vs_time.csv` horizon and regenerate the
long-time analysis figures before replacing these stand-ins.
