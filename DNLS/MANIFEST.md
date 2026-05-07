# DNLS V4 manifest

Checksums use SHA-256. Paths are relative to `DNLS/`.

| Path | SHA-256 | Provenance |
|---|---|---|
| `LICENSE-CODE` | `15c92596aca12afd3b391bbb49d908b60581236cdd894482d1f6b11fdf8e6e0e` | Imported from the split-license DNLS update branch. |
| `LICENSE-PROSE` | `6275df7d855defe5624ddc499bc5f3d8d68d4a1b42e137238c367319bd4967f8` | Imported from the split-license DNLS update branch. |
| `README.md` | `7d25de1a1049390d70a075618475884fed74fb0e1a30675ecb61fb59693aba8f` | Authored in this PR as the top-level V4 deposit guide. |
| `V4_DESCRIPTION.md` | `355c11ad3483d113844bc01276eb3fafde114fe4f912b2ed6640937c0d01455a` | Imported from the V4 metadata branches and trimmed to the canonical shipped assets. |
| `code/analyze_long_time.py` | `64f1863efd917e3a7a1660b45a50126fd1b67306411467cab5b4b00cdea7ec13` | Mirrored from repo-root analyze_long_time.py with DNLS-local default data/figure paths. |
| `code/dnls_long_time.py` | `d427c0b22776b85c610c71c8a66aa7b2f860e8cfa6f2dc33dd63af65c3f11f0e` | Mirrored from repo-root dnls_long_time.py with default output redirected into DNLS/data. |
| `code/dnls_nbonacci.py` | `84610ede34504d23568511de782ed6442afb17daef1b930bc463fe16722ea73e` | Moved from DNLS root with history preserved. |
| `code/generate_figures.py` | `a701e0ab5c038f7f98edf51eadf02ab9b32908594dbca050c35defa9afb18188` | Moved from DNLS root with output paths retargeted to DNLS/figures. |
| `code/paper_figures.py` | `8152faa702a652d49433f9439717b2d2f100442434a12a37e6c9f4577560b6d6` | Added in this PR as the deposit-facing wrapper for d2_natural_lengths.py and fss_analyze.py. |
| `data/ipr_vs_time.csv` | `4d49e31fa89b3c5e3c0414ba533628b9a4ccf6d42a1ef29e6fb97c4d53dd73b2` | Mirrored from the repo-root long-time run output produced by dnls_long_time.py. |
| `data/spreading_exponents.csv` | `f972a55dae34fc47deb29bd7fc2eebf381739777a98c4bad761d05e7ab3ac1e5` | Mirrored from the repo-root analyze_long_time.py output. |
| `figures/fig1_chain_structure.pdf` | `1cc7e484d83b9ce1a06295c156c459ba44a6b7e2d70a6122c6d0f4922f1ae55f` | Generated in this PR by code/generate_figures.py after the folder move. |
| `figures/fig2_eigenstates.pdf` | `8b32602f4d585992584e1f28b9d9091cad67e4c5cdbab2128756b614bbe78821` | Generated in this PR by code/generate_figures.py after the folder move. |
| `figures/fig3_ipr_vs_lambda.pdf` | `2a5114f6dd337d9a7f18a10001f6116b7571ff738a1a934af5ca610deac72679` | Generated in this PR by code/generate_figures.py after the folder move. |
| `figures/fig4_ipr_ratio.pdf` | `ac8f194eeb457e467c0df8430c65dbaad51e7486ebd7b0ba0431f72ecf663023` | Generated in this PR by code/generate_figures.py after the folder move. |
| `figures/fig5_substitution_tree.pdf` | `b7459548fe87daeab87b50fedfe08c1901324e84ceaa382667f96ddc48b7d145` | Generated in this PR by code/generate_figures.py after the folder move. |
| `figures/fig6_d2_natural.pdf` | `824f9e3628aa0ec7a0321b4f043354263df26bdfa508a2192ff46fa2e64073e2` | Imported from the deposit-facing figures branch (rendered from d2_natural_lengths.py). |
| `figures/fig7_fss_T1e4.pdf` | `4d4ac56bed9965160c88ce33710c8cfcf9e83c8e241315a6f070731bc4a6b3eb` | Imported from the deposit-facing figures branch (rendered from fss_analyze.py). |
| `figures/fig8_ratio_collapse.pdf` | `583caf5d9165e1c239d7f3033c0de9bf4c96e5f51cff3f069af9ee487c5abf5e` | Imported from the deposit-facing figures branch (rendered from fss_analyze.py). |
| `figures/fig9_T1e6_saturation.pdf` | `72c608268deb730a43cf91e48f8e0a2f1566cb95893ec612701240f4b284afd5` | Imported from the deposit-facing figures branch; regeneration matched the canonical long-time PNG output. |
| `figures/fig_A_d2_scaling.png` | `ce316b55477bafc40dd47e3262512498dda877df459784efd5dd885982ae7fa9` | Imported from the deposit-facing figures branch; byte-identical to figures/d2_natural_lengths.png. |
| `figures/fig_B_nstability.png` | `6750a65abcacc0fd829787a181c705f90841ff7d318c1a157fdcd18b3f271b6a` | Imported from the deposit-facing figures branch; byte-identical to fig_fss_alpha_vs_N.png. |
| `figures/fig_C_inversion.png` | `b192267d2462f8a0bec307493131293fc902c80f33c02ffe1091b804dd34cb90` | Imported from the deposit-facing figures branch; byte-identical to fig_fss_D2.png. |
| `figures/fig_D_homogenization.png` | `904d01d0b7057398147ad124c9876528523c7b1dce28466d596d948a01b5381d` | Imported from the deposit-facing figures branch; byte-identical to fig_fss_tsat_vs_N.png. |
| `lean/FoldEvents.lean` | `d4803c2cf85567f6f425e2dfadb5c1da0ac4a10521aee1cfd1f4d25948cf5e84` | Moved from DNLS root; canonical upstream remains TOTOGT/AXLE. |
| `lean/TribonacciDNLS.lean` | `749f9aefc1f1c55e8f70ad878ba96b16d4e4ebaa9e3c30ea464febf8914b2487` | Moved from DNLS root; canonical upstream remains TOTOGT/AXLE. |
| `paper/nbonacci_dnls_paper_v3.pdf` | `7ac68c1f544e591b96bfb532f78a58c7179b13d05a5409980235d736f3019af7` | Moved/renamed from DNLS/nbonacci_dnls_paper_with_figs.pdf. |
| `paper/nbonacci_dnls_paper_v3.tex` | `f504b3d1b589a4b21f2acbfd536adf6db28aad83cf3094a1e395adaba5fbd3d8` | Moved/renamed from DNLS/nbonacci_dnls_paper_with_figs.tex with graphicspath updated for paper/. |
| `paper/refs.bib` | `5f68ff7767df54da091a00b1f48a692eecf10f03e44123b9f59430d94595331b` | Added in this PR as a bibliography-placeholder companion because the TeX source keeps an inline thebibliography block. |
| `MANIFEST.md` | `self-referential` | Authored in this PR; checksum omitted because editing the file changes its own digest. |
