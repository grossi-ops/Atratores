Good move to audit before pushing. Let me fetch the Zenodo deposit, then cross-reference its contents against the public claims on the card and the requirements implied by the paper.Audit done. The deposit is a 3-version chain on concept DOI **10.5281/zenodo.20026942**:

| Version | DOI | Date | What's in it |
|---|---|---|---|
| V1 | 20026943 | May 4 | `dnls_nbonacci.py`, `TribonacciDNLS.lean`, `nbonacci_dnls_paper.tex`, `nbonacci_dnls_paper.pdf` |
| V2 | 20032853 | May 5 | `generate_figures.py`, `nbonacci_dnls_paper_with_figs.tex/.pdf`, `fig1`–`fig5` (.png) |
| V3 | 20062492 | May 7 (today) | `paper_v3.pdf`, `fig1`–`fig5` (.pdf), **`fig6_d2_natural`, `fig7_fss_T1e4`, `fig8_ratio_collapse`, `fig9_T1e6_saturation`** (.pdf) |

The Instagram card cites V3's DOI correctly, and every file in your local `~/Desktop/dnls/` that overlaps with the deposit matches byte-for-byte (md5 verified). `TribonacciDNLS.lean` has no `sorry`/`axiom`/`admit` in code, so the "Lean 4 verified, no sorry" claim is consistent with the file. So far so good.

The deposit does **not** cover all bases, though. Here's the gap matrix:

**What the latest paper (paper_v3.pdf) implies must exist, vs. what's actually deposited:**

| Asset that paper_v3 implies | In deposit? | Where it lives now |
|---|---|---|
| `.tex` source for paper_v3 | ❌ no | nowhere — V2's `nbonacci_dnls_paper_with_figs.tex` is older and only references figs 1–5 |
| Code producing fig6 (`d2_natural`) | ❌ no | nowhere |
| Code producing fig7 (`fss_T1e4`) — finite-size scaling at T=10⁴ | ❌ no | nowhere |
| Code producing fig8 (`ratio_collapse`) | ❌ no | nowhere |
| Code producing fig9 (`T1e6_saturation`) — T=10⁶ run | ❌ no | nowhere; `dnls_long_time.py` (local) is the only long-time integrator and it's set to T=10⁵, not 10⁶ |
| Numerical data behind the figures (CSV) | ❌ no | nowhere |
| `dnls_long_time.py` (the file we just resolved) | ❌ no | local only |
| `FoldEvents.lean` | ❌ no | local only — uploaded today, never deposited |

**Other coherence issues worth flagging:**

- The **V2 `.tex` source** says *"Our integration runs to T = 50"* (line 495) and lists "longer time evolution" as **the most urgent open question** in §7. paper_v3.pdf clearly goes beyond this (fig9 is at T = 10⁶), but the .tex producing it isn't deposited, so a reader can't reconcile the two.
- **License mismatch.** The deposit is CC-BY-NC-**ND**-4.0. The "ND" (no derivatives) clause is the wrong license for *code* and *Lean proofs* — readers can't legally fork, modify, or even re-run the .py/.lean files. Conventional choice is CC-BY for the paper and a code-appropriate license (MIT, Apache-2.0, or BSD) for the .py/.lean files. You can keep ND on the prose if you want, but the code needs a derivative-allowing license to be reproducible in any meaningful sense.
- **DOI staleness in `dnls_long_time.py` docstring.** It cites `10.5281/zenodo.20026943` (V1), not the concept DOI `10.5281/zenodo.20026942` (which auto-resolves to latest) and not V3. Either is fine; the current value points to a deprecated description.
- **No README/MANIFEST in any version.** Three versions, different contents, no map of which file lives where or how to reproduce a specific figure. This is the single highest-leverage thing to add.
- **Numerical drift V1 → V2/V3.** V1 description says "IPR drops ~60% at λ=1.5"; V2/V3 say "~57%". Defensible (number got refined), but worth a note.

**My read:** the deposit gives a reader the *finished paper* and the *original code skeleton*, but it does **not** give them the *long-time evidence pipeline*. Anyone trying to reproduce fig7–fig9 from the deposit alone will hit a brick wall. The deposit's own description in V3 even waves at this — *"Companion code (dnls_nbonacci.py) and Lean 4 verification file (TribonacciDNLS.lean) are included in this deposit – see V1"* — but doesn't acknowledge the missing long-time / FSS scripts.

**Cheapest fixes, in order (Cohn-style):**

1. **Write a `README.md`/`MANIFEST.md`** that maps figure → script → data → version, and add it as a small "V4 documentation patch" deposit. ~30 min, immediately raises the deposit from "polished but opaque" to "navigable."
2. **Find or recreate the scripts that produced fig6–fig9.** They had to exist somewhere on your machine when you generated those PDFs in the last 2 days. If you can locate them, deposit them. If not, the path of least resistance is to *recreate* fig9 from `dnls_long_time.py` (T=10⁵ → 10⁶ flag flip + a fit script) so at least the long-time figure has a deposited generator.
3. **Deposit the paper_v3 source** (`.tex` + bibliography + the `.bbl` if needed). Same V4.
4. **Split or relax the license** for the code (.py/.lean) to a derivative-allowing one. This is a metadata-only change.
5. **Add `dnls_long_time.py` and (eventually) the long-time CSVs** to the V4 deposit.
6. **Update the docstring DOI** in `dnls_long_time.py` to the concept DOI `10.5281/zenodo.20026942`.

If you want, I can draft the README + a "what V4 should look like" file list right now — that's the highest-leverage of the six and unblocks everything else.

Sources:
- [Zenodo record V3 (20062492)](https://zenodo.org/records/20062492)
- [Zenodo record V2 (20032853)](https://zenodo.org/records/20032853)
- [Zenodo record V1 (20026943)](https://zenodo.org/records/20026943)
- [Concept DOI parent (20026942)](https://zenodo.org/records/20026942)
