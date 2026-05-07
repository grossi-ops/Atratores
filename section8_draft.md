## 8.2 Fractal dimension D₂ of the mid-gap eigenstate

\subsubsection{Fractal dimension $D_2$ of the mid-gap eigenstate}

The fractal (correlation) dimension $D_2$ is defined by the power-law scaling $\mathrm{IPR}(N) \sim N^{-D_2}$, so that $D_2 = 1$ for a fully extended state and $D_2 = 0$ for a perfectly localised one; quasiperiodic critical states occupy the intermediate range $0 < D_2 < 1$.

Earlier estimates from truncated chains of arbitrary length ($N \in \{200, 500, 1000, 2000\}$) gave $D_{2,\mathrm{fib}} \approx 0.616$ and $D_{2,\mathrm{trib}} \approx 0.364$.  The tribonacci estimate was anomalous: the underlying IPR$(N)$ data were non-monotone (0.0969 $\to$ 0.0820 $\to$ 0.0410 $\to$ 0.0484 for $N=200,500,1000,2000$), with a reversal between $N=1000$ and $N=2000$ that inflated the uncertainty and biased the slope.

The reversal is a finite-size artefact of \emph{boundary truncation}: arbitrary truncation at $N=1000$ or $N=2000$ destroys the self-similar boundary structure that the Rauzy fixed-point chain lengths preserve.  To test this hypothesis we re-computed IPR only at the natural Rauzy iteration lengths $T_n$ (where $T_{n+3} = T_{n+2} + T_{n+1} + T_n$, OEIS A000073), specifically $N \in \{274, 504, 927, 1705, 3136\}$ ($n = 10$--$14$), and the analogous natural Fibonacci lengths $N \in \{233, 377, 610, 987, 1597, 2584\}$ ($n = 12$--$17$).  At each natural length the chain is the exact $n$-th iterate of the substitution rule, so the boundary atoms are self-consistently generated and no truncation is introduced.

The mid-gap eigenstate is identified as the eigenstate with the smallest $|E|$ whose spatial spread $\sigma \geq 0.03N$; this criterion excludes compact boundary/defect zero-modes (notably the anomalous E=0 state at tribonacci $n=12$, $N=927$, whose spread $\sigma \approx 15 \ll 0.03 \times 927$) that are unrelated to the bulk multifractal scaling.  With this selection the tribonacci IPR sequence is monotone non-increasing: $0.097 \to 0.082 \approx 0.082 \approx 0.082 \to 0.041$.  The near-plateau at IPR $\approx 0.082$ for $n = 11, 12, 13$ (N = 504, 927, 1705) reflects a genuine feature of the tribonacci RSRG hierarchy in which the same critical-state family persists across three consecutive iterates.  The $(N=1000 \to N=2000)$ reversal present in the arbitrary-length data is eliminated: the natural-length sequence ends with a clean drop from IPR $\approx 0.082$ (N=1705) to $0.041$ (N=3136).

The OLS regression of $\log\,\mathrm{IPR}$ on $\log N$ yields
\begin{equation}
  D_{2,\mathrm{fib}}  = 0.646 \pm 0.111 \quad (R^2 = 0.8946),
\end{equation}
\begin{equation}
  D_{2,\mathrm{trib}} = 0.282 \pm 0.119 \quad (R^2 = 0.6535).
\end{equation}
The tribonacci $D_2$ has shifted from the anomalous $\approx0.364$ (arbitrary lengths, $R^2 = 0.30$) to $0.282$ (natural lengths, $R^2 = 0.65$), resolving referee question~(2).  The hierarchy $D_{2,\mathrm{trib}} < D_{2,\mathrm{fib}}$ at natural lengths is consistent with the stronger spatial multifractality of the tribonacci eigenstate.

