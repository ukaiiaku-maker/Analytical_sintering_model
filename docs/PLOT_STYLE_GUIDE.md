# Plot Style Guide

All manuscript figures use `plot_style.py` and DejaVu Sans with embedded
TrueType PDF fonts. Main labels are 11 pt, ticks and legends are 9 pt, panel
letters are 13 pt bold, and plotted lines are 2.2 pt. Final files are emitted
as vector PDF and 600-dpi PNG using deterministic stems.

Protocol colors are fixed: slow heating blue, fast heating vermillion,
high-temperature isothermal magenta, and two-step green. Positive, neutral,
harmful, and unattainable classes use green, gray, vermillion, and purple.
The q0/q1 colors are light blue/orange. TJ modes use blue, teal, orange, and
purple consistently across main and supplementary figures.

Notation follows the manuscript convention: relative density $\rho$, grain
size $G$ [nm], connected mean pore radius $\bar r_p^{\,c}$ [nm], connected
fine-pore fraction $f_{\mathrm{fine}}^c$, PR work $W_{\mathrm{PR}}$, persistent
junction state $X_J$, and TJ activity ratio
$\Lambda_{\mathrm{TJ}}/K_{\mathrm{TJ}}$. TJ population superscripts distinguish
total, pore, structural, constraint, relaxed, and pinned fractions.

Two-column figures use a 7.2-inch canvas; focused one-column figures use
approximately 4.3–5.5 inches. Multi-panel figures use constrained layout and
panel letters. Legends remain inside only when they do not cover data;
otherwise they move outside. Dimensionless bounded fractions use explicit
`[-]` labels where helpful. Logarithmic heating-rate coordinates are used only
for the decade-spanning rate axis.

Regenerate the package with:

```bash
python3 generate_paper_figures.py
python3 generate_supplement_figures.py
```
