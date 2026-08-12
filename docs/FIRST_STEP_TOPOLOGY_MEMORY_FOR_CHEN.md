# First-Step Topology Memory for Chen

This checkpoint registers six physical memory families and evaluates 768 parameter candidates against exact matched-density E0021/E0142 first-step states. The bounded state set covers G0=100/300/600 nm, T1=1300/1500 C, and switch density 0.80/0.88. There are 9,216 candidate/state rows; 3,816 pass diagnostic divergence/persistence thresholds.

This is **not a dynamic family validation**. Candidate memory is currently a transparent projection of measured pore/grain/PR/XJ differences with a named exponential persistence time. Because all six families share the same observable base state, their divergence distributions are identical. The screen identifies parameter regions worth implementing, but cannot rank the families physically.

No new Tier A/B claim is made. Existing E0142 Tier B and E0021 Tier C windows are carried only as baselines. A next implementation must add one family-specific state at a time to the integrator and show that it changes migration without changing shared-state densification.
