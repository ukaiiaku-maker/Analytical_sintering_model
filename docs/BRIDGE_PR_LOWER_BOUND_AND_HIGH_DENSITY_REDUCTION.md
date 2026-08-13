# Bridge between PR lower bound and high-density reduction

Inputs were recovered selectively from `results/Archive.zip`; no historical result tree was restored. Both anchors were deterministically reconstructed from candidate IDs, seed 20260813, and the recorded decoder.

Seventeen block transplants, 90 block morphs, and 6,400 seeded optimizer evaluations were run. Coarse screening appeared to produce bridge candidates, but mandatory 30-minute-step reconfirmation rejected all 38 apparent Tier B cases. The coarse 2-hour integrator inflated anchor A's median reduction from 11.3% to about 63%.

The exact result remains split: anchor 155976 retains the complete 50 C window but only 11.3% reduction; anchor 4412 retains strong reduction but no lower boundary. No Tier A/B bridge was accepted. This is a numerical and mechanistic negative result, not validation.

Bridge-specific tests pass. The repository-wide suite reports 60 missing-fixture failures because older tracked result folders were intentionally archived and deleted; those artifacts were not restored or included in this commit.
