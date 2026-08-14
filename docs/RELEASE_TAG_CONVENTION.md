# Release tag convention

This repository uses annotated, immutable Git tags for tested scientific
handoff milestones. Tags identify a reproducible repository state; they do not
constitute model validation or a paper-ready calibration.

## Naming

The final synthesis series uses this namespace:

```text
final-sintering-mechanism-synthesis-YYYY-MM-DD
```

Use the UTC handoff date. If a milestone needs review before final publication,
append `-rcN`. If a published milestone requires a correction, create a new
`-revN` tag rather than moving or reusing the original tag:

```text
final-sintering-mechanism-synthesis-2026-08-14-rc1
final-sintering-mechanism-synthesis-2026-08-14-rev1
```

## Tag requirements

- Create an annotated tag, never a lightweight tag.
- Tag a commit that has already been committed and pushed on its provenance
  branch.
- Require the applicable tests and compilation checks to pass before tagging.
- Include the handoff report, branch index, compact evidence manifests, and
  non-claims applicable to the milestone in the tagged commit.
- Summarize the evidence scope and principal non-claim in the annotation.
- Never force-update, delete, or reuse a tag after it has been pushed. Publish a
  new dated or `-revN` tag for corrections.
- Push tags explicitly; do not use `git push --tags`, which can publish
  unrelated local tags.

## Publication and verification

For a tested handoff commit at `HEAD`:

```bash
git tag -a final-sintering-mechanism-synthesis-YYYY-MM-DD \
  -m "Final exact mechanism synthesis and source-data handoff; conditional Tier B, not validation"
git push origin final-sintering-mechanism-synthesis-YYYY-MM-DD
git show --no-patch final-sintering-mechanism-synthesis-YYYY-MM-DD
git ls-remote --tags origin final-sintering-mechanism-synthesis-YYYY-MM-DD \
  'final-sintering-mechanism-synthesis-YYYY-MM-DD^{}'
```

An annotated tag produces both a tag-object reference and a peeled commit
reference in the remote listing. The peeled commit must equal the intended,
tested handoff commit.

## Current milestone

The first tag in this series is
`final-sintering-mechanism-synthesis-2026-08-14`. It marks the final exact
mechanism synthesis, equation-audit QC, figure source-data package, repository
handoff, and this convention. Candidate 693168 remains conditional Tier B, not
validation.
