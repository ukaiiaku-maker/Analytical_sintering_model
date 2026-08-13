# Archived result fixture policy

Historical result folders may be compressed into `results/Archive.zip` and intentionally removed from the working tree. Tests that read those persisted artifacts are marked `requires_archived_results`; they remain runnable after selective restoration, but are excluded with `-m "not requires_archived_results"` during current science development. Unit tests and current-branch tests are never marked by this policy. Archives are not extracted wholesale and archived deletions are not mixed into science commits.
