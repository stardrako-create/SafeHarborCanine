# 05_SHIP — Safe Harbor Integration Prioritization

Per-locus scoring that combines the ATAC-seq, RRBS, and Hi-C tracks from
`04_tracks_processadas/` into a single safe-harbor suitability score genome-
wide. The Hi-C layer enters this score at a deliberately reduced weight
relative to the ATAC/RRBS Mother Tracks, since it derives from a single
individual of a different breed/tissue rather than the 71-dog cohort (see
the top-level README, "Data layers").

Not yet populated — this phase begins once the Hi-C contact matrix and TAD
boundary calls in `04_tracks_processadas/ROS_Cfam_1.0/HiC/` are complete.
