# Sensitivity report (redesigned 2026-09-11 - see script docstring for why)

## 1. Threshold sensitivity (survivor SET per threshold, baseline weights)

| threshold | # survivors |
|---|---|
| p10 | 22 |
| p15 | 19 |
| p20 | 17 |
| p25 | 13 |
| p35 | 6 |
| p45 | 4 |
| p55 | 1 |

## 2. Ranking stability under weight perturbation, held at p25 (26 weight configs)

Threshold is fixed here, so the survivor SET is fixed too (weights cannot change hard_veto) - this isolates whether the SCORE ORDERING among that fixed field is stable under weight uncertainty, which p55's single-survivor field could not test.

| locus | genes | mean rank | rank range | #1 in N/total configs |
|---|---|---|---|---|
| NC_051833.1:23112498-23172320 | LOC111093212/LOC111093135 | 1.2 | 1-2 | 22/26 |
| NC_051826.1:29235878-29299159 | LOC100684181/LOC111091654 | 2.0 | 1-3 | 4/26 |
| NC_051805.1:7072137-7132579 | LOC111090579/LOC100685067 | 3.3 | 2-4 | 0/26 |
| NC_051823.1:40699370-40753962 | LOC119864393/LOC119864525 | 3.7 | 3-5 | 0/26 |
| NC_051843.1:104147730-104215864 | LOC119868710/IGSF1 | 5.3 | 4-6 | 0/26 |
| NC_051820.1:51445768-51508808 | LOC119874730/LOC119874609 | 5.6 | 5-6 | 0/26 |
| NC_051811.1:25873242-25943461 | LOC111096660/LOC111096885 | 7.0 | 7-7 | 0/26 |
| NC_051815.1:7585907-7653876 | LOC100687242/LOC106559451 | 8.2 | 8-9 | 0/26 |
| NC_051826.1:55859005-55918918 | LOC119865267/LOC119865154 | 8.8 | 8-9 | 0/26 |
| NC_051831.1:39677324-39739751 | ANO2/NTF3 | 10.5 | 10-11 | 0/26 |
| NC_051806.1:17055808-17110317 | PTCHD3/ANKRD26 | 10.5 | 10-11 | 0/26 |
| NC_051827.1:1109045-1172665 | SYNDIG1/LOC100688361 | 12.0 | 12-12 | 0/26 |
| NC_051837.1:7790832-7842018 | LOC111093705/IMPG2 | 13.0 | 13-13 | 0/26 |