# Selection-strength check: plan written before the runs (2026-11-16 21:39 UTC)

Question. In the deposited alternative-model analysis, brain size (B) held under the accelerating (superlinear)
social return and declined under linear and saturating returns. In those two models the average agent sits at or
above the 0.999 survival cap. Does B still decline under linear and saturating returns when survival is kept
below the cap?

Model. Deposited code, release v1.0.0 (childhood_first_simulation.py), with the 4 edits in patch.diff. With the
default switches the model is unchanged; this is verified by reproducing published replicates exactly.

Manipulation. One constant, baseline_offset = 3.5, subtracted from every agent's net fitness in every condition.
Chosen in advance by this rule: baseline mortality of about 10% per generation at the published Simulation 1
control endpoint (net fitness 5.66 - 3.5 = 2.16, survival 0.90), which is of the order of the mortality in the
Muthukrishna et al. model (6.5% at their default lambda = 1, b = 42, full knowledge), and which keeps survival
below the cap for all three forms over the whole range of L (largest net fitness about 7.7 - 3.5 = 4.2, survival 0.985).
No other value will be tried before these results are reported.

Design. Simulation 1 protocol unchanged (10,000-generation burn-in, care introduced, end at generation 100,000),
N = 1,000. Seeds: the published seeds of replicates 1-20 (control and treatment) from Data_S1_Summary.csv.
  Arm A (offset 3.5): 20 control runs; 20 treatment runs for each of superlinear, linear, saturating.
  Arm B (offset 0, same code and seeds): 20 treatment runs for linear and for saturating. Control and superlinear
  treatment at offset 0 are the published Data_S1 rows for replicates 1-20 (identical seeds and code).
Outcomes, as in the paper: population means at the final generation; delta L % and delta B % (treatment vs
control); two-sided t-test on L; TOST on B with +/-2% bounds at alpha = 0.05; decoupling criterion
delta L >= 10% and |delta B| < 2%. Also: mean survival and share of agents at the survival cap in the last generation.

Predictions (from the equations, before running).
  P1. Offset 0 arm reproduces the published pattern: B declines under linear and saturating returns.
  P2. Offset 3.5: B does not differ between treatment and control under any of the three forms (|delta B| < 2%).
  P3. Offset 3.5: childhood (L) behaves differently by form: saturating returns give a stable, moderate extension
      (interior optimum near L = 7.25); linear returns push L toward its upper bound (20); superlinear returns
      are bistable around L of about 5.8 (below it no extension, above it L runs toward the bound).
If P2 fails, the published reading (shape of the return matters for brain size) is supported and my concern is wrong.

## Amendment (2026-11-16 21:50 UTC), made after 20 runs of the first attempt and before any treatment comparison was read as a result

Problem found. With the offset applied from generation 0, survival at the arbitrary starting values is about 1%,
and 2 of the first 4 control populations collapsed during the burn-in to a degenerate state (B < 1, L about 0.1).
This is a burn-in artifact of the manipulation, unrelated to the question. Those 20 runs are kept in
attempt1_offset_from_generation0.jsonl and are not used.

Change. The offset is switched on at the start of Phase 2 (generation 10,001), in control and treatment alike.
Phase 1 is therefore exactly the published burn-in. Offset value (3.5), seeds, replicate count, outcomes, and
predictions are unchanged. All 120 runs are redone with the amended code (patch.diff), including the offset-0 arm.
