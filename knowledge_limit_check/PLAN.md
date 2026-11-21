# Knowledge-limit analysis: plan written before the runs (2026-11-20 20:37 UTC)

Question (Reviewer 1). Knowledge in the model is capped at 0.6 x brain size (K <= 0.6B). Does the main result of
Simulation 1, that cooperative care lengthens childhood (L) without changing brain size (B), hold when this limit is
tightened (0.4B) or relaxed (0.8B)?

Disclosure. A pilot with 5 replicates per condition was run on 2026-11-16 (saturating return only; limits 0.4B and
0.8B; constants 0 and 7.8). Its results are known to me: control L 3.52 and 6.39, care L 7.64 and 8.16, brain size
under care +1.0% and -0.5% relative to controls. The superlinear (default) return has not been run at either limit.
The pilot files are kept in pilot_2026-11-16/ and are not used in the analysis.

Model. selection_strength_check/sim_patched.py with one further edit (make_patched_cap.py): the coefficient 0.6 in the
storage rule becomes a parameter, knowledge_cap (default 0.6). With the default the code is unchanged; this is checked
by rerunning 4 runs of the selection-strength analysis and comparing them with the stored results.

Design. Simulation 1 protocol (10,000-generation burn-in, care introduced, end at generation 100,000), N = 1,000,
recorded seeds of replicates 1-20 of Simulation 1. For each limit (0.4B, 0.8B): 20 control runs, 20 care runs with the
superlinear (default) return, 20 care runs with the saturating return. 120 runs. The limit applies throughout the run,
including the burn-in. The limit 0.6B is covered by the selection-strength analysis (same seeds).

Selection strength. As in the selection-strength analysis, a constant is subtracted from net fitness in both
conditions from the start of Phase 2 so that baseline mortality in control populations is about 10% per generation.
The constant is fixed in advance from the equations: net fitness at the predicted control equilibrium minus ln(9)
(the net fitness at which survival is 0.90), rounded to 0.1. Predicted control equilibrium: B* where the marginal
return on brain size equals the marginal penalty near the ceiling, and L* where asocial learning fills capacity
(a L^q = cap). This gives B* = 39.1, L* = 3.3, constant -0.3 at 0.4B (that is, 0.3 is added), and B* = 43.9, L* = 6.1,
constant 7.9 at 0.8B. (The same formula gives 3.7 at 0.6B; the selection-strength analysis used 3.5, derived from the
observed Simulation 1 control endpoint, and realized control mortality of 9%.) No other values will be tried.

Outcomes, as in the paper: population means at the final generation; delta L % and delta B % (care vs control);
two-sided t-test on L; TOST on B with +/-2% bounds at alpha = 0.05; decoupling criterion (delta L >= 10% and
|delta B| < 2%); mean survival, share of agents at the survival bound, and K relative to its limit.

Predictions (from the equations and, where stated, the pilot).
  P1. Brain size. At both limits and under both forms, B under care is equivalent to B in controls (|delta B| < 2%):
      once knowledge is at its limit, a longer childhood adds no knowledge and does not change selection on B.
  P2. Baseline childhood. Control L settles close to the duration at which learning fills capacity: about 3.5 at 0.4B
      and about 6.4 at 0.8B (pilot), against 5.0 at 0.6B.
  P3. Saturating return. L under care moves to the interior optimum near 7.25 at both limits (pilot: 7.6 and 8.2).
  P4. Superlinear return (not piloted). The marginal return exceeds the temporal cost only above L = 5.78.
      At 0.8B the baseline (about 6.4) is above this threshold, so L should run toward its upper limit (20).
      At 0.4B the baseline (about 3.5) is far below it, so little or no lengthening of childhood is expected.
If P4 holds at 0.4B, the manuscript must say that under an accelerating return the lengthening of childhood depends on
the baseline duration, which the knowledge limit sets. If P1 fails, the claim that the brain-size result does not
depend on the value of the limit is wrong and must be withdrawn.
