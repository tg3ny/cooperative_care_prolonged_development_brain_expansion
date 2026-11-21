"""Build a patched copy of the deposited simulation (release v1.0.0) with two switches added.
With the defaults (baseline_offset = 0.0, social_form = 'superlinear') the model is unchanged."""
import sys, hashlib
src_path, out_path = sys.argv[1], sys.argv[2]
src = open(src_path, encoding='utf-8').read()
print('source sha256', hashlib.sha256(src.encode()).hexdigest())
edits = [
 # 1. three new parameters (defaults leave the model unchanged)
 ("    L_max: float = 20.0              # Maximum childhood duration\n",
  "    L_max: float = 20.0              # Maximum childhood duration\n"
  "    baseline_offset: float = 0.0     # CHECK: constant subtracted from net fitness (baseline mortality, Phase 2 on)\n"
  "    social_form: str = 'superlinear' # CHECK: 'superlinear' | 'linear' | 'saturating'\n"
  "    saturation_K: float = 5.0        # CHECK: half-saturation constant, as in alternative_models_analysis.py\n"),
 # 2. alternative benefit forms, expressions copied from alternative_models_analysis.py
 ("        social_benefit = (self.Sub_care * p.social_rate * np.power(L_normalized, p.social_exponent) * p.L_max)\n",
  "        if p.social_form == 'superlinear':\n"
  "            social_benefit = (self.Sub_care * p.social_rate * np.power(L_normalized, p.social_exponent) * p.L_max)\n"
  "        elif p.social_form == 'linear':\n"
  "            social_benefit = self.Sub_care * p.social_rate * self.L\n"
  "        elif p.social_form == 'saturating':\n"
  "            social_benefit = self.Sub_care * p.social_rate * p.L_max * self.L / (self.L + p.saturation_K)\n"
  "        else:\n"
  "            raise ValueError(p.social_form)\n"),
 # 3. baseline mortality
 ("        net_fitness = total_benefit - total_cost\n",
  "        net_fitness = total_benefit - total_cost - self.active_offset\n"),
 # 3b. the offset is inactive during the burn-in and switched on, in both groups, when Phase 2 starts
 ("        self.Sub_care = 0.0\n",
  "        self.Sub_care = 0.0\n"
  "        self.active_offset = 0.0         # CHECK: no baseline mortality during the burn-in\n"),
 ("    pop.Sub_care = params.Sub_care_treatment if group == 'treatment' else 0.0\n",
  "    pop.Sub_care = params.Sub_care_treatment if group == 'treatment' else 0.0\n"
  "    pop.active_offset = params.baseline_offset   # CHECK: applies from Phase 2 on, in control and treatment alike\n"),
 # 4. two read-only diagnostics returned per replicate (no random numbers consumed)
 ("        'final_K': stats['K_mean'],\n    }\n",
  "        'final_K': stats['K_mean'],\n"
  "        'final_surv_mean': float(np.mean(pop.prev_fitness)),\n"
  "        'final_frac_at_cap': float(np.mean(pop.prev_fitness >= 0.999)),\n"
  "        'final_B_sd': stats['B_std'],\n"
  "        'final_L_sd': stats['L_std'],\n    }\n"),
]
for old, new in edits:
    assert src.count(old) == 1, ('expected exactly one match', old[:60], src.count(old))
    src = src.replace(old, new)
open(out_path, 'w', encoding='utf-8').write(src)
print('wrote', out_path)
