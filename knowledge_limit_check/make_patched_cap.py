"""Add one parameter to selection_strength_check/sim_patched.py: the coefficient of the storage rule K <= 0.6 B.
With the default (knowledge_cap = 0.6) the model is unchanged.
Usage: python knowledge_limit_check/make_patched_cap.py selection_strength_check/sim_patched.py knowledge_limit_check/sim_patched_cap.py"""
import sys
src_path, out_path = sys.argv[1], sys.argv[2]
src = open(src_path, encoding='utf-8').read()
edits = [
 ("    saturation_K: float = 5.0        # CHECK: half-saturation constant, as in alternative_models_analysis.py\n",
  "    saturation_K: float = 5.0        # CHECK: half-saturation constant, as in alternative_models_analysis.py\n"
  "    knowledge_cap: float = 0.6       # CHECK: coefficient of the storage rule K <= knowledge_cap * B\n"),
 ("        self.K = np.minimum(K_total, self.B * 0.6)\n",
  "        self.K = np.minimum(K_total, self.B * p.knowledge_cap)\n"),
]
for old, new in edits:
    assert src.count(old) == 1, ('expected exactly one match', old[:60], src.count(old))
    src = src.replace(old, new)
open(out_path, 'w', encoding='utf-8').write(src)
print('wrote', out_path)
