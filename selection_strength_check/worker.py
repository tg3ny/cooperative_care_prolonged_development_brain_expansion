"""Run a list of replicate tasks with the patched simulation and append one JSON line per finished task."""
import sys, json, time, importlib.util, os
from dataclasses import replace
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('sim_patched', os.path.join(here, 'sim_patched.py'))
sim = importlib.util.module_from_spec(spec); sys.modules['sim_patched'] = sim; spec.loader.exec_module(sim)

task_file, out_file, worker, n_workers = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
tasks = [json.loads(l) for l in open(task_file)]
done = set()
if os.path.exists(out_file):
    for l in open(out_file):
        try: done.add(json.loads(l)['task_id'])
        except Exception: pass
base = sim.SimulationParameters()          # deposited defaults: T_epoch1 = 10,000; T_epoch2 = 100,000
for i, t in enumerate(tasks):
    if i % n_workers != worker or t['task_id'] in done: continue
    p = replace(base, baseline_offset=t['offset'], social_form=t['form'])
    t0 = time.time()
    r = sim.run_replicate(t['rep'], t['group'], t['seed'], p)   # the deposited replicate runner, unmodified loop
    h = r.pop('history')
    rec = dict(t, **{k: r[k] for k in ('final_B', 'final_L', 'final_S', 'final_K', 'final_surv_mean', 'final_frac_at_cap', 'final_B_sd', 'final_L_sd')})
    rec['phase1_end'] = {k: h[99][k] for k in ('generation', 'B_mean', 'L_mean', 'S_mean', 'K_mean')}   # sample at generation 9,900
    rec['traj'] = [[x['generation'], round(x['B_mean'], 4), round(x['L_mean'], 4), round(x['S_mean'], 4), round(x['K_mean'], 4)] for x in h[::10]]
    rec['seconds'] = round(time.time() - t0, 1)
    with open(out_file, 'a') as f: f.write(json.dumps(rec) + '\n')
