#!/bin/bash
# Queued Z12 pin-scan: tests whether lowering the edge-degree pin below R's
# native q-bar (5.10423) unlocks Z12 solubility (curvature-lowering dopant).
# Waits for the machine to free up, launches the 5-point pin ladder at mu=4,
# then analyzes n_excess(pin), q-bar(pin), legality(pin).
#   setsid nohup bash scripts/r_pinscan_queue.sh &
cd /home/aaron/Desktop/Discrete-Differential-Geometry
SP=/tmp/claude-1000/-home-aaron-Desktop-Discrete-Differential-Geometry/93453a5d-72d0-489e-b664-2065eda88aa4/scratchpad
exec > "$SP/r_pinscan.log" 2>&1
echo "### pin-scan queued $(date) ###"

# wait until the machine has real headroom (current Z12 isotherm + large tier drained)
while [ "$(python3 -c 'import sys;sys.path.insert(0,"tools");from seed_utils import get_free_memory_gb;print(int(get_free_memory_gb()))')" -lt 14 ]; do
  sleep 120
done
echo "headroom available; launching pin-scan $(date)"

python3 scripts/r_solubility.py --species Z12 \
  --edge-target 5.1042 5.100 5.097 5.094 5.091 --mu 4 --replicas 3 \
  --sweeps 20000 --out-dir data/r_pinscan --base-seed 33000 \
  --mem-floor-gb 9 --stagger 12

# wait for all 15 finals
while [ "$(ls data/r_pinscan/r_sol_Z12_*_final.mfd 2>/dev/null | wc -l)" -lt 15 ]; do
  sleep 60
done
echo "pin-scan complete $(date); analyzing"

python3 - <<'PY'
import csv, glob, re
import numpy as np
from collections import defaultdict
NATIVE_Z12 = 2187
g = defaultdict(list)
for p in glob.glob("data/r_pinscan/r_sol_Z12_et*_mu4_r*.csv"):
    m = re.search(r"_et([0-9p]+)_mu", p)
    et = float(m.group(1).replace("p", "."))
    g[et].append(p)
def equil(rows, col):
    v=[float(r[col]) for r in rows if float(r["sweeps"])>=12000]
    return np.mean(v) if v else float("nan")
print(f"{'edge_pin':>9} {'d_from_native':>13} {'Z12_excess':>11} {'qbar':>9} {'pure56':>7}")
rows_out=[]
for et in sorted(g, reverse=True):
    exc,qb,p56=[],[],[]
    for p in g[et]:
        r=[x for x in csv.DictReader(o for o in open(p) if not o.startswith("#"))]
        if not r: continue
        exc.append(equil(r,"n_dop")-NATIVE_Z12); qb.append(equil(r,"mean_edeg")); p56.append(equil(r,"pure56"))
    print(f"{et:9.4f} {et-5.10423:>+13.5f} {np.mean(exc):>11.1f} {np.mean(qb):>9.5f} {np.mean(p56):>7.3f}")
    rows_out.append((et,np.mean(exc),np.mean(qb),np.mean(p56)))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
et=[r[0] for r in rows_out]; exc=[r[1] for r in rows_out]; qb=[r[2] for r in rows_out]; p56=[r[3] for r in rows_out]
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].plot(et,exc,"o-",color="#762a83"); ax[0].axvline(5.10423,color="k",ls="--",lw=.8)
ax[0].set_xlabel("edge-degree pin target"); ax[0].set_ylabel("excess Z12 (uptake)")
ax[0].set_title("Does lowering the pin unlock Z12?"); ax[0].invert_xaxis(); ax[0].grid(alpha=.3)
ax[1].plot(et,qb,"o-",color="#2166ac",label="achieved q̄"); ax[1].plot(et,et,":",color="gray",label="pin target")
ax[1].axhline(5.10423,color="green",ls="--",lw=.8,label="native")
ax[1].set_xlabel("edge-degree pin target"); ax[1].set_ylabel("achieved q̄"); ax[1].invert_xaxis()
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_title("Below-flat R via Z12?")
fig.suptitle("Z12 solubility vs edge-pin displacement (mu=4)"); fig.tight_layout()
fig.savefig("out/r_pinscan.png",dpi=125); print("wrote out/r_pinscan.png")
PY
echo "### pin-scan analysis done $(date) ###"
