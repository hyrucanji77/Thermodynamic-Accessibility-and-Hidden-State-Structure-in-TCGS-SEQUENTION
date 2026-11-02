#!/usr/bin/env python3
"""No-branch-fit transfer test of the Figure 4 source-data column labels.

Run after ``python code/reproduce.py``. Original workbooks are read, never edited.
All 24 column-to-branch permutations are compared by two descriptive scores:
absolute-work curve MSE and squared log-ratios of inferred normalization constants.
Neither score is a likelihood or a test based on independent time samples.
Model parameters, threshold, delay, and source speeds are not optimized here.
"""
from __future__ import annotations
import csv
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import zipfile
from xml.etree import ElementTree as ET
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from accessibility import Parameters, binary_condition, constant_speed_work, accessible, riccati

ROOT = Path(__file__).resolve().parents[1]
OUT, FIG = ROOT / "results", ROOT / "figures"
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

def read_workbook(path: Path) -> dict[str, list[list]]:
    """Read cached OOXML values without modifying spreadsheet bytes."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            strings = ["".join(e.itertext()) for e in ET.fromstring(z.read("xl/sharedStrings.xml"))]
        rels = {r.attrib["Id"]: r.attrib["Target"] for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        book = ET.fromstring(z.read("xl/workbook.xml"))
        result = {}
        for sheet in book.find("m:sheets", NS):
            target = rels[sheet.attrib[RNS + "id"]]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            rows = []
            for row in ET.fromstring(z.read(target)).findall("m:sheetData/m:row", NS):
                values = {}
                for cell in row:
                    letters = "".join(c for c in cell.attrib["r"] if c.isalpha())
                    col = 0
                    for c in letters: col = 26 * col + ord(c) - 64
                    elem = cell.find("m:v", NS)
                    typ = cell.attrib.get("t")
                    if typ == "inlineStr": value = "".join(cell.find("m:is", NS).itertext())
                    elif elem is None: value = None
                    elif typ == "s": value = strings[int(elem.text)]
                    elif typ in ("str", "e"): value = elem.text
                    else: value = float(elem.text)
                    values[col - 1] = value
                rows.append([values.get(c) for c in range(max(values, default=-1) + 1)])
            result[sheet.attrib["name"]] = rows
        return result

def write_csv(name, header, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f); writer.writerow(header); writer.writerows(rows)

checks = []
def check(name: str, residual: float, tolerance: float) -> None:
    residual = float(residual)
    ok = bool(np.isfinite(residual) and residual <= tolerance)
    checks.append(dict(name=name, residual=residual, tolerance=tolerance, passed=ok))
    if not ok: raise AssertionError(checks[-1])

p = Parameters()
path = ROOT / "data" / "41567_2026_3442_MOESM6_ESM.xlsx"
hash_before = hashlib.sha256(path.read_bytes()).hexdigest()
sheets = read_workbook(path)
columns = (1, 3, 5, 7)
column_names = ("B", "D", "F", "H")
states = ("00", "01", "10", "11")
# Use all 300 paired positive-time samples in the declared 0--3 s protocol.
def data_rows(name):
    return np.array([r[:9] for r in sheets[name][1:]
                     if len(r) >= 9 and isinstance(r[0], (int, float))
                     and 0 < r[0] <= 3], dtype=float)
a, n = data_rows("Inset_4b"), data_rows("Figure_4b")
check("paired_time_grids", np.max(abs(a[:, 0] - n[:, 0])), 1e-14)
check("sample_count_300", abs(len(a) - 300), 0)
times = a[:, 0]
absolute, normalized = a[:, columns], n[:, columns]
implied = np.zeros(4)
for j, col in enumerate(column_names):
    valid = abs(normalized[:, j]) > 1e-8
    ratios = absolute[valid, j] / normalized[valid, j]
    implied[j] = np.median(ratios)
    check(f"constant_normalization_{col}", np.max(abs(ratios / implied[j] - 1)), 2e-13)

branches = []
curves = np.empty((len(times), 4))
for j, state in enumerate(states):
    prob, mu, _ = binary_condition(1, *map(int, state), threshold=70, p=p)
    speed = 25. if state == "10" else 45.
    baseline = .5 * p.k * mu[0] ** 2 / p.kBT
    curves[:, j] = [constant_speed_work(mu, t, speed, p)[0] / p.kBT for t in times]
    # Independent mean-equation and work-integral calculation of the ramp.
    def rhs(t, z):
        x, y, _ = z
        return [(-p.k * (x-speed*t) - p.kb * (x-y)) / p.gamma,
                p.kb * (x-y) / p.gammab, p.k * (x-speed*t) * speed / p.kBT]
    ode = solve_ivp(rhs, (0, 3), [*mu, 0], t_eval=times, rtol=2e-11, atol=2e-11)
    if not ode.success: raise RuntimeError(ode.message)
    check(f"ramp_closed_form_vs_ODE_{state}", np.max(abs(ode.y[2] - curves[:, j])), 2e-9)
    candidates = [(0., 0.), (3., curves[-1, j] / baseline)]
    opt = minimize_scalar(lambda t: -constant_speed_work(mu, t, speed, p)[0]/p.kBT/baseline,
                          bounds=(0., 3.), method="bounded", options={"xatol": 1e-12})
    candidates.append((float(opt.x), float(-opt.fun)))
    pt, pv = max(candidates, key=lambda q: q[1])
    jump = float(mu[0] + (p.kb-.5*riccati(3,p)) * (mu[0]-mu[1]) / p.k)
    ideal = float(accessible(mu, 3, p) / p.kBT)
    check(f"ramp_value_below_ideal_{state}", max(0., curves[-1,j]-ideal), 2e-12)
    branches.append(dict(state=state, probability=prob, mean_x_nm=float(mu[0]),
                         mean_y_nm=float(mu[1]), B_kBT=baseline, speed_nm_s=speed,
                         W_3s_kBT=float(curves[-1,j]), W_ideal_3s_kBT=ideal,
                         model_peak_time_s=pt, model_peak_ratio=pv,
                         initial_jump_nm=jump))

baselines = np.array([b["B_kBT"] for b in branches])
assignments = []
for perm in itertools.permutations(range(4)):
    work_score = float(np.mean((absolute-curves[:,perm])**2))
    denominator_score = float(np.sum(np.log(implied/baselines[list(perm)])**2))
    assignments.append(dict(mapping=[states[j] for j in perm],
                            absolute_work_mse=work_score, log_denominator_score=denominator_score))
work_order = sorted(assignments, key=lambda r: r["absolute_work_mse"])
b_order = sorted(assignments, key=lambda r: r["log_denominator_score"])
best = work_order[0]["mapping"]
check("two_scores_select_same_permutation", float(best != b_order[0]["mapping"]), 0)
check("selected_permutation_regression", float(best != ["11", "00", "01", "10"]), 0)
check("work_score_unique_minimum", float(work_order[1]["absolute_work_mse"] <= work_order[0]["absolute_work_mse"]), 0)
check("denominator_score_unique_minimum", float(b_order[1]["log_denominator_score"] <= b_order[0]["log_denominator_score"]), 0)
check("branch_10_jump_negative", max(0., branches[2]["initial_jump_nm"]), 0)
check("branch_10_initial_jump_value", abs(branches[2]["initial_jump_nm"] + 24.560441263094), 1e-9)
# Positivity invoked in Proposition 9.1, checked over a logarithmic test grid.
check("Metzler_offdiagonal_positivity", max(0., p.A[0,1], p.A[1,0]), 0)
vals = [expm(-p.A*t)[0,1] for t in np.geomspace(1e-7,100,100)]
check("positive_semigroup_12", max(0., -min(vals)), 0)
check("input_file_unchanged", float(hashlib.sha256(path.read_bytes()).hexdigest() != hash_before), 0)

rows = []
for c, col in enumerate(column_names):
    state = best[c]; branch = branches[states.index(state)]
    peak = int(np.argmax(normalized[:,c]))
    rows.append(dict(column=col, stored_state=states[c], inferred_state=state,
                     inferred_B_kBT=float(implied[c]), model_B_kBT=branch["B_kBT"],
                     relative_B_residual=float(implied[c]/branch["B_kBT"]-1),
                     implied_abs_mean_x_nm=float(np.sqrt(2*p.kBT*implied[c]/p.k)),
                     model_mean_x_nm=branch["mean_x_nm"],
                     model_abs_mean_x_nm=abs(branch["mean_x_nm"]),
                     abs_mean_offset_nm=float(np.sqrt(2*p.kBT*implied[c]/p.k)-abs(branch["mean_x_nm"])),
                     denominator_score_component=float(np.log(implied[c]/branch["B_kBT"])**2),
                     workbook_peak_time_s=float(times[peak]), workbook_peak_ratio=float(normalized[peak,c]),
                     model_peak_time_s=branch["model_peak_time_s"], model_peak_ratio=branch["model_peak_ratio"],
                     absolute_curve_RMSE_kBT=float(np.sqrt(np.mean((absolute[:,c]-curves[:,states.index(state)])**2)))))
# Baseline inversion identifies a magnitude, not the sign of a conditional mean.
# No source values or fitted parameters are changed by this re-expression.
for row in rows:
    mu_abs = row["implied_abs_mean_x_nm"]
    model_abs = row["model_abs_mean_x_nm"]
    check(f"mean_inversion_round_trip_{row['column']}",
          abs(.5*p.k*mu_abs**2/p.kBT-row["inferred_B_kBT"]), 2e-14)
    check(f"baseline_vs_mean_residual_identity_{row['column']}",
          abs((mu_abs/model_abs)**2-1-row["relative_B_residual"]), 2e-14)
components = [row["denominator_score_component"] for row in rows]
check("denominator_score_component_sum", abs(sum(components)-b_order[0]["log_denominator_score"]), 2e-14)
swapped = [best[2],best[1],best[0],best[3]]
swap_record = next(r for r in assignments if r["mapping"] == swapped)
swap_terms = np.log(implied/np.array([baselines[states.index(m)] for m in swapped]))**2
check("BF_swap_preserves_D_H_components", max(abs(swap_terms[j]-components[j]) for j in (1,3)), 2e-14)
check("BF_swap_accounts_for_full_denominator_margin",
      abs(sum(swap_terms[j]-components[j] for j in (0,2))
          -(b_order[1]["log_denominator_score"]-b_order[0]["log_denominator_score"])), 2e-14)
check("BF_swap_is_denominator_runner_up", float(swapped != b_order[1]["mapping"]), 0)
score_diagnostics = dict(
    column_order=list(column_names),
    denominator_components=components,
    denominator_D_fraction=float(components[1]/sum(components)),
    denominator_runner_up_gap=float(b_order[1]["log_denominator_score"]-b_order[0]["log_denominator_score"]),
    denominator_runner_up_gap_relative_to_minimum=float(b_order[1]["log_denominator_score"]/b_order[0]["log_denominator_score"]-1),
    BF_swapped_mapping=swapped,
    BF_swapped_denominator_components=swap_terms.tolist(),
    BF_swapped_work_score=swap_record["absolute_work_mse"],
    BF_swap_work_gap=float(swap_record["absolute_work_mse"]-work_order[0]["absolute_work_mse"]),
    min_denominator_score_with_wrong_D=min(r["log_denominator_score"] for r in assignments if r["mapping"][1] != best[1]),
    min_denominator_score_with_wrong_H=min(r["log_denominator_score"] for r in assignments if r["mapping"][3] != best[3]),
    interpretation="Descriptive separation conditional on fixed parameters; not two statistically independent confirmations."
)
write_csv("label_resolution.csv", list(rows[0]), [[r[k] for k in rows[0]] for r in rows])
write_csv("label_score_components.csv", ["column","selected_branch","selected_denominator_component","BF_swapped_branch","BF_swapped_denominator_component"],
          [[col,best[j],components[j],swapped[j],float(swap_terms[j])] for j,col in enumerate(column_names)])
write_csv("label_assignment_scores.csv", ["rank_by_work", "B_branch", "D_branch", "F_branch", "H_branch", "absolute_work_mse", "log_denominator_score"],
          [[k+1,*r["mapping"],r["absolute_work_mse"],r["log_denominator_score"]] for k,r in enumerate(work_order)])
write_csv("constant_speed_trace_comparison.csv", ["time_s"]+[f"source_{c}_kBT" for c in column_names]+[f"model_{s}_kBT" for s in states],
          [[float(t),*absolute[j].tolist(),*curves[j].tolist()] for j,t in enumerate(times)])
summary = dict(method="Retrospective no-branch-fit parameter transfer; two descriptive scores, no statistical independence assumed",
               workbook_sha256=hash_before, sample_count=len(times), window_s=[0,3],
               branches=branches, columns=rows, selected_mapping=best,
               best_work_score=work_order[0]["absolute_work_mse"],
               runner_up_work_score=work_order[1]["absolute_work_mse"],
               runner_up_work_mapping=work_order[1]["mapping"],
               best_denominator_score=b_order[0]["log_denominator_score"],
               runner_up_denominator_score=b_order[1]["log_denominator_score"],
               runner_up_denominator_mapping=b_order[1]["mapping"],
               score_diagnostics=score_diagnostics,
               mean_inversion_assumptions="Transferred trap stiffness and temperature; baseline inversion determines only |mean_x|, not its sign.",
               new_fitted_parameters=0,
               limitations=["Source parameters were calibrated upstream, not parameter-free physics",
                            "Retrospective assignment is not a prospective or independent replication",
                            "Column D normalization differs by approximately 94%, corresponding to a +4.8032 nm mean-magnitude offset under transferred stiffness",
                            "The two scores reuse the same data and model; their minima are not independent statistical confirmations",
                            "Bandwidth has not been modeled as a cause of the residuals",
                            "Pointwise error bars do not determine uncertainty of curve maxima"])
(OUT/"label_resolution.json").write_text(json.dumps(summary,indent=2), encoding="utf-8")
(OUT/"label_resolution_verification.json").write_text(json.dumps(dict(
    python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
    checks=checks,n_checks=len(checks),all_passed=all(c["passed"] for c in checks)),indent=2))
# Absolute-work comparison: one panel, default color cycle, no source redrawing.
fig, ax = plt.subplots(figsize=(6.4,4.0))
for c, col in enumerate(column_names):
    state = best[c]; idx = states.index(state)
    ax.plot(times, curves[:,idx], label=f"Model {state}")
# The cycle repeats in the same branch order for sampled source markers.
ax.set_prop_cycle(None)
for c, col in enumerate(column_names):
    ax.plot(times[::12],absolute[::12,c],linestyle="none",marker="o",markersize=3.1,
            fillstyle="none",label=f"Workbook {col} -> {best[c]}")
ax.axhline(0, linewidth=.6, linestyle=":")
ax.set_xlabel("Extraction time (s)"); ax.set_ylabel(r"Extracted work / $k_{\mathrm{B}}T$")
ax.legend(ncol=2,fontsize=8); ax.grid(alpha=.15)
fig.tight_layout(); fig.savefig(FIG/"label_resolution.pdf"); plt.close(fig)
print(json.dumps({"new_checks":len(checks),"all_passed":True,"mapping":best,
                  "best_work_score":work_order[0]["absolute_work_mse"],
                  "jump_10_nm":branches[2]["initial_jump_nm"]},indent=2))
