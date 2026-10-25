# -*- coding: utf-8 -*-
"""多模型候选：找出能复现 EC50_summary.csv 的那一个"""
import os, re, warnings, numpy as np, pandas as pd
from scipy.optimize import curve_fit
warnings.filterwarnings("ignore")
BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUT  = os.path.join(BASE, "教学示范图")
L = []

def norm(s):
    if s is None: return ""
    s = re.sub(r"\s+", "", str(s))
    for a, b in [("（", "("), ("）", ")")]: s = s.replace(a, b)
    return s
ALIAS = {"二甲基三硫醚": "二甲基三硫",
         "反，顺-2，6-壬二稀醇": "反，顺-2,6-壬二烯醛",
         "反，顺-2,6-壬二稀醇": "反，顺-2,6-壬二烯醛"}
def base_name(s):
    b = re.sub(r"\(C\d+\)$", "", norm(s), flags=re.I)
    return ALIAS.get(b, b)
def level(s):
    m = re.search(r"\(C(\d+)\)$", norm(s), flags=re.I)
    return int(m.group(1)) if m else None

sh5 = pd.read_excel(os.path.join(BASE, "单体混合感官特征理化数据.xlsx"), sheet_name="Sheet5")
sh5.columns = ["key", "C"]
sh5["key"] = sh5["key"].map(norm); sh5["C"] = pd.to_numeric(sh5["C"], errors="coerce")
sh5["cmpd"] = sh5["key"].map(base_name); sh5["lvl"] = sh5["key"].map(level)

tab = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"), sheet_name="单体数据平均值")
tab["key"] = tab["Sample"].map(norm).map(
    lambda k: "%s(C%d)" % (ALIAS.get(base_name(k), base_name(k)), level(k))
    if level(k) is not None else k)
tab["y"] = pd.to_numeric(tab["气味强度"], errors="coerce")
Y = tab.set_index("key")["y"].to_dict()

raw = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"), sheet_name="单体数据原始数据")
raw = raw.rename(columns={raw.columns[2]: "I"})
raw["key"] = raw["Sample"].ffill().map(norm)
raw["I"] = pd.to_numeric(raw["I"], errors="coerce")
raw["Subject"] = raw["Subject"].astype(str).str.strip()
raw = raw[raw["Subject"].str.match(r"^S\d+$", na=False) & raw["I"].notna()]
raw["occ"] = raw.groupby(["key", "Subject"]).cumcount()
raw = raw[raw["occ"] == 0].copy()
raw["cmpd"] = raw["key"].map(base_name)

ref = pd.read_csv(os.path.join(OUT, "EC50_summary.csv"), encoding="utf-8-sig")
ref["Compound"] = ref["Compound"].map(lambda s: ALIAS.get(s, s))
R = dict(zip(ref["Compound"], ref.itertuples(index=False)))

def hill3(C, Imax, n, EC50):
    C = np.asarray(C, float); return Imax * C ** n / (EC50 ** n + C ** n)
def hill4(C, d, Imax, n, EC50):
    C = np.asarray(C, float); return d + (Imax - d) * C ** n / (EC50 ** n + C ** n)

def fit_models(x, y):
    out = {}
    for tag, fn, p0s, bnd in [
        ("3PL", hill3,
         [[max(y.max() * 1.1, 1), 1, np.median(x)],
          [max(y.max(), 1), 0.3, x.min() * 10],
          [max(y.max(), 1), 3, x.max() / 10]],
         ([1e-3, 0.05, x.min() / 1e3], [60, 20, x.max() * 1e3])),
        ("4PL", hill4,
         [[y.min(), max(y.max() * 1.1, y.min() + 1), 1, np.median(x)],
          [max(y.min() - 1, -20), y.max() + 2, 0.5, np.median(x)],
          [y.min(), y.max() + 1, 3, np.median(x)]],
         ([-30, -30, 0.05, x.min() / 1e3], [60, 60, 20, x.max() * 1e3])),
    ]:
        best = None
        for p0 in p0s:
            try:
                p, _ = curve_fit(fn, x, y, p0=p0, bounds=bnd, maxfev=300000)
                r = y - fn(x, *p)
                ss = 1 - np.sum(r ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12)
                if best is None or ss > best[1]: best = (p, ss)
            except Exception:
                pass
        if best: out[tag] = best
    return out

def compare(tag, rows):
    A = pd.DataFrame(rows)
    j = A.dropna(subset=["Imax_ref"])
    L.append("\n--- %s vs EC50_summary.csv (%d 可比) ---" % (tag, len(j)))
    for col, rc in [("Imax", "Imax_ref"), ("n", "n_ref"), ("EC50", "EC50_ref")]:
        r = (j[col] / j[rc]).replace([np.inf, -np.inf], np.nan).dropna()
        L.append("   %-5s 中位倍数=%.3f 几何均值=%.3f 范围 %.3f–%.3f |误差|<10%%: %d/%d"
                 % (col, r.median(), np.exp(np.log(r.clip(1e-12)).mean()),
                    r.min(), r.max(), int(((r > 0.9) & (r < 1.1)).sum()), len(r)))
    L.append("   R²: 我=%.4f±%.4f  ref=%.4f±%.4f"
             % (j["R2"].mean(), j["R2"].std(), j["R2_ref"].mean(), j["R2_ref"].std()))
    return A

L.append("模型候选复现测试（数据 = 论文《单体数据平均值》+ Sheet5 浓度）")

for model in ["3PL", "4PL"]:
    rows = []
    for cmpd, d in sh5.groupby("cmpd"):
        d = d.dropna(subset=["C"]).sort_values("lvl")
        y = np.array([Y.get(k, np.nan) for k in d["key"]], float)
        m = np.isfinite(y)
        if m.sum() < 4: continue
        x, yy = d["C"].to_numpy(float)[m], y[m]
        res = fit_models(x, yy)
        if model not in res: continue
        p, r2 = res[model]
        rf = R.get(cmpd)
        rows.append(dict(cmpd=cmpd, n_pts=len(x), Imax=p[1] if model == "4PL" else p[0],
                         n=p[2] if model == "4PL" else p[1],
                         EC50=p[3] if model == "4PL" else p[2], R2=r2,
                         Imax_ref=getattr(rf, "Imax", np.nan),
                         n_ref=getattr(rf, "Hill_n", np.nan),
                         EC50_ref=getattr(rf, "EC50_ppb", np.nan),
                         R2_ref=getattr(rf, "R2", np.nan)))
    compare(model, rows)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "refit_%s.csv" % model),
                              index=False, encoding="utf-8-sig")

# 附：3PL/4PL 逐行明细 (4PL)
L.append("\n" + "=" * 90)
L.append("4PL 逐行明细")
L.append("=" * 90)
A = pd.read_csv(os.path.join(OUT, "refit_4PL.csv"), encoding="utf-8-sig")
L.append("%-16s %3s | %8s %8s %12s %7s | %8s %8s %12s %7s"
         % ("化合物", "n", "Imax", "n", "EC50", "R²", "Imax_ref", "n_ref", "EC50_ref", "R2_ref"))
for _, r in A.iterrows():
    L.append("%-16s %3d | %8.3f %8.3f %12.4g %7.4f | %8.3f %8.3f %12.4g %7.4f"
             % (r["cmpd"], r["n_pts"], r["Imax"], r["n"], r["EC50"], r["R2"],
                r["Imax_ref"], r["n_ref"], r["EC50_ref"], r["R2_ref"]))

open(r"C:\Users\24252\AppData\Local\Temp\opencode\refit2.txt", "w", encoding="utf-8").write("\n".join(L))
print("ok")
