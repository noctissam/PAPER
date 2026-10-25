# -*- coding: utf-8 -*-
"""独立复现 3 参数 Hill 拟合，与 EC50_summary.csv 逐行比对"""
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

# ---- 浓度 ----
sh5 = pd.read_excel(os.path.join(BASE, "单体混合感官特征理化数据.xlsx"), sheet_name="Sheet5")
sh5.columns = ["key", "C"]
sh5["key"] = sh5["key"].map(norm)
sh5["cmpd"] = sh5["key"].map(base_name)
sh5["lvl"] = sh5["key"].map(level)
sh5["C"] = pd.to_numeric(sh5["C"], errors="coerce")
L.append("Sheet5: %d 行, 化合物 %d, C 范围 %.4g ~ %.4g ppb"
         % (len(sh5), sh5["cmpd"].nunique(), sh5["C"].min(), sh5["C"].max()))
ser = sh5.groupby("cmpd").apply(lambda d: pd.Series({
    "n_lvl": d["lvl"].nunique(), "Cmin": d["C"].min(), "Cmax": d["C"].max(),
    "fold": d["C"].max() / d["C"].min()}), include_groups=False)
L.append("\n各化合物浓度阶梯:")
L.append(ser.round(4).to_string())

# ---- 响应 ----
def hill(C, Imax, n, EC50):
    C = np.asarray(C, float)
    return Imax * C ** n / (EC50 ** n + C ** n)

def fit(x, y):
    best = None
    for n0 in [0.3, 1.0, 3.0]:
        for e0 in [np.sqrt(x.min() * x.max()), x.min() * 10, x.max() / 10]:
            try:
                p, _ = curve_fit(hill, x, y,
                                 p0=[max(y.max() * 1.15, 1e-3), n0, e0],
                                 bounds=([1e-3, 0.05, x.min() / 1e3],
                                         [60.0, 20.0, x.max() * 1e3]),
                                 maxfev=200000)
                resid = y - hill(x, *p)
                ss = 1 - np.sum(resid ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12)
                if best is None or ss > best[1]: best = (p, ss)
            except Exception:
                pass
    return best

tab = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"), sheet_name="单体数据平均值")
tab["key"] = tab["Sample"].map(norm).map(
    lambda k: "%s(C%d)" % (ALIAS.get(base_name(k), base_name(k)), level(k))
    if level(k) is not None else k)
tab["y"] = pd.to_numeric(tab["气味强度"], errors="coerce")

clean = pd.read_csv(os.path.join(OUT, "clean_assessor_level.csv"), encoding="utf-8-sig")
my = clean.groupby("key")["I"].mean()

ref = pd.read_csv(os.path.join(OUT, "EC50_summary.csv"), encoding="utf-8-sig")
ref["Compound"] = ref["Compound"].map(lambda s: ALIAS.get(s, s))
R = dict(zip(ref["Compound"], ref.itertuples(index=False)))

def run(src, label):
    L.append("\n" + "=" * 88)
    L.append("【%s】3 参数 Hill 拟合  I(C) = Imax·C^n/(EC50^n + C^n)" % label)
    L.append("=" * 88)
    L.append("%-16s %3s %8s %8s %12s %8s | %8s %8s %12s"
             % ("化合物", "n", "Imax", "n(Hill)", "EC50_ppb", "R²",
                "Imax_ref", "n_ref", "EC50_ref"))
    rows = []
    for cmpd, d in sh5.groupby("cmpd"):
        d = d.dropna(subset=["C"]).sort_values("lvl")
        y = np.array([src.get(k, np.nan) for k in d["key"]], float)
        m = np.isfinite(y)
        if m.sum() < 4: continue
        x, y = d["C"].to_numpy(float)[m], y[m]
        r = fit(x, y)
        if r is None: continue
        (Imax, n, EC50), r2 = r
        ref_ = R.get(cmpd)
        rows.append(dict(cmpd=cmpd, n_pts=len(x), Imax=Imax, Hill_n=n, EC50=EC50, R2=r2,
                         Imax_ref=getattr(ref_, "Imax", np.nan),
                         n_ref=getattr(ref_, "Hill_n", np.nan),
                         EC50_ref=getattr(ref_, "EC50_ppb", np.nan)))
        L.append("%-16s %3d %8.3f %8.3f %12.1f %8.4f | %8s %8s %12s"
                 % (cmpd, len(x), Imax, n, EC50, r2,
                    ("%8.3f" % ref_.Imax) if ref_ else "       -",
                    ("%8.3f" % ref_.Hill_n) if ref_ else "       -",
                    ("%12.1f" % ref_.EC50_ppb) if ref_ else "           -"))
    A = pd.DataFrame(rows)
    if "Imax_ref" in A and A["Imax_ref"].notna().any():
        j = A.dropna(subset=["Imax_ref"])
        L.append("\n与 EC50_summary.csv 一致性（%d 个可比）:" % len(j))
        for col, rc in [("Imax", "Imax_ref"), ("Hill_n", "n_ref"), ("EC50", "EC50_ref")]:
            ratio = j[col] / j[rc]
            L.append("   %-8s 倍数中位数=%.3f  几何均值=%.3f  范围 %.3f–%.3f  |倍数-1|<5%%: %d/%d"
                     % (col, ratio.median(), np.exp(np.log(ratio).mean()),
                        ratio.min(), ratio.max(), int((abs(ratio - 1) < 0.05).sum()), len(j)))
    return A

A1 = run(tab.set_index("key")["y"].to_dict(), "用论文《单体数据平均值》")
A2 = run(my.to_dict(), "用我重算的干净均值")
A1.to_csv(os.path.join(OUT, "refit_from_table.csv"), index=False, encoding="utf-8-sig")
A2.to_csv(os.path.join(OUT, "refit_from_clean.csv"), index=False, encoding="utf-8-sig")

open(r"C:\Users\24252\AppData\Local\Temp\opencode\refit.txt", "w", encoding="utf-8").write("\n".join(L))
print("ok")
