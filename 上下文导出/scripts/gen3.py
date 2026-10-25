# -*- coding: utf-8 -*-
"""三张示范图：气泡图 / σ-τ 等效线图 / Ridgeline 山脊图"""
import os, re, math, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from matplotlib.collections import PathCollection
from scipy.stats import gaussian_kde

BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUT  = os.path.join(BASE, "教学示范图")
os.makedirs(OUT, exist_ok=True)
L = []

def norm(s):
    if s is None: return ""
    s = re.sub(r"\s+", "", str(s))
    for a, b in [("（", "("), ("）", ")")]: s = s.replace(a, b)
    return s

def base_name(s): return re.sub(r"\(C\d+\)$", "", norm(s), flags=re.I)
def level(s):
    m = re.search(r"\(C(\d+)\)$", norm(s), flags=re.I)
    return int(m.group(1)) if m else None

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "font.size": 9,
    "axes.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "figure.facecolor": "white", "axes.facecolor": "white"})

def logaxis(ax, lo, hi, axis="x"):
    e0, e1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    span = e1 - e0
    step = 1 if span <= 4 else (2 if span <= 7 else 3)
    ticks = [10.0 ** k for k in range(e0, e1 + 1, step)]
    if axis == "x":
        ax.set_xscale("log"); ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: "%g" % v))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=100))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        ax.set_yscale("log"); ax.set_yticks(ticks)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "%g" % v))
        ax.yaxis.set_minor_formatter(NullFormatter())

# =====================================================================
# 图 1 —— 气泡图：EC50 × Imax × n
# =====================================================================
ec = pd.read_csv(os.path.join(OUT, "EC50_summary.csv"), encoding="utf-8-sig")
ec = ec.dropna(subset=["EC50_ppb", "Imax", "Hill_n"])
L.append("### 气泡图数据 %d 行" % len(ec))
L.append(ec.to_string(index=False))

cmap_v = plt.get_cmap("viridis")
norm_v = plt.Normalize(ec["Hill_n"].min(), ec["Hill_n"].max())

fig, ax = plt.subplots(figsize=(6.4, 4.6))
sc = ax.scatter(ec["EC50_ppb"], ec["Imax"], s=ec["Hill_n"] * 130 + 60,
                c=ec["R2"], cmap="RdYlGn", vmin=0.9, vmax=1.0,
                edgecolors="#333333", linewidths=0.7, alpha=0.88, zorder=3)
for _, r in ec.iterrows():
    ax.annotate(r["Compound"], (r["EC50_ppb"], r["Imax"]),
                xytext=(0, math.sqrt(r["Hill_n"] * 130 + 60) / 2 + 4),
                textcoords="offset points", ha="center", fontsize=7)
logaxis(ax, ec["EC50_ppb"].min() * 0.5, ec["EC50_ppb"].max() * 2.0)
ax.set_xlim(ec["EC50_ppb"].min() * 0.5, ec["EC50_ppb"].max() * 2.6)
ax.set_ylim(0, max(14.5, ec["Imax"].max() * 1.22))
ax.set_xlabel("EC$_{50}$ — potency (ppb, log scale)   ←  越左越强")
ax.set_ylabel("Imax — efficacy (0–12)")
ax.set_title("气泡图：效价 × 效能 × 陡峭度 三参数一次呈现", fontsize=11, pad=10)
ax.grid(True, which="major", lw=.35, color="#E2E2E2", zorder=0)
ax.set_axisbelow(True)
# 气泡大小图例
for k, val in zip(["弱", "中", "强"], [0.5, 1.5, 3.0]):
    ax.scatter([], [], s=val * 130 + 60, c="#888888", edgecolors="#333333",
               linewidths=.7, alpha=.88, label="n = %s" % k)
ax.legend(title="Hill n（陡峭度）", fontsize=7.5, title_fontsize=8,
          loc="lower left", frameon=True, framealpha=.9)
cb = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.045)
cb.set_label("R²（拟合优度）", fontsize=8); cb.ax.tick_params(labelsize=7.5)
fig.savefig(os.path.join(OUT, "fig_bubble.png"), dpi=300)
plt.close(fig)
L.append("\n[图1 气泡图] 已生成")

# =====================================================================
# 图 2 —— σ-τ 等效线图（基于 CATA 频次）
# =====================================================================
z = pd.read_excel(os.path.join(BASE, "单体混合感官特征理化数据.xlsx"), sheet_name="组合数据")
zcol = list(z.columns)
cnames = [norm(c) for c in zcol[2:27]]
anames = [str(c) for c in zcol[27:43]]
KEY = ["二甲基三硫", "乙酸异戊酯", "丁酸乙酯", "乙缩醛", "己酸", "正己醇"]
for k in KEY:
    assert k in cnames, "missing %s" % k

sub = z[z[zcol[0]].astype(str).str.startswith(("样品1-", "母液"))].copy()
base_rows = sub[sub[zcol[0]].astype(str).str.startswith("母液1")]
assert len(base_rows) == 1
c0 = {k: float(base_rows.iloc[0][zcol[2 + cnames.index(k)]]) for k in KEY}
L.append("\n### 基线（母液1）浓度: " + ", ".join("%s=%.4g" % (k, v) for k, v in c0.items()))

def flags(row):
    return {k: float(row[zcol[2 + cnames.index(k)]]) > 1.3 * c0[k] for k in KEY}

recs = []
for i in range(len(sub)):
    code = str(sub.iloc[i, 0])
    row = sub.iloc[i]
    f = flags(row)
    hi = [k for k in KEY if f[k]]
    cata = {a: (float(row[zcol[27 + j]]) if pd.notna(row[zcol[27 + j]]) else np.nan)
            for j, a in enumerate(anames)}
    recs.append(dict(code=code, n=len(hi), hi=tuple(hi), cata=cata))

byp = {r["code"]: r for r in recs}
L.append("### 高水平化合物数分布: " + str({n: sum(1 for r in recs if r["n"] == n) for n in range(7)}))
singles = {r["hi"][0]: r for r in recs if r["n"] == 1}
bins    = [r for r in recs if r["n"] == 2]
base    = [r for r in recs if r["n"] == 0][0]
L.append("### 单因子样品: " + ", ".join(sorted(singles)))
L.append("### 二元样品数: %d" % len(bins))

# 加和判定：效应 e = p - p_baseline（CATA 频次增量）
pts = []
for b in bins:
    a1, a2 = b["hi"]
    for at in anames:
        if at == "其他": continue
        p0 = base["cata"][at]
        e1 = singles[a1]["cata"][at] - p0
        e2 = singles[a2]["cata"][at] - p0
        e12 = b["cata"][at] - p0
        if not all(np.isfinite([p0, e1, e2, e12])): continue
        if e1 + e2 < 0.06: continue          # 效应太小，σ 不稳定
        tau = e1 / (e1 + e2)
        sig = e12 / (e1 + e2)
        pts.append(dict(A=a1, B=a2, attr=at, tau=tau, sigma=sig, e_sum=e1 + e2, e12=e12))
pt = pd.DataFrame(pts)
L.append("\n### σ-τ 有效点数: %d (来自 %d 个二元组合 × %d 属性)" % (len(pt), len(bins), pt["attr"].nunique()))
L.append(pt.groupby("attr").agg(n=("sigma", "size"), mean_sigma=("sigma", "mean"),
                                sd=("sigma", "std")).round(3).to_string())

attrs = sorted(pt["attr"].unique())
palette = plt.get_cmap("tab20")(np.linspace(0, 1, 20))[:len(attrs)]
colmap = dict(zip(attrs, palette))

fig, ax = plt.subplots(figsize=(6.6, 5.4))
ax.axhspan(1, max(2.6, pt["sigma"].max() * 1.05), color="#FDECEC", zorder=0)
ax.axhspan(0, 1, color="#EAF2FB", zorder=0)
for at in attrs:
    d = pt[pt["attr"] == at]
    ax.scatter(d["tau"], d["sigma"], s=34, color=colmap[at], alpha=.85,
               edgecolors="white", linewidths=.5, zorder=3, label=at)
ax.axhline(1.0, color="#B2182B", lw=1.4, ls="--", zorder=2)
ax.text(0.985, 1.012, "加和线 σ = 1", transform=ax.get_yaxis_transform(),
        ha="right", va="bottom", fontsize=8.5, color="#B2182B")
ax.text(0.985, 1.72, "协同 / 增强  σ > 1", transform=ax.get_yaxis_transform(),
        ha="right", va="center", fontsize=8.5, color="#7A1420")
ax.text(0.985, 0.52, "掩蔽 / 减弱  σ < 1", transform=ax.get_yaxis_transform(),
        ha="right", va="center", fontsize=8.5, color="#12406E")
ax.set_xlabel(r"$\tau_A = e_A/(e_A+e_B)$   — 组分 A 的贡献份额")
ax.set_ylabel(r"$\sigma = e_{AB}/(e_A+e_B)$   — 实测/加和预测")
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0, max(2.6, pt["sigma"].max() * 1.05))
ax.set_title("σ–τ 等效线图：二元混合的加和性判定（响应 = CATA 频次）", fontsize=11, pad=10)
ax.grid(True, which="major", lw=.3, color="#DDDDDD", zorder=0); ax.set_axisbelow(True)
ax.legend(fontsize=6.6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.115),
          frameon=False, title="属性", title_fontsize=7.5)
fig.savefig(os.path.join(OUT, "fig_sigma_tau.png"), dpi=300)
plt.close(fig)
L.append("[图2 σ-τ] 已生成")

# =====================================================================
# 图 3 —— Ridgeline 山脊图（单体逐人强度分布）
# =====================================================================
raw = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"), sheet_name="单体数据原始数据")
raw = raw.rename(columns={raw.columns[2]: "I"})
raw["key"] = raw["Sample"].ffill()
raw["cmpd"] = raw["key"].map(base_name)
raw["lvl"] = raw["key"].map(level)
raw["I"] = pd.to_numeric(raw["I"], errors="coerce")
raw = raw[raw["cmpd"].ne("") & raw["I"].notna()].copy()
L.append("\n### Ridgeline 用到逐人记录: %d 条, 化合物 %d 个" % (len(raw), raw["cmpd"].nunique()))

cnt = raw.groupby("cmpd")["lvl"].nunique().sort_values(ascending=False)
L.append("### 各化合物级数: " + ", ".join("%s=%d" % (k, v) for k, v in cnt.items()))
picks = list(cnt[cnt >= 6].index)[:6]
if len(picks) < 4: picks = list(cnt.index[:6])
L.append("### 入选化合物: " + " | ".join(picks))

grid_x = np.linspace(0, 12, 300)
nrows = len(picks)
fig, axes = plt.subplots(nrows, 1, figsize=(7.2, 1.75 * nrows + 0.6), sharex=True)
axes = np.atleast_1d(axes)
for r_i, g in enumerate(picks):
    ax = axes[r_i]
    d = raw[raw["cmpd"] == g]
    lv = sorted(d["lvl"].dropna().unique())
    base_y = 0.0
    for j, k in enumerate(lv):
        v = d.loc[d["lvl"] == k, "I"].to_numpy(float)
        if len(v) < 3: continue
        try:
            kde = gaussian_kde(v, bw_method=0.35)
            dens = kde(grid_x)
        except Exception:
            dens = np.zeros_like(grid_x)
        dens = dens / max(dens.max(), 1e-9)
        cc = plt.get_cmap("viridis")(j / max(len(lv) - 1, 1))
        ax.fill_between(grid_x, base_y, base_y + dens * 0.82, color=cc,
                        alpha=.72, lw=.6, edgecolor="white", zorder=3)
        ax.plot(grid_x, base_y + dens * 0.82, lw=.6, color="#333333", zorder=4)
        if j % 2 == 0:
            ax.text(12.35, base_y + 0.40, "C%d" % k, fontsize=6.4, va="center", ha="left")
        base_y += 1.02
    ax.set_ylim(0, base_y + 0.15)
    ax.set_yticks([])
    ax.set_ylabel(g, rotation=0, ha="right", va="center", fontsize=8.5, labelpad=8)
    for s in ["top", "right", "left"]: ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_bounds(0, 12)
    ax.tick_params(axis="y", length=0)
axes[-1].set_xlim(0, 12)
axes[-1].set_xlabel("Perceived odour intensity (0–12)")
fig.suptitle("Ridgeline：每个浓度下 18 位品评员的强度分布（自下而上 = C1→C7）",
             fontsize=11, y=1.002)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(os.path.join(OUT, "fig_ridgeline.png"), dpi=300)
plt.close(fig)
L.append("[图3 Ridgeline] 已生成")

open(r"C:\Users\24252\AppData\Local\Temp\opencode\gen3_log.txt", "w", encoding="utf-8").write("\n".join(L))
print("done")
