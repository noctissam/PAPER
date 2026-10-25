# -*- coding: utf-8 -*-
"""A-2: 用清理后的数据重画全部 14 个化合物的 Ridgeline（分 2 张）"""
import os, math, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUT  = os.path.join(BASE, "教学示范图")
L = []

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "font.size": 9,
    "axes.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "figure.facecolor": "white", "axes.facecolor": "white"})

df = pd.read_csv(os.path.join(OUT, "clean_assessor_level.csv"), encoding="utf-8-sig")
L.append("清理后：%d 条（每人每浓度 1 值）, 化合物 %d, key %d, 品评员 %d"
         % (len(df), df["cmpd"].nunique(), df["key"].nunique(), df["Subject"].nunique()))

stats = (df.groupby("cmpd").agg(n_lvl=("lvl", "nunique"), n=("I", "size"),
                                mean=("I", "mean"), mx=("I", "max"))
           .sort_values(["n_lvl", "cmpd"], ascending=[False, True]))
L.append("\n各化合物：")
L.append(stats.to_string())
L.append("\n越界(已截断到12)记录: %d 条" % int(df["oob"].sum()))

order = list(stats.index)
half = math.ceil(len(order) / 2)
GRID_X = np.linspace(-0.3, 12.3, 400)
STEP, AMP = 1.0, 0.86

def draw(picks, fname, part):
    m = len(picks)
    fig, axes = plt.subplots(m, 1, figsize=(7.6, 2.05 * m + 0.7), sharex=True)
    axes = np.atleast_1d(axes)
    anom = []
    for ri, g in enumerate(picks):
        ax = axes[ri]
        d = df[df["cmpd"] == g]
        lv = sorted(d["lvl"].dropna().unique())
        base_y = 0.0
        for j, k in enumerate(lv):
            v = d.loc[d["lvl"] == k, "I"].to_numpy(float)
            if len(v) < 3: continue
            try:
                dens = gaussian_kde(v, bw_method=0.13)(GRID_X)
            except Exception:
                dens = np.zeros_like(GRID_X)
            dens = np.clip(dens, 0, None)
            dens = dens / max(dens.max(), 1e-9)
            cc = plt.get_cmap("viridis")(j / max(len(lv) - 1, 1))
            ax.fill_between(GRID_X, base_y, base_y + dens * AMP, color=cc,
                            alpha=.80, lw=.7, edgecolor="white", zorder=3)
            ax.plot(GRID_X, base_y + dens * AMP, lw=.7, color="#2B2B2B", zorder=4)
            mval = v.mean()
            ax.plot([mval, mval], [base_y, base_y + AMP * 0.99], lw=1.0, ls="--",
                    color="#B2182B", zorder=5)
            ax.text(12.45, base_y + 0.06, "C%d  n=%d  $\\bar{x}$=%.1f" % (k, len(v), mval),
                    fontsize=6.2, va="bottom", ha="left")
            if len(v) != 18: anom.append("%s C%d n=%d" % (g, k, len(v)))
            base_y += STEP
        ax.set_ylim(0, base_y + 0.55); ax.set_yticks([])
        ax.set_ylabel(g, rotation=0, ha="right", va="center", fontsize=9, labelpad=10)
        for s_ in ["top", "right", "left"]: ax.spines[s_].set_visible(False)
        ax.spines["bottom"].set_bounds(0, 12)
    axes[0].plot([], [], color="#B2182B", ls="--", lw=1, label="该浓度下的均值 $\\bar{x}$")
    axes[0].legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 1.16))
    axes[-1].set_xlim(0, 12)
    axes[-1].set_xlabel("Perceived odour intensity (0–12)")
    fig.suptitle("Ridgeline %s：单条山脊 = 一个浓度下 18 位品评员的强度分布（下→上 = C1→最高级）"
                 % part, fontsize=11.5, y=0.999)
    fig.subplots_adjust(left=0.155, right=0.775, top=0.945, bottom=0.05, hspace=0.42)
    if anom:
        fig.text(0.995, 0.006, "样本数异常：" + "；".join(anom),
                 ha="right", va="bottom", fontsize=7.5, color="#B2182B")
    fig.savefig(os.path.join(OUT, fname), dpi=300)
    plt.close(fig)
    L.append("[Ridgeline %s] %d 化合物 -> %s ; 异常=%s"
             % (part, m, fname, anom if anom else "无"))

draw(order[:half], "fig_ridgeline_A.png", "1/2")
draw(order[half:], "fig_ridgeline_B.png", "2/2")

# 双峰指标：Sarle 偏度-峰度系数
from scipy.stats import skew, kurtosis
rows = []
for (g, k), d in df.groupby(["cmpd", "lvl"]):
    v = d["I"].to_numpy(float)
    if len(v) < 5: continue
    sk = skew(v, bias=False); ku = kurtosis(v, fisher=False, bias=False); n = len(v)
    bc = (sk ** 2 + 1) / (ku + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    rows.append(dict(cmpd=g, lvl=k, n=n, mean=v.mean(), sd=v.std(ddof=1),
                     skew=sk, kurt=ku, bimodality=bc))
B = pd.DataFrame(rows).sort_values("bimodality", ascending=False)
B.to_csv(os.path.join(OUT, "bimodality_by_level.csv"), index=False, encoding="utf-8-sig")
_nb = int((B["bimodality"] > 0.555).sum())
L.append("\n双峰系数 BC > 0.555（Sarle 经验阈值）的浓度点: %d / %d (%.0f%%)"
         % (_nb, len(B), 100.0 * _nb / max(len(B), 1)))
L.append("BC 前 15 名:")
L.append(B.head(15).to_string(index=False))

open(os.path.join(OUT, "QC_双峰与命名.txt"), "w", encoding="utf-8").write("\n".join(L))
open(r"C:\Users\24252\AppData\Local\Temp\opencode\ridge_log.txt", "w",
     encoding="utf-8").write("\n".join(L))
print("done")
