# -*- coding: utf-8 -*-
"""三张图 v2：气泡图(标签避让+量表外推区) / 加和性判定 / Ridgeline"""
import os, re, math, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from scipy.stats import gaussian_kde, pearsonr

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

def pow10fmt(v, p):
    if v <= 0: return ""
    k = int(round(math.log10(v)))
    if abs(10.0 ** k - v) / v > 0.05: return "%g" % v
    return r"$10^{%d}$" % k

def logaxis(ax, lo, hi, axis="x"):
    e0, e1 = int(math.floor(math.log10(lo))), int(math.ceil(math.log10(hi)))
    span = e1 - e0
    step = 1 if span <= 4 else (2 if span <= 7 else 3)
    ticks = [10.0 ** k for k in range(e0, e1 + 1, step)]
    if axis == "x":
        ax.set_xscale("log"); ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FuncFormatter(pow10fmt))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=100))
        ax.xaxis.set_minor_formatter(NullFormatter())

def rect_overlap(a, b, pad=2.0):
    return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0 or
                a.y1 + pad < b.y0 or b.y1 + pad < a.y0)

def declutter(ax, pts, texts, fontsize=6.8, dy=11, dx=8, obstacles=None, radii=None):
    """贪心避让：给每个标签找一个不与已放置框(含气泡本身)重叠的偏移"""
    from matplotlib.transforms import Bbox
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    placed, kept = [], []
    if obstacles: placed.extend(list(obstacles))
    if radii is not None and len(radii) == len(pts):
        dp = ax.transData.transform(list(pts))
        r_px = np.sqrt(np.asarray(radii, float) / np.pi) * fig.dpi / 72.0
        for (x, y), rr in zip(dp, r_px):
            m = rr + 5.0                      # 额外留白 5px，避免标签贴着气泡边
            placed.append(Bbox([[x - m, y - m], [x + m, y + m]]))
    offs = [(0, dy), (0, -dy), (dx, dy), (-dx, dy), (dx, -dy), (-dx, -dy),
            (0, 2 * dy), (0, -2 * dy), (2 * dx, dy), (-2 * dx, dy),
            (2 * dx, -dy), (-2 * dx, -dy), (0, 3 * dy), (0, -3 * dy),
            (3 * dx, 2 * dy), (-3 * dx, 2 * dy)]
    for (x, y), txt in zip(pts, texts):
        chosen = None
        for (ox, oy) in offs:
            t = ax.annotate(txt, (x, y), xytext=(ox, oy), textcoords="offset points",
                            ha="center", va="center", fontsize=fontsize, zorder=6)
            bb = t.get_window_extent(renderer=rend)
            if not any(rect_overlap(bb, p) for p in placed):
                chosen = (t, bb); break
            t.remove()
        if chosen is None:
            # 兜底：远处找空位 + 引线
            import matplotlib.transforms as mtrans
            x0, y0 = ax.transData.transform((x, y))
            for dist in (30, 40, 50, 62, 76, 92):
                for ang in range(0, 360, 15):
                    ox = dist * math.cos(math.radians(ang))
                    oy = dist * math.sin(math.radians(ang))
                    t = ax.annotate(txt, (x, y),
                                    xytext=(ox, oy), textcoords="offset points",
                                    ha="center", va="center", fontsize=fontsize - 0.4,
                                    zorder=7,
                                    arrowprops=dict(arrowstyle="-", lw=0.45,
                                                    color="#777777", shrinkA=1, shrinkB=3))
                    bb = t.get_window_extent(renderer=rend)
                    if not any(rect_overlap(bb, p, pad=1.0) for p in placed):
                        chosen = (t, bb); break
                    t.remove()
                if chosen: break
            if chosen is None:
                t = ax.annotate(txt, (x, y), xytext=(0, dy), textcoords="offset points",
                                ha="center", va="center", fontsize=fontsize - 0.6,
                                alpha=0.55, zorder=6)
                chosen = (t, t.get_window_extent(renderer=rend))
        placed.append(chosen[1]); kept.append(chosen[0])
    return kept

# =====================================================================
# 图 1 —— 气泡图
# =====================================================================
ec = pd.read_csv(os.path.join(OUT, "EC50_summary.csv"), encoding="utf-8-sig").dropna(
    subset=["EC50_ppb", "Imax", "Hill_n"])

fig, ax = plt.subplots(figsize=(9.2, 6.3))
IMAX_SCALE = 12.0
ax.axhspan(IMAX_SCALE, max(23, ec["Imax"].max() * 1.18), color="#FDECEC", zorder=0)
ax.scatter(ec["EC50_ppb"], ec["Imax"], s=ec["Hill_n"] * 95 + 45,
           c=ec["R2"], cmap="RdYlGn", vmin=0.94, vmax=1.0,
           edgecolors="#333333", linewidths=0.7, alpha=0.9, zorder=3)
ax.axhline(IMAX_SCALE, ls=":", lw=1.1, color="#B2182B", zorder=2)
t1 = ax.text(0.985, IMAX_SCALE, "量表上限 12", transform=ax.get_yaxis_transform(),
             va="bottom", ha="right", fontsize=7.5, color="#B2182B", zorder=6)
t2 = ax.text(0.015, 1.5, "注：红色区域为外推区 —— Imax > 12 不可能真实发生",
             transform=ax.get_yaxis_transform(), va="bottom", ha="left",
             fontsize=8, color="#8A1520", style="italic", zorder=6)
# 关键：必须先确定坐标范围，标签避让算出的像素位置才与最终渲染一致
logaxis(ax, ec["EC50_ppb"].min() * 0.5, ec["EC50_ppb"].max() * 2.2)
ax.set_xlim(ec["EC50_ppb"].min() * 0.5, ec["EC50_ppb"].max() * 3.0)
ax.set_ylim(0, max(23, ec["Imax"].max() * 1.18))
fig.canvas.draw()
_rend = fig.canvas.get_renderer()
obs = [t1.get_window_extent(renderer=_rend), t2.get_window_extent(renderer=_rend)]
declutter(ax, list(zip(ec["EC50_ppb"], ec["Imax"])), list(ec["Compound"]),
          fontsize=7.2, obstacles=obs, radii=ec["Hill_n"] * 95 + 45)
ax.set_xlabel(r"EC$_{50}$ — potency (ppb, log scale)    越靠左 = 越强")
ax.set_ylabel("Imax — efficacy (0–12)")
ax.set_title("气泡图：效价 × 效能 × 陡峭度，三参数一张图", fontsize=11.5, pad=10)
ax.grid(True, which="major", lw=.35, color="#E4E4E4", zorder=0); ax.set_axisbelow(True)
for k, lab in zip([0.4, 1.5, 4.0], ["n = 0.4 (缓)", "n = 1.5", "n = 4 (陡)"]):
    ax.scatter([], [], s=k * 95 + 45, c="#888888", edgecolors="#333333",
               linewidths=.7, alpha=.9, label=lab)
ax.legend(title="气泡大小 = Hill n", fontsize=7.5, title_fontsize=8,
          loc="lower right", frameon=True, framealpha=.95, borderpad=0.7)
sc = ax.collections[0]
cb = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.042)
cb.set_label("R²（拟合优度）", fontsize=8); cb.ax.tick_params(labelsize=7.5)
fig.savefig(os.path.join(OUT, "fig_bubble.png"), dpi=300)
plt.close(fig)
L.append("[图1] 气泡图 ok；Imax>12 的化合物: " +
         ", ".join(ec.loc[ec["Imax"] > 12, "Compound"]))

# =====================================================================
# 图 2 —— 加和性判定（预测 vs 实测）
# =====================================================================
z = pd.read_excel(os.path.join(BASE, "单体混合感官特征理化数据.xlsx"), sheet_name="组合数据")
zcol = list(z.columns)
cnames = [norm(c) for c in zcol[2:27]]
anames = [str(c) for c in zcol[27:43]]
KEY = ["二甲基三硫", "乙酸异戊酯", "丁酸乙酯", "乙缩醛", "己酸", "正己醇"]
sub = z[z[zcol[0]].astype(str).str.startswith(("样品1-", "母液"))].copy()
brow = sub[sub[zcol[0]].astype(str).str.startswith("母液1")].iloc[0]
c0 = {k: float(brow[zcol[2 + cnames.index(k)]]) for k in KEY}

def parse(row, code):
    hi = tuple(k for k in KEY if float(row[zcol[2 + cnames.index(k)]]) > 1.3 * c0[k])
    cata = {a: (float(row[zcol[27 + j]]) if pd.notna(row[zcol[27 + j]]) else np.nan)
            for j, a in enumerate(anames)}
    return dict(code=code, n=len(hi), hi=hi, cata=cata)

recs = [parse(sub.iloc[i], str(sub.iloc[i, 0])) for i in range(len(sub))]
base = [r for r in recs if r["n"] == 0][0]
singles = {r["hi"][0]: r for r in recs if r["n"] == 1}
bins = [r for r in recs if r["n"] == 2]
p0 = {a: base["cata"][a] for a in anames}

rows = []
for b in bins:
    a1, a2 = b["hi"]
    for at in anames:
        if at == "其他": continue
        x = [p0[at], singles[a1]["cata"][at], singles[a2]["cata"][at], b["cata"][at]]
        if not all(np.isfinite(x)): continue
        pred = min(max(p0[at] + (x[1] - p0[at]) + (x[2] - p0[at]), 0.0), 1.0)
        obs = x[3]
        rows.append(dict(pair="%s+%s" % (a1[:3], a2[:3]), attr=at,
                         pred=pred, obs=obs, dev=obs - pred))
D = pd.DataFrame(rows)
r, p = pearsonr(D["pred"], D["obs"])
L.append("\n[图2] 点数=%d  r=%.3f (p=%.3g)  平均偏差=%+.4f  |偏差|均值=%.4f  RMSE=%.4f"
         % (len(D), r, p, D["dev"].mean(), D["dev"].abs().mean(),
            math.sqrt((D["dev"] ** 2).mean())))
summ = D.groupby("attr").agg(n=("dev", "size"), mean_dev=("dev", "mean"),
                             sd=("dev", "std")).sort_values("mean_dev")
L.append(summ.round(4).to_string())

fig = plt.figure(figsize=(10.4, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.0], wspace=0.06)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axA.set_facecolor("#FAFAFA")
lim = 1.0
axA.plot([0, lim], [0, lim], ls="--", lw=1.5, color="#B2182B", zorder=2)
axA.text(0.97, 0.985, "完美加和 y = x", transform=axA.transAxes, ha="right", va="top",
         fontsize=8.5, color="#B2182B", rotation=0)
attrs = sorted(D["attr"].unique())
pal = plt.get_cmap("tab20")(np.linspace(0, 1, 20))[:len(attrs)]
for at, cc in zip(attrs, pal):
    d = D[D["attr"] == at]
    axA.scatter(d["pred"], d["obs"], s=30, color=cc, alpha=.85,
                edgecolors="white", linewidths=.5, zorder=3, label=at)
# 分箱趋势
D["bin"] = pd.cut(D["pred"], bins=np.linspace(0, 1, 11))
tr = D.groupby("bin", observed=True).agg(x=("pred", "mean"), y=("obs", "mean"),
                                          s=("obs", "size")).dropna()
tr = tr[tr["s"] >= 5]          # 样本太少的分箱不画，避免末端假性下折
axA.plot(tr["x"], tr["y"], "o-", ms=6, lw=2, color="#111111", zorder=5,
         label="分箱均值（实测）")
axA.axhline(1, lw=.6, color="#999999"); axA.axvline(1, lw=.6, color="#999999")
axA.set_xlim(-0.02, 1.04); axA.set_ylim(-0.02, 1.04)
axA.set_xlabel("加和预测频率  $p_0 + (p_A-p_0) + (p_B-p_0)$")
axA.set_ylabel("实测混合频率 $p_{AB}$   (CATA, 18 人)")
axA.set_title("A  加和性判定：二元混合的实测 vs 加和预测", fontsize=11, pad=8)
axA.grid(True, lw=.3, color="#DDDDDD", zorder=0); axA.set_axisbelow(True)
axA.legend(fontsize=6.6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.16),
           frameon=False, title="属性 / 趋势", title_fontsize=7.5)
axA.annotate("n = %d 个点（15 组二元 × 15 属性）\nr = %.3f   RMSE = %.3f"
             % (len(D), r, math.sqrt((D["dev"] ** 2).mean())),
             xy=(0.03, 0.96), xycoords="axes fraction", va="top", ha="left",
             fontsize=8.5, bbox=dict(boxstyle="round,pad=0.35", fc="white",
                                     ec="#BBBBBB", lw=.6, alpha=.92))

# 右：各属性系统偏差
s = summ.reset_index().sort_values("mean_dev")
cols = ["#B2182B" if v > 0 else "#2166AC" for v in s["mean_dev"]]
axB.barh(np.arange(len(s)), s["mean_dev"], color=cols, alpha=.85,
         xerr=s["sd"] / np.sqrt(s["n"]), error_kw=dict(lw=.8, ecolor="#555555", capsize=2))
axB.axvline(0, color="#111111", lw=1.2)
axB.set_yticks(np.arange(len(s))); axB.set_yticklabels(s["attr"], fontsize=8)
axB.set_xlabel("平均偏差  实测 − 加和预测")
axB.set_title("B  哪些属性系统性偏离加和", fontsize=11, pad=8)
axB.grid(True, axis="x", lw=.3, color="#DDDDDD"); axB.set_axisbelow(True)
axB.text(0.03, 0.04, "→ 协同/增强\n← 掩蔽/减弱", transform=axB.transAxes,
         fontsize=8, ha="left", va="bottom")
fig.suptitle("σ–τ 等效线思路的可行替代：基于 CATA 频次的加和性判定（响应非强度，见正文说明）",
             fontsize=11.5, y=1.02)
fig.savefig(os.path.join(OUT, "fig_additivity.png"), dpi=300)
plt.close(fig)
L.append("[图2] 加和性判定图 ok")

# =====================================================================
# 图 3 —— Ridgeline
# =====================================================================
raw = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"), sheet_name="单体数据原始数据")
raw = raw.rename(columns={raw.columns[2]: "I"})
raw["key"] = raw["Sample"].ffill()
raw["cmpd"] = raw["key"].map(base_name)
raw["lvl"] = raw["key"].map(level)
raw["I"] = pd.to_numeric(raw["I"], errors="coerce")
raw = raw[raw["cmpd"].ne("") & raw["I"].notna()].copy()
cnt = raw.groupby("cmpd")["lvl"].nunique().sort_values(ascending=False)
picks = list(cnt[cnt >= 6].index)[:6]
L.append("[图3] 入选: " + " | ".join(picks))

grid_x = np.linspace(-0.3, 12.3, 400)
n = len(picks)
anom = []
fig, axes = plt.subplots(n, 1, figsize=(7.6, 2.05 * n + 0.7), sharex=True)
axes = np.atleast_1d(axes)
STEP, AMP = 1.0, 0.86
for ri, g in enumerate(picks):
    ax = axes[ri]
    d = raw[raw["cmpd"] == g]
    lv = sorted(d["lvl"].dropna().unique())
    base_y = 0.0
    for j, k in enumerate(lv):
        v = d.loc[d["lvl"] == k, "I"].to_numpy(float)
        if len(v) < 3: continue
        try:
            dens = gaussian_kde(v, bw_method=0.13)(grid_x)
        except Exception:
            dens = np.zeros_like(grid_x)
        dens = np.clip(dens, 0, None)
        dens = dens / max(dens.max(), 1e-9)
        cc = plt.get_cmap("viridis")(j / max(len(lv) - 1, 1))
        ax.fill_between(grid_x, base_y, base_y + dens * AMP, color=cc,
                        alpha=.80, lw=.7, edgecolor="white", zorder=3)
        ax.plot(grid_x, base_y + dens * AMP, lw=.7, color="#2B2B2B", zorder=4)
        m = v.mean()
        if len(v) != 18: anom.append("%s C%d 重复 %d 条" % (g, k, len(v)))
        ax.plot([m, m], [base_y, base_y + AMP * 0.99], lw=1.0, ls="--",
                color="#B2182B", zorder=5)
        ax.text(12.45, base_y + 0.06, "C%d  n=%d  $\\bar{x}$=%.1f" % (k, len(v), m),
                fontsize=6.2, va="bottom", ha="left")
        base_y += STEP
    ax.set_ylim(0, base_y + 0.55); ax.set_yticks([])
    ax.set_ylabel(g, rotation=0, ha="right", va="center", fontsize=9, labelpad=10)
    for s_ in ["top", "right", "left"]: ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_bounds(0, 12)
ax = axes[0]
ax.plot([], [], color="#B2182B", ls="--", lw=1, label="该浓度下的均值 $\\bar{x}$")
ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0.06, 1.16))
axes[-1].set_xlim(0, 12); axes[-1].set_xlabel("Perceived odour intensity (0–12)")
fig.suptitle("Ridgeline：单条山脊 = 一个浓度下 18 位品评员的强度分布（下→上 = C1→最高级）",
             fontsize=11.5, y=0.999)
fig.subplots_adjust(left=0.155, right=0.775, top=0.945, bottom=0.05, hspace=0.42)
if anom:
    fig.text(0.995, 0.006, "数据异常（自动检出）：" + "；".join(anom),
             ha="right", va="bottom", fontsize=7.5, color="#B2182B")
L.append("[图3] 样本数异常: " + ("；".join(anom) if anom else "无"))
fig.savefig(os.path.join(OUT, "fig_ridgeline.png"), dpi=300)
plt.close(fig)
L.append("[图3] Ridgeline ok")

open(r"C:\Users\24252\AppData\Local\Temp\opencode\gen4_log.txt", "w", encoding="utf-8").write("\n".join(L))
print("done")
