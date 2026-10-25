# -*- coding: utf-8 -*-
"""B-v2: 修正 ΔBIC 符号 + 天花板过滤 + 双特征集聚类 + 真双峰定位"""
import os, math, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUT  = os.path.join(BASE, "教学示范图")
L = []
np.random.seed(0)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "font.size": 9,
    "axes.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "figure.facecolor": "white", "axes.facecolor": "white"})

df = pd.read_csv(os.path.join(OUT, "clean_assessor_level.csv"), encoding="utf-8-sig")
subj = sorted(df["Subject"].unique()); cmps = sorted(df["cmpd"].unique())

# ================= 特征集 1：平均响应 =================
P = df.pivot_table(index="Subject", columns="cmpd", values="I", aggfunc="mean")
P = P.reindex(index=subj, columns=cmps)

# ================= 特征集 2：检出阈值（每人每化合物首次≥6的最低浓度级） =================
hi = df[df["I"] >= 6]
T = hi.groupby(["Subject", "cmpd"])["lvl"].min().unstack()
full = df.groupby("cmpd")["lvl"].max()
T = T.reindex(index=subj, columns=cmps)
for c in cmps:
    T[c] = T[c].fillna(full[c] + 1)          # 从未达阈 = 最高级+1
T = T.astype(float)

Zr = StandardScaler().fit_transform(P)        # response
Zt = StandardScaler().fit_transform(T)        # threshold

gm_mean = P.mean(axis=1)

def sil_scan(Z, tag, rows):
    for k in range(2, 7):
        lab = KMeans(n_clusters=k, n_init=100, random_state=0).fit_predict(Z)
        rows.append(dict(features=tag, k=k,
                          silhouette=silhouette_score(Z, lab),
                          sizes=",".join(str(int((lab == i).sum())) for i in range(k))))
    return rows

L.append("=" * 78)
L.append("【B1】特征集比较 —— silhouette（>0.5 才算有真实簇；本数据 n=18）")
L.append("=" * 78)
rows = []
sil_scan(Zr, "平均响应", rows); sil_scan(Zt, "检出阈值", rows)
S = pd.DataFrame(rows)
L.append(S.to_string(index=False))
best = S.loc[S["silhouette"].idxmax()]
L.append(">>> 最优: 特征=%s, k=%d, silhouette=%.3f" % (best["features"], best["k"], best["silhouette"]))
L.append(">>> 判读: silhouette 全部 %.3f–%.3f，均远低于 0.5 → 聚类结构很弱/不存在"
         % (S["silhouette"].min(), S["silhouette"].max()))

# 用最优特征集做正式聚类
Zb = Zt if best["features"] == "检出阈值" else Zr
best_k = int(best["k"])
lab = KMeans(n_clusters=best_k, n_init=200, random_state=0).fit_predict(Zb)
pca = PCA(n_components=2).fit(Zb); X = pca.transform(Zb)

# ================= B2: GMM BIC（符号语义修正） =================
L.append("\n" + "=" * 78)
L.append("【B2】18 人全局平均强度的 GMM 组分数选择")
L.append("=" * 78)
x = gm_mean.to_numpy().reshape(-1, 1)
br = []
for k in range(1, 7):
    g = GaussianMixture(k, n_init=50, random_state=0).fit(x)
    br.append(dict(k=k, BIC=g.bic(x)))
BI = pd.DataFrame(br)
# 关键修正：ΔBIC = BIC1 − BICk，**正值**才表示 k 组分更优
BI["dBIC_favors_k"] = BI.loc[BI["k"] == 1, "BIC"].iloc[0] - BI["BIC"]
L.append(BI.round(2).to_string(index=False))
kb = int(BI.loc[BI["BIC"].idxmin(), "k"])
d2 = float(BI.loc[BI["k"] == 2, "dBIC_favors_k"].iloc[0])
L.append(">>> BIC 最优组分数 = %d" % kb)
L.append(">>> ΔBIC = BIC(k=1) − BIC(k=2) = %+.2f" % d2)
L.append("    语义: 正值→2组分更优；负值→1组分更优。Kass 规则: >10 强, 6–10 中, <6 弱")
L.append("    本数据 = %+.2f  →  %s"
         % (d2, "支持双峰" if d2 > 6 else ("弱证据" if d2 > 0 else "不支持双峰（1组分更优）")))
L.append("    全局均值范围仅 %.2f–%.2f（sd=%.2f），量程占比极小 → 天然难以分群"
         % (gm_mean.min(), gm_mean.max(), gm_mean.std(ddof=1)))

# ================= B3: 逐浓度点双峰 + 天花板过滤 =================
L.append("\n" + "=" * 78)
L.append("【B3】逐浓度点双峰检验 —— 必须先剔除天花板效应")
L.append("=" * 78)
rows = []
for (g, lv), d in df.groupby(["cmpd", "lvl"]):
    v = d["I"].to_numpy().reshape(-1, 1)
    if len(v) < 8 or np.unique(v).size < 3: continue
    g1 = GaussianMixture(1, n_init=50, random_state=0).fit(v)
    g2 = GaussianMixture(2, n_init=50, random_state=0).fit(v)
    ceil = float((v.ravel() == 12).mean())
    rows.append(dict(cmpd=g, lvl=lv, mean=v.mean(), sd=v.std(ddof=1),
                     ceil_frac=ceil,
                     dBIC=g1.bic(v) - g2.bic(v),   # 正值 = 2 组分更优
                     w_min=min(g2.weights_),
                     vals=",".join("%g" % t for t in v.ravel())))
BD = pd.DataFrame(rows)
# 天花板判定：>=50% 观测顶格 12
BD["ceiling"] = BD["ceil_frac"] >= 0.5
BD = BD.sort_values("dBIC", ascending=False)
BD.to_csv(os.path.join(OUT, "bimodality_GMM_by_level.csv"), index=False, encoding="utf-8-sig")

n_all = int((BD["dBIC"] > 10).sum())
cl = BD[BD["ceiling"]]
nn = BD[~BD["ceiling"]]
n_cl = int((cl["dBIC"] > 10).sum()); n_nn = int((nn["dBIC"] > 10).sum())
L.append("全部 87 点中 ΔBIC>10 = %d (%.0f%%)" % (n_all, 100 * n_all / len(BD)))
L.append("  ├─ 天花板点 (≥50%% 观测=12): %d 个，其中 ΔBIC>10 = %d (%.0f%%)  ← 伪双峰"
         % (len(cl), n_cl, 100 * n_cl / max(len(cl), 1)))
L.append("  └─ 非天花板点            : %d 个，其中 ΔBIC>10 = %d (%.0f%%)  ← 真双峰候选"
         % (len(nn), n_nn, 100 * n_nn / max(len(nn), 1)))
L.append("\n【真双峰候选（非天花板）按 ΔBIC 排序】")
for _, r in nn[nn["dBIC"] > 10].head(12).iterrows():
    L.append("   %-14s C%-2d mean=%5.2f sd=%4.2f ceil=%.0f%% ΔBIC=%+7.2f w_min=%.2f"
             % (r["cmpd"], r["lvl"], r["mean"], r["sd"], 100 * r["ceil_frac"], r["dBIC"], r["w_min"]))
    L.append("      18 个原始值: " + r["vals"])

# 最强的一例做直方图
top = nn[nn["dBIC"] > 10].iloc[0] if n_nn > 0 else nn.iloc[0]
L.append("\n>>> 代表例: %s C%d  ΔBIC=%+.2f" % (top["cmpd"], top["lvl"], top["dBIC"]))

# ---- B3b 严格双峰判据（同时排除地板/天花板/无间隙）----
L.append("\n" + "=" * 78)
L.append("【B3b】严格双峰判据：5 条同时满足才判为真双峰")
L.append("=" * 78)
L.append("  (1) dBIC > 10            —— 2组分显著优于1组分")
L.append("  (2) 0.20 <= w <= 0.80    —— 小组分不是零星离群点")
L.append("  (3) |mu1-mu2| >= 3 量表单位 —— 两群真正分开")
L.append("  (4) 3 <= mean <= 9       —— 远离地板0与天花板12")
L.append("  (5) 两组分之间存在步长>=3的连续空档 —— 观测中确有断层(>=2个空箱)")
strict = []
for _, r in BD.iterrows():
    vv_ = np.array([float(t) for t in r["vals"].split(",")]).reshape(-1, 1)
    g2s = GaussianMixture(2, n_init=50, random_state=0).fit(vv_)
    mu = sorted(g2s.means_.ravel()); w = min(g2s.weights_)
    u = np.unique(vv_.ravel())
    gap = False; gap_max = 0
    if len(u) >= 2:
        d = np.diff(u)
        i = int(np.argmax(d)); gap_max = int(d[i])
        # 最大空档需落在两个组分均值之间
        gap = bool(gap_max >= 3 and u[i] < mu[1] and u[i + 1] > mu[0])
    ok = bool((r["dBIC"] > 10) and (0.20 <= w <= 0.80) and (abs(mu[1] - mu[0]) >= 3)
              and (3 <= r["mean"] <= 9) and gap)
    strict.append(dict(cmpd=r["cmpd"], lvl=r["lvl"], mean=r["mean"], sd=r["sd"],
                       dBIC=r["dBIC"], w=w, mu_diff=abs(mu[1] - mu[0]),
                       gap=gap, gap_span=gap_max, strict=ok,
                       c1=bool(r["dBIC"] > 10), c2=bool(0.20 <= w <= 0.80),
                       c3=bool(abs(mu[1] - mu[0]) >= 3), c4=bool(3 <= r["mean"] <= 9),
                       c5=bool(gap), vals=r["vals"]))
ST = pd.DataFrame(strict).sort_values("dBIC", ascending=False)
ST.to_csv(os.path.join(OUT, "bimodality_strict.csv"), index=False, encoding="utf-8-sig")
n_strict = int(ST["strict"].sum())
L.append("")
L.append("严格判据通过: %d / %d 个浓度点 (%.0f%%)" % (n_strict, len(ST), 100 * n_strict / max(len(ST), 1)))
_n_all = int((BD["dBIC"] > 10).sum())
L.append("从 %d 个 dBIC>10 的点中淘汰原因（可重叠）:" % _n_all)
L.append("   权重失衡 (w<0.2 或 >0.8)  : %d"
         % int(((BD["dBIC"] > 10) & ((BD["w_min"] < 0.2) | (BD["w_min"] > 0.8))).sum()))
L.append("   均值落在地板/天花板区     : %d"
         % int(((BD["dBIC"] > 10) & ((BD["mean"] < 3) | (BD["mean"] > 9))).sum()))
_l = ST[ST["c1"]]
L.append("\n逐条通过率（在 %d 个 dBIC>10 的点内）:" % len(_l))
for cc, nm in [("c1", "(1) dBIC>10"), ("c2", "(2) 权重 0.2-0.8"),
               ("c3", "(3) 组分均值差>=3"), ("c4", "(4) 3<=mean<=9"),
               ("c5", "(5) 存在>=3的空档")]:
    L.append("   %-22s %2d / %2d" % (nm, int(_l[cc].sum()), len(_l)))
_l2 = _l[_l[["c2", "c3", "c4", "c5"]].all(axis=1)]
_l3 = _l[_l[["c2", "c3", "c4"]].all(axis=1)]
L.append("   同时满足 (2)(3)(4)(5)  : %d / %d" % (len(_l2), len(_l)))
L.append("   满足 (2)(3)(4) 但无空档 : %d / %d  ← 连续偏态而非双峰"
         % (len(_l3) - len(_l2), len(_l)))
L.append("\n【通过严格判据的点】")
for _, r in ST[ST["strict"]].iterrows():
    L.append("   %-14s C%-2d mean=%5.2f sd=%4.2f dBIC=%+7.2f mu差=%.1f 空档跨度=%d"
             % (r["cmpd"], r["lvl"], r["mean"], r["sd"], r["dBIC"], r["mu_diff"], r["gap_span"]))
    L.append("      18 个原始值: " + r["vals"])
if n_strict > 0:
    top = ST[ST["strict"]].iloc[0]
    L.append("\n>>> 严格判据下最强例: %s C%d  dBIC=%+.2f" % (top["cmpd"], top["lvl"], top["dBIC"]))

# ================= 图 =================
fig = plt.figure(figsize=(14.0, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.0, 1.0], hspace=0.48, wspace=0.30)

# A 热图
axA = fig.add_subplot(gs[0, 0])
order_idx = np.argsort(gm_mean.to_numpy())[::-1]
im = axA.imshow(P.to_numpy()[order_idx, :], aspect="auto", cmap="viridis",
                vmin=0, vmax=12, interpolation="nearest")
axA.set_yticks(range(len(subj)))
axA.set_yticklabels([gm_mean.index[i] for i in order_idx], fontsize=6.2)
axA.set_xticks(range(len(cmps)))
axA.set_xticklabels(cmps, rotation=60, ha="right", fontsize=6.2)
axA.set_title("A  品评员 × 化合物 平均强度\n（行按全局均值降序；红线=聚类边界）", fontsize=9.5, pad=6)
cb = fig.colorbar(im, ax=axA, fraction=0.035, pad=0.02)
cb.ax.tick_params(labelsize=6.2); cb.set_label("强度", fontsize=7)
cuts = [i for i in range(1, len(order_idx)) if lab[order_idx[i]] != lab[order_idx[i - 1]]]
for c in cuts: axA.axhline(c - 0.5, color="#E8352A", lw=1.4)

# B PCA
axB = fig.add_subplot(gs[0, 1])
cols = ["#2166AC", "#B2182B", "#1B7837", "#762A83", "#E08214"]
for c in range(best_k):
    m = lab == c
    axB.scatter(X[m, 0], X[m, 1], s=72, color=cols[c % 5], alpha=.88,
                edgecolors="white", lw=.9, label="簇 %d (n=%d)" % (c + 1, m.sum()))
    for s_, r_ in zip(np.array(subj)[m], X[m]):
        axB.annotate(s_, (r_[0], r_[1]), xytext=(4, 3), textcoords="offset points",
                     fontsize=5.8, color="#333333")
axB.set_xlabel("PC1 (%.0f%%)" % (100 * pca.explained_variance_ratio_[0]))
axB.set_ylabel("PC2 (%.0f%%)" % (100 * pca.explained_variance_ratio_[1]))
axB.set_title("B  品评员 PCA + k-means（特征=%s）\nsilhouette=%.2f → 结构很弱"
              % (best["features"], best["silhouette"]), fontsize=9.5, pad=6)
axB.legend(fontsize=7, frameon=True); axB.grid(True, lw=.3, color="#E0E0E0")
axB.set_axisbelow(True)

# C silhouette 曲线（两特征集）
axC = fig.add_subplot(gs[0, 2])
for tag, cc, mk in [("平均响应", "#B2182B", "o"), ("检出阈值", "#2166AC", "s")]:
    d_ = S[S["features"] == tag]
    axC.plot(d_["k"], d_["silhouette"], mk + "-", color=cc, lw=1.6, ms=5, label=tag)
axC.axhline(0.5, ls="--", lw=1.1, color="#777777")
axC.text(2.05, 0.505, "0.5 = 有真实簇的常用下限", fontsize=6.8, color="#555555")
axC.set_xlabel("簇数 k"); axC.set_ylabel("silhouette")
axC.set_title("C  聚类可分性检验（n=18）\n全部 <0.5 → 不存在可靠的敏感/不敏感分组", fontsize=9.5, pad=6)
axC.set_xticks(range(2, 7)); axC.legend(fontsize=7.5); axC.grid(True, lw=.3, color="#E0E0E0")
axC.set_axisbelow(True)

# D GMM BIC
axD = fig.add_subplot(gs[1, 0])
axD.plot(BI["k"], BI["BIC"], "o-", ms=6, color="#B2182B", lw=1.7, label="BIC")
axD.axvline(kb, ls=":", lw=1.3, color="#2166AC")
axD.set_xlabel("高斯组分数 k"); axD.set_ylabel("BIC（越小越好）")
axD.set_title("D  18 人全局均值的组分数选择", fontsize=9.5, pad=6)
axD.set_xticks(list(BI["k"])); axD.legend(fontsize=7.5)
axD.grid(True, lw=.3, color="#E0E0E0"); axD.set_axisbelow(True)
axD.text(0.97, 0.94, "最优 k = %d\nΔBIC = BIC(1) − BIC(2) = %+.2f\n→ %s"
         % (kb, d2, "不支持双峰" if d2 <= 6 else "支持双峰"),
         transform=axD.transAxes, ha="right", va="top", fontsize=8,
         bbox=dict(boxstyle="round,pad=0.35", fc="#FFF6F6", ec="#B2182B", lw=.7))

# E 逐浓度点 ΔBIC vs 均值，天花板标色
axE = fig.add_subplot(gs[1, 1])
axE.scatter(BD.loc[BD["ceiling"], "mean"], BD.loc[BD["ceiling"], "dBIC"],
            s=42, facecolors="#E0E0E0", edgecolors="#999999", lw=.7,
            label="天花板点 (≥50%%=12)  n=%d" % len(cl), zorder=3)
axE.scatter(BD.loc[~BD["ceiling"], "mean"], BD.loc[~BD["ceiling"], "dBIC"],
            s=42, facecolors="#2166AC", edgecolors="white", lw=.7,
            label="非天花板点  n=%d" % len(nn), zorder=4)
axE.axhline(10, ls="--", lw=1.1, color="#B2182B")
axE.text(0.02, 0.97, "ΔBIC=10 强证据线", transform=axE.transAxes, va="top",
         fontsize=7, color="#B2182B")
axE.set_xlabel("该浓度点的平均强度"); axE.set_ylabel("ΔBIC = BIC(1) − BIC(2)")
axE.set_title("E  逐浓度点双峰检验：87 点中 %d 点 dBIC>10，\n加严格判据后仅 %d 点为真双峰"
              % (n_all, n_strict), fontsize=9.5, pad=6)
axE.legend(fontsize=7, frameon=True, loc="upper center"); axE.grid(True, lw=.3, color="#E0E0E0")
axE.set_axisbelow(True)

# F 代表例直方图
axF = fig.add_subplot(gs[1, 2])
vv = np.array([float(t) for t in top["vals"].split(",")])
axF.hist(vv, bins=np.arange(1.5, 13.5, 1.0), color="#C6DBEF", edgecolor="white",
         lw=.8, zorder=2)
g2f = GaussianMixture(2, n_init=50, random_state=0).fit(vv.reshape(-1, 1))
xx = np.linspace(-0.5, 12.5, 300)
axF.plot(xx, np.exp(g2f.score_samples(xx.reshape(-1, 1))), lw=2.0, color="#B2182B", zorder=5,
         label="2 组分拟合")
for mu, sg, w in zip(g2f.means_.ravel(), np.sqrt(g2f.covariances_.ravel()), g2f.weights_):
    axF.plot(xx, w * np.exp(-0.5 * ((xx - mu) / sg) ** 2) / (sg * math.sqrt(2 * math.pi)),
             ls=":", lw=1.3, color="#B2182B", alpha=.8, zorder=4)
    axF.axvline(mu, ls="--", lw=1.0, color="#2166AC", alpha=.8)
    axF.text(mu, axF.get_ylim()[1] * .97, "μ=%.1f" % mu, ha="center", fontsize=7, color="#2166AC")
axF.set_xlabel("强度评分"); axF.set_ylabel("频数")
axF.set_title("F  严格判据下最强真双峰例\n%s C%d  ΔBIC=%+.1f" % (top["cmpd"], top["lvl"], top["dBIC"]),
              fontsize=9.5, pad=6)
axF.legend(fontsize=7.5, frameon=False); axF.grid(True, lw=.3, color="#E0E0E0")
axF.set_axisbelow(True)

fig.suptitle("品评员聚类与双峰性严格检验（18 人 × 87 浓度点 = 1,566 条清理后记录）",
             fontsize=12.5, y=0.995)
fig.savefig(os.path.join(OUT, "fig_cluster_bimodality.png"), dpi=300)
plt.close(fig)

# ================= 明细 =================
L.append("\n" + "=" * 78)
L.append("【B4】聚类明细（特征=%s, k=%d）" % (best["features"], best_k))
L.append("=" * 78)
cl_df = pd.DataFrame(dict(Subject=subj, cluster=lab + 1, global_mean=gm_mean.to_numpy(),
                          PC1=X[:, 0], PC2=X[:, 1]))
for c in range(best_k):
    m = cl_df["cluster"] == c + 1
    L.append("簇 %d: n=%d 全局均值 %.2f±%.2f | %s"
             % (c + 1, m.sum(), cl_df.loc[m, "global_mean"].mean(),
                cl_df.loc[m, "global_mean"].std(ddof=1),
                ", ".join(cl_df.loc[m, "Subject"])))
if best_k >= 2:
    from scipy.stats import ttest_ind, mannwhitneyu
    a = cl_df.loc[cl_df["cluster"] == 1, "global_mean"]
    b = cl_df.loc[cl_df["cluster"] == 2, "global_mean"]
    t, p = ttest_ind(a, b, equal_var=False)
    pu = mannwhitneyu(a, b).pvalue
    L.append("簇1 vs 簇2 全局均值: Welch t=%.2f p=%.3f ; Mann-Whitney p=%.3f"
             % (t, p, pu))
    L.append("!! 循环论证警告: 簇本身就是用同一批强度变量算出来的，组间必然显著，")
    L.append("   此处 p 值不构成独立证据。真正的独立判据是 silhouette = %.3f (<0.5 无结构)"
             % float(best["silhouette"]))
cl_df.to_csv(os.path.join(OUT, "assessor_clusters.csv"), index=False, encoding="utf-8-sig")

L.append("\n" + "=" * 78)
L.append("【结论】")
L.append("=" * 78)
L.append("1. 品评员整体不呈双峰：GMM 最优 k=%d，ΔBIC(BIC1−BIC2)=%+.2f（<=6 即弱/不支持）" % (kb, d2))
L.append("2. silhouette 全部 %.2f–%.2f（<0.5）→ 敏感者/不敏感者分群无法从本数据证实"
         % (S["silhouette"].min(), S["silhouette"].max()))
L.append("3. 单点双峰: 87 点中 %d 点 dBIC>10，但经 5 条严格判据（权重/量程/空档）后"
         "仅 %d 点 (%.0f%%) 成立" % (n_all, n_strict, 100 * n_strict / max(len(BD), 1)))
L.append("4. 顶部条目（二甲基三硫C5、异丁酸乙酯C8 等）观测值多为 12 顶格 → 量表上限"
         "截断造成的伪双峰，非感知异质性")
L.append("5. 低浓度端同理：0 值膨胀（地板效应）也会被 GMM 误判为双峰")
L.append("6. 若要真正证明敏感者/不敏感者双峰：需在中间浓度补梯度、改用无上限量表"
         "（0–100 VAS）或 3-AFC 检出/未检出二值设计 + 更大样本(n>=30)")
L.append("7. 现有数据只能支持『强度响应连续、个体差异存在』，不支持『存在两类人群』")

open(os.path.join(OUT, "QC_聚类与双峰检验.txt"), "w", encoding="utf-8").write("\n".join(L))
open(r"C:\Users\24252\AppData\Local\Temp\opencode\cluster_log2.txt", "w",
     encoding="utf-8").write("\n".join(L))
print("done")
