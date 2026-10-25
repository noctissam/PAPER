# -*- coding: utf-8 -*-
"""A-1: 建立干净数据层 + 数据质量报告"""
import os, re, math, json, numpy as np, pandas as pd
BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUT  = os.path.join(BASE, "教学示范图")
os.makedirs(OUT, exist_ok=True)
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

SCALE_MIN, SCALE_MAX = 0.0, 12.0

raw = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"),
                    sheet_name="单体数据原始数据")
raw = raw.rename(columns={raw.columns[2]: "I"})
raw["key"] = raw["Sample"].ffill()
raw["Subject"] = raw["Subject"].astype(str).str.strip()
raw["I"] = pd.to_numeric(raw["I"], errors="coerce")
raw = raw[raw["key"].notna() & raw["Subject"].str.match(r"^S\d+$", na=False)].copy()
raw["cmpd"] = raw["key"].map(base_name)
raw["lvl"] = raw["key"].map(level)
raw = raw[raw["cmpd"].ne("") & raw["I"].notna()].copy()
raw["key"] = raw["key"].map(norm)
raw["key"] = raw.apply(
    lambda r: "%s(C%d)" % (ALIAS.get(r["cmpd"], r["cmpd"]), r["lvl"])
    if r["lvl"] is not None else r["key"], axis=1)
raw["I_raw"] = raw["I"]

L.append("=" * 78)
L.append("【1】原始数据层")
L.append("=" * 78)
L.append("有效记录 %d 条, 化合物 %d, key 数 %d, 品评员 %d"
         % (len(raw), raw["cmpd"].nunique(), raw["key"].nunique(), raw["Subject"].nunique()))

# --- 越界 ---
oob = raw[(raw["I"] < SCALE_MIN) | (raw["I"] > SCALE_MAX)].copy()
L.append("\n越界记录 %d 条:" % len(oob))
for _, r in oob.iterrows():
    L.append("   %s / %s = %.1f  (量表 0–12)" % (r["key"], r["Subject"], r["I"]))
raw["oob"] = (raw["I"] < SCALE_MIN) | (raw["I"] > SCALE_MAX)
raw["I_clip"] = raw["I"].clip(SCALE_MIN, SCALE_MAX)

# --- 重复检测 ---
dup_mask = raw.duplicated(subset=["key", "Subject"], keep=False)
dup_keys = sorted(raw.loc[dup_mask, "key"].unique())
L.append("\n重复 (key, Subject) 记录 %d 条，涉及 %d 个 key:" % (dup_mask.sum(), len(dup_keys)))
for k in dup_keys:
    L.append("   %s  n=%d" % (k, (raw["key"] == k).sum()))

# --- 两轮重测信度（苯乙醇）---
L.append("\n" + "=" * 78)
L.append("【2】第二轮数据（被平均值表丢弃的 18×6 = 108 条）")
L.append("=" * 78)
rep_keys = [k for k in dup_keys if (raw["key"] == k).sum() == 36]
testcorr = {}
for k in rep_keys:
    d = raw[raw["key"] == k].reset_index(drop=True)
    r1 = d["I_clip"].iloc[:18].to_numpy()
    r2 = d["I_clip"].iloc[18:].to_numpy()
    if len(d) == 36 and np.std(r1) > 0 and np.std(r2) > 0:
        testcorr[k] = float(np.corrcoef(r1, r2)[0, 1])
L.append("两轮重测相关 (Pearson r):")
for k, v in sorted(testcorr.items(), key=lambda x: x[1]):
    L.append("   %-22s r = %+.3f" % (k, v))
if testcorr:
    L.append("   均值 r = %+.3f  (理想 >0.7；说明第二轮数据是否可用)"
             % np.mean(list(testcorr.values())))

# 两轮均值对比 vs 平均值表
m = pd.read_excel(os.path.join(BASE, "感官数据_混合化合物.xlsx"),
                  sheet_name="单体数据平均值")
m["key"] = m["Sample"].map(norm).map(
    lambda k: "%s(C%d)" % (ALIAS.get(base_name(k), base_name(k)), level(k))
    if level(k) is not None else k)
m["I_mean_tab"] = pd.to_numeric(m["气味强度"], errors="coerce")
tab = dict(zip(m["key"], m["I_mean_tab"]))
L.append("\n三套均值对比 (平均值表 / 第一轮 / 第二轮 / 36条合并):")
L.append("   %-22s %8s %8s %8s %8s %8s" % ("key", "平均表", "第1轮", "第2轮", "36合并", "差值"))
diff_rows = []
for k in rep_keys:
    d = raw[raw["key"] == k].reset_index(drop=True)
    r1 = d["I_clip"].iloc[:18].mean(); r2 = d["I_clip"].iloc[18:].mean()
    a36 = d["I_clip"].mean(); t = tab.get(k, np.nan)
    L.append("   %-22s %8.3f %8.3f %8.3f %8.3f %+8.3f" % (k, t, r1, r2, a36, a36 - t))
    diff_rows.append(dict(key=k, tab=t, round1=r1, round2=r2, all36=a36, delta=a36 - t))

# --- 按 Subject 聚合（统计正确：每人 1 个独立单位）---
# 关键：苯乙醇第二轮重测 r≈0（见上），不可与第一轮合并 → 主分析只用第一轮
raw["occ"] = raw.groupby(["key", "Subject"]).cumcount()
primary = raw[raw["occ"] == 0].copy()
agg = (primary.groupby(["cmpd", "lvl", "key", "Subject"], as_index=False)
         .agg(I=("I_clip", "mean"), n_rep=("I_clip", "size"),
              oob=("oob", "any")))
agg_all = (raw.groupby(["cmpd", "lvl", "key", "Subject"], as_index=False)
             .agg(I=("I_clip", "mean"), n_rep=("I_clip", "size"), oob=("oob", "any")))
L.append("\n" + "=" * 78)
L.append("【3】聚合后（每人 1 值 → 18 个独立单位）")
L.append("=" * 78)
cc = agg.groupby("key").size()
L.append("key 数 %d；每个 key 的人数分布: %s" % (len(cc), str(cc.value_counts().to_dict())))
L.append("总记录 %d 条（原 %d 条 → 去掉重复轮次 %d 条）"
         % (len(agg), len(raw), len(raw) - len(agg)))

# --- 与平均值表的全面对比 ---
L.append("\n" + "=" * 78)
L.append("【4】聚合均值 vs 论文《单体数据平均值》逐格核对")
L.append("=" * 78)
g = agg.groupby("key")["I"].mean()
cmp_rows = []
for k, v in g.items():
    t = tab.get(k, np.nan)
    cmp_rows.append(dict(key=k, mine=v, table=t,
                         delta=(v - t) if pd.notna(t) else np.nan))
C = pd.DataFrame(cmp_rows)
have = C[C["table"].notna()]
L.append("可比对 key = %d / %d" % (len(have), len(C)))
L.append("零误差 (|Δ|<1e-6) = %d" % int((have["delta"].abs() < 1e-6).sum()))
L.append("有差异 = %d" % int((have["delta"].abs() >= 1e-6).sum()))
L.append("\n差异明细:")
for _, r in have[have["delta"].abs() >= 1e-6].sort_values("delta", key=abs, ascending=False).iterrows():
    L.append("   %-22s 我=%7.3f  表=%7.3f  Δ=%+7.3f" % (r["key"], r["mine"], r["table"], r["delta"]))
only_tab = sorted(set(tab) - set(g.index))
only_mine = sorted(set(g.index) - set(tab))
L.append("\n平均值表有、原始数据无的 key (%d):" % len(only_tab))
L.append("   " + " | ".join(only_tab))
L.append("\n原始数据有、平均值表无的 key (%d):" % len(only_mine))
L.append("   " + " | ".join(only_mine))

# --- 保存 ---
agg.to_csv(os.path.join(OUT, "clean_assessor_level.csv"), index=False, encoding="utf-8-sig")
agg_all.to_csv(os.path.join(OUT, "clean_assessor_level_两轮合并.csv"),
               index=False, encoding="utf-8-sig")
C.to_csv(os.path.join(OUT, "qc_mean_vs_table.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(diff_rows).to_csv(os.path.join(OUT, "qc_round1_vs_round2.csv"),
                               index=False, encoding="utf-8-sig")
raw.to_csv(os.path.join(OUT, "clean_raw_flagged.csv"), index=False, encoding="utf-8-sig")

open(os.path.join(OUT, "QC_数据质量报告.txt"), "w", encoding="utf-8").write("\n".join(L))
open(r"C:\Users\24252\AppData\Local\Temp\opencode\clean_log.txt", "w",
     encoding="utf-8").write("\n".join(L))
print("done")
