# -*- coding: utf-8 -*-
import os, re, pandas as pd, numpy as np
BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
L=[]
def norm(s):
    if s is None: return ""
    s=re.sub(r"\s+","",str(s))
    for a,b in [("（","("),("）",")")]: s=s.replace(a,b)
    return s

# sklearn 可用性
try:
    import sklearn; L.append("sklearn = %s" % sklearn.__version__)
except Exception as e: L.append("sklearn 不可用: %s" % e)
try:
    import scipy; L.append("scipy = %s" % scipy.__version__)
except Exception as e: L.append("scipy 不可用")

raw = pd.read_excel(os.path.join(BASE,"感官数据_混合化合物.xlsx"), sheet_name="单体数据原始数据")
cols = list(raw.columns)
L.append("列: %s" % str(cols))
raw = raw.rename(columns={cols[2]:"I"})
raw["key"] = raw["Sample"].ffill()
raw["Subject"] = raw["Subject"].astype(str)
raw["I"] = pd.to_numeric(raw["I"], errors="coerce")

# 苯乙醇全部行
L.append("\n### 苯乙醇 相关行 (key 含 苯乙醇)")
sel = raw[raw["key"].astype(str).str.contains("苯乙醇", na=False)]
L.append("  行数 = %d" % len(sel))
L.append("  key 取值 = %s" % str(sel["key"].unique().tolist()))
L.append("  Subject 前缀分布: %s" % str(sel["Subject"].value_counts().head(8).to_dict()))
L.append("\n  逐 level 分块看（前 3 个 level）:")
for lv in sorted(sel["key"].unique())[:3]:
    d = sel[sel["key"]==lv]
    L.append("   [%s] n=%d" % (lv, len(d)))
    L.append("      序列1: %s" % str(d["I"].head(20).tolist()))
    L.append("      序列2: %s" % str(d["I"].tail(20).tolist()))
    L.append("      Subject: %s" % str(d["Subject"].tolist()))
    # 是否前18 = 后18
    a = d["I"].head(18).to_numpy(); b = d["I"].tail(18).to_numpy()
    if len(d)==36:
        L.append("      前18==后18 ? %s ; 相关系数=%.3f ; 前18均值=%.3f 后18均值=%.3f"
                 % (bool(np.array_equal(a,b)), np.corrcoef(a,b)[0,1], np.nanmean(a), np.nanmean(b)))

# 全表：每个 key 的行数分布
L.append("\n### 全表每个 key 的记录数分布")
cnt = raw[raw["key"].ne("") & raw["I"].notna()].groupby("key").size()
L.append("  分布: %s" % str(cnt.value_counts().to_dict()))
L.append("  非18的: %s" % str(cnt[cnt!=18].to_dict()))

# Subject 命名规则
L.append("\n### Subject 取值样例: %s" % str(raw["Subject"].dropna().unique()[:30].tolist()))

# 平均值表里的苯乙醇
m = pd.read_excel(os.path.join(BASE,"感官数据_混合化合物.xlsx"), sheet_name="单体数据平均值")
L.append("\n### 单体数据平均值 shape=%s 列=%s" % (str(m.shape), str(list(m.columns))))
mm = m[m.iloc[:,0].astype(str).str.contains("苯乙醇", na=False)]
L.append("  苯乙醇行数 = %d" % len(mm))
for i in range(len(mm)):
    L.append("    " + " | ".join(str(v) for v in mm.iloc[i,:4].tolist()))

# 强度超出 0-12 的
L.append("\n### 原始强度越界检查")
L.append("  I 范围 = %.3f ~ %.3f ; 越界条数 = %d" % (raw["I"].min(), raw["I"].max(),
         int(((raw["I"]<0)|(raw["I"]>12)).sum())))
bad = raw[(raw["I"]<0)|(raw["I"]>12)]
if len(bad):
    L.append("  越界明细: %s" % str(bad[["key","Subject","I"]].head(20).values.tolist()))

open(r"C:\Users\24252\AppData\Local\Temp\opencode\dup.txt","w",encoding="utf-8").write("\n".join(L))
print("ok")
