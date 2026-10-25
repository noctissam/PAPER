# -*- coding: utf-8 -*-
import os, pandas as pd, numpy as np
BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
L = []
f = os.path.join(BASE, "感官数据_混合化合物.xlsx")
xl = pd.ExcelFile(f)
L.append("感官数据_混合化合物.xlsx sheets: %s" % str(xl.sheet_names))
for s in xl.sheet_names:
    d = pd.read_excel(f, sheet_name=s, nrows=4)
    L.append("\n[%s] shape(cols)=%d 列名: %s" % (s, d.shape[1], str(list(d.columns)[:30])))

f2 = os.path.join(BASE, "单体混合感官特征理化数据.xlsx")
xl2 = pd.ExcelFile(f2)
L.append("\n\n单体混合感官特征理化数据.xlsx sheets: %s" % str(xl2.sheet_names))
for s in xl2.sheet_names:
    d = pd.read_excel(f2, sheet_name=s, nrows=4)
    L.append("\n[%s] 列数=%d 列名: %s" % (s, d.shape[1], str(list(d.columns)[:45])))
    if d.shape[1] <= 5:
        L.append("   前几行: %s" % str(d.head(4).values.tolist()))

# 找含浓度的表
for path in [f, f2, os.path.join(BASE, "重构分析.xlsx"), os.path.join(BASE, "data-analysis.xlsx")]:
    if not os.path.exists(path): continue
    xl = pd.ExcelFile(path)
    for s in xl.sheet_names:
        d = pd.read_excel(path, sheet_name=s, nrows=60)
        cols = [str(c) for c in d.columns]
        hit = [c for c in cols if any(k in c.lower() for k in
               ["浓度", "conc", "ppb", "μg", "ug/l", "mg/l", "阈值", "threshold", "content"])]
        if hit:
            L.append("\n*** 含浓度字段 %s :: %s -> %s" % (os.path.basename(path), s, str(hit)))

open(r"C:\Users\24252\AppData\Local\Temp\opencode\conc.txt", "w", encoding="utf-8").write("\n".join(L))
print("ok")
