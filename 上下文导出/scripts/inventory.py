# -*- coding: utf-8 -*-
"""生成文件清单（真实大小/时间/用途）"""
import os, hashlib, datetime, csv

BASE = r"C:\Users\24252\Desktop\PPT\课题\26年\26.10\英文"
OUTDIR = os.path.join(BASE, "上下文导出")
os.makedirs(OUTDIR, exist_ok=True)

PURPOSE = {
    "单体混合感官特征理化数据.xlsx": ("原始数据", "主表：组合数据(156x43) / Sheet5 浓度阶梯(91x2) / 样品设计(65x26)", "只读"),
    "感官数据_混合化合物.xlsx": ("原始数据", "单体原始 1674 行 / 单体平均 96 / 混合平均 65(强度全空) / 混合原始 1261(无强度)", "只读"),
    "重构分析.xlsx": ("原始数据", "含第三套阈值 Threshold 列", "只读（曾存在 ~$ 锁文件）"),
    "data-analysis.xlsx": ("原始数据", "29 化合物 x 170 列文献浓度来源", "只读"),
    "杨泽镔-毕业论文_20250506000923.pdf": ("论文原文", "单体香气化合物感知特性", "只读"),
    "华铖凯-毕业论文.pdf": ("论文原文", "混合香气化合物感知特性", "只读"),
    "教程_浓度强度剂量响应分析.md": ("教程", "R值机制/模型选择/软件对比/15篇文献", "交付"),
    "教程_浓度响应图制作.md": ("教程", "7图层拆解 + Prism/R 流程 + 图注模板", "交付"),
    "英文图注与Methods.md": ("交付文档", "英文 Figure Legends + Methods 2.1-2.8 + 数字溯源表 + 核对清单", "交付"),
}

rows = []
for root, dirs, files in os.walk(BASE):
    dirs[:] = [d for d in dirs if d not in ("__pycache__",)]
    for fn in files:
        p = os.path.join(root, fn)
        rel = os.path.relpath(p, BASE)
        try:
            st = os.stat(p)
        except OSError:
            continue
        try:
            h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
        except OSError:
            h = ""
        cat, desc, mode = PURPOSE.get(fn, ("产出/中间文件", "", ""))
        if rel.startswith("教学示范图"):
            cat, desc, mode = "重分析产出", "图/表/QC报告", "新增"
        elif rel.startswith("上下文导出"):
            if fn == "文件清单.csv":
                continue
            cat = "上下文导出"
            desc = {"context.json": "机读版上下文快照",
                    "00_会话上下文.md": "人类可读的完整交接文档",
                    "README.md": "本包说明与阅读顺序"}.get(fn, "复现脚本" if rel.startswith("上下文导出\\scripts") else "")
            mode = "新增"
        rows.append(dict(相对路径=rel, 类别=cat, 用途=desc, 读写=mode,
                         大小KB=round(st.st_size / 1024, 1),
                         修改时间=datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
                         SHA256前16=h))

rows.sort(key=lambda r: (r["类别"], r["相对路径"]))
csvp = os.path.join(OUTDIR, "文件清单.csv")
with open(csvp, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print("files=%d" % len(rows))
for r in rows:
    print("%-58s %8.1f KB  %s  [%s]" % (r["相对路径"], r["大小KB"], r["类别"], r["读写"]))
