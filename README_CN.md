# PhyloSlicer

[![tests](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/tests.yml/badge.svg)](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/tests.yml)
[![CodeQL](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/codeql.yml/badge.svg)](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/codeql.yml)
[![build](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/build.yml/badge.svg)](https://github.com/ZengZichao/PhyloSlicer/actions/workflows/build.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23057620.svg)](https://doi.org/10.5281/zenodo.23057620)

**面向宏生态学的系统发育树时间切片与多样性累积速率分析工具（Python）**

版本 0.1.1 · BSD-3-Clause · Python ≥ 3.10

仓库地址：<https://github.com/ZengZichao/PhyloSlicer>

---

## 简介

PhyloSlicer 独立实现了 R 包 **treesliceR** 开创的"切片系统发育树以推断进化
时间维度上系统发育格局"范式（Araujo et al. 2025, *Ecography* 2025, e07364,
<https://doi.org/10.1111/ecog.07364>），并扩展了后验树不确定性量化、稳健拟合
诊断、空间工作流工具和命令行界面。

"快"是结构性的结果，不是宣传用语。切片只是一张 `E × k` 的区间相交贡献矩阵，
而参考设计要为 `k` 个切片各复制一棵完整的树：唯一随 `k` 增长的只有这一个数组，
那 `k` 份整树复制则完全不存在。Python 侧峰值内存已实测
（`benchmarks/benchmark_memory.py` → `benchmarks/results_memory_py.csv`）：
在 2000 尖处，单张 `(E, k)` 数组在 k = 25 / 50 / 100 时峰值为
1.93 / 3.49 / 6.69 MB，而 `k` 份数组级整树复制为 2.81 / 5.62 / 11.23 MB，
节省 1.45～1.68 倍且随 `k` 扩大。由于本包的矩阵本身也随 `k` 增长，诚实的表述
针对的是**对象数量与复制结构**，而不是一个会发散的内存比值；treesliceR 的 R 侧
峰值**未实测**，文档同样保留这一限定。

### 主要特性

1. **基于数组的切片** — 切片表示为区间相交贡献矩阵
   （`C[e, j] = max(0, min(d_e, t_j) − max(b_e, t_{j−1}))`），不复制任何树。
2. **稀疏矩阵指数** — 通过单次稀疏矩阵乘 `A @ C` 计算逐样区逐切片
   PD/PE/PB(+RW)。
3. **稳健速率拟合** — 有界多起点优化，每次拟合返回收敛状态、R²、残差和失败原因。
4. **后验树不确定性** — `posterior_rate` / `posterior_origin` 可信区间和
   `bootstrap_sites` 稳定性分析。
5. **空间工作流** — rook/queen/knn/radius 邻接构建器和接受 GeoDataFrame 的速率
   地图（geopandas 可选）。
6. **命令行界面** — `phyloslicer slice | rates | posterior | validate`，附带 JSON
   运行摘要。

## 文档

- **英文使用手册（正文）：** [docs/USER_GUIDE.md](docs/USER_GUIDE.md)
- **中文使用手册：** [docs/USER_GUIDE_CN.md](docs/USER_GUIDE_CN.md)
- **变更日志：** [CHANGELOG_CN.md](CHANGELOG_CN.md) · 英文：[CHANGELOG.md](CHANGELOG.md)
- **贡献指南：** [CONTRIBUTING_CN.md](CONTRIBUTING_CN.md) · 英文：[CONTRIBUTING.md](CONTRIBUTING.md)
- **发布流程：** [RELEASE_CN.md](RELEASE_CN.md) · 英文：[RELEASE.md](RELEASE.md)
- **CI 到底跑什么（审稿人可见镜像）：** [docs/CI_WORKFLOWS_CN.md](docs/CI_WORKFLOWS_CN.md) · 英文：[docs/CI_WORKFLOWS.md](docs/CI_WORKFLOWS.md)
- **稿件字数统计口径（可复现）：** [docs/WORD_COUNT_CN.md](docs/WORD_COUNT_CN.md) · 英文：[docs/WORD_COUNT.md](docs/WORD_COUNT.md)
- **Zenodo 存档核对清单（仅作者可执行）：** [docs/ZENODO_DEPOSIT_CHECKLIST_CN.md](docs/ZENODO_DEPOSIT_CHECKLIST_CN.md) · 英文：[docs/ZENODO_DEPOSIT_CHECKLIST.md](docs/ZENODO_DEPOSIT_CHECKLIST.md)

每一对文档中，英文版是权威版本，中文版是它的翻译，两者保持同一份章节结构；
`tests/test_release_hygiene.py` 会在两份标题结构出现偏离时让构建失败，因此中文版本
不会悄悄退化成节译。

## 安装

```bash
pip install -e .            # 核心包（numpy、scipy、pandas、matplotlib）
pip install -e ".[io]"      # + DendroPy / Biopython 适配器
pip install -e ".[spatial]" # + GeoPandas
pip install -e ".[accel]"   # + numba 批量加速内核
pip install -e ".[dev]"     # + pytest、ruff、mypy、pre-commit
```

## 快速上手

先运行内置示例——无需任何外部数据、下载或配置，它正是使用手册逐步演示的流程：

```bash
python examples/quickstart.py
```

下面的代码片段是把同一流程用在**你自己的**文件上：

```python
import phyloslicer as ps

tree = ps.TreeArray.from_newick("birds.tre")
mat = ps.io.presence_matrix(long_table, site_col="grid", species_col="species")
adj = ps.spatial.adjacency_from_coords(x, y, method="knn", k=8)

rates = ps.cpd_rate(tree, mat, n_slices=100)      # CpD, PD, pDO + 诊断
post = ps.posterior_rate("posterior_trees/", mat, n_slices=100)
beta = ps.cpb_rate(tree, mat, adj, component="sorensen")
```

命令行：

```bash
phyloslicer rates input.tre sites.csv --rate cpd --n-slices 100 --out rates.csv
phyloslicer slice input.tre --rootward 5.0 --out crown.tre
```

## 测试与验证

```bash
python -m pytest tests/                     # 单元 + 回归测试
python validation/compare.py                # Python 与 treesliceR 等价比对
```

比对集合为 `{PD, PE, CpD, CpE}` 加上 beta 族（`CpB`、`CpB_RW`、
`r_phylo(index = "PB")`）——14 个 β「组分 × 方法」家族中有 11 个已与
treesliceR 1.1.0 做数值交叉核对。完整缺口清单与容差依据见
[`validation/NOTES_CN.md`](validation/NOTES_CN.md)。

## 基准测试

`benchmarks/` 收录四套彼此独立的协议及其存档表格。引用其中数字前请先读
[`benchmarks/README_CN.md`](benchmarks/README_CN.md)：同名列在不同文件里含义不同，
且只有固定 k 的缩放表是跨语言比较。`python benchmarks/recompute_ratios.py`
会从 CSV 重算文档引用的每一个加速比区间，CI 会执行它。

## 引用

使用 PhyloSlicer 时请同时引用本软件与方法学来源。本包的引用信息见
[`CITATION.cff`](CITATION.cff)；全部版本通用的 Zenodo 概念 DOI 为
[10.5281/zenodo.23057620](https://doi.org/10.5281/zenodo.23057620)。

> Araujo, M.L., Ferreira, L.G.S.S., Nakamura, G., Coelho, M.T.P., Rangel, T.F.
> (2025) 'treesliceR': a package for slicing phylogenies and inferring
> phylogenetic patterns over evolutionary time. *Ecography* 2025, e07364.
> https://doi.org/10.1111/ecog.07364

## 许可

代码：BSD-3-Clause（见 [LICENSE](LICENSE)）。文档：CC-BY 4.0。
第三方署名，包括 `treesliceR` 参考值的来源说明：[NOTICE_CN.md](NOTICE_CN.md) ·
英文：[NOTICE.md](NOTICE.md)。

### 为什么是 BSD-3-Clause

*Ecography* 的 Software Notes 页面推荐"完全开放的许可证（例如 MIT）"，同时说明
copyleft 类许可证（例如 GPL）也接受。BSD-3-Clause 不是其中任一范例，但它确实在
推荐范围之内：OSI 认证、宽松（permissive）、不含 copyleft 条件，BSD-3-Clause 的
作品可以并入 MIT 或 GPL 项目而不产生再许可纠纷。它也是本包所依赖的科学计算栈的
许可证——numpy、scipy 与 pandas 均为 BSD-3-Clause；用于交叉验证的参考实现
treesliceR 是 MIT，同样宽松、同样兼容。与 MIT 唯一的实质差别是第三条：它要求使用者
不得以版权持有人名义为衍生作品背书。这是有意选择而非障碍，也正是部分期刊更偏爱
MIT 的原因。此处不存在任何使用限制：本软件自由且开源，随附的 `LICENSE` 是
BSD-3-Clause 原文。

### 开放获取

*Ecography* 为完全开放获取期刊：每篇文章都在知识共享（Creative Commons）许可证下
发表，所有作者都需向北欧生态学会（Nordic Society Oikos）授予出版权，并且需缴纳
文章出版费（APC）。本仓库不对金额、也不对任何人的减免资格作任何主张：通讯作者应在
录用、签署 Wiley 作者许可协议（WALS）时向编辑部核实当期费用与本机构的协议情况。
减免与折扣确实存在（Wiley 面向符合条件国家的 Research4Life 计划，以及 NSO 资助的
豁免），因此没有出版经费的团队应当主动询问，而不是默认这笔费用无法避免。请注意
边界划分：知识共享许可证管的是*文章*本身；软件仍然是 BSD-3-Clause，本仓库的文字
仍然是 CC-BY 4.0（见 [NOTICE_CN.md](NOTICE_CN.md)）。

## 作者

曾子超（Zichao Zeng）— 上海交通大学生命科学技术学院 ·
<zengzichao@sjtu.edu.cn> ·
ORCID [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)
