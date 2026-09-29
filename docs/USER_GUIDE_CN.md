# PhyloSlicer — 用户手册

**面向宏生态学的系统发育时间切片与累积多样性速率分析（Python）**

这里的"快"是结构性的结果，而不是宣传用语：切片只是一张 `E x k` 的区间相交贡献矩阵，而不是参考设计所物化出来的那 `k` 份完整树副本，因此运行时间与内存都不会被切片数量放大。这句话所支撑的就是这一主张；实测证据位于 `benchmarks/`，其中的跨语言加速比是从固定 k 的缩放表里读出的（`results_scaling_py.csv` 对 `results_scaling_r.csv`，仅含切片内核，输入树字节级相同，每个点取 15 次计时重复的中位数），而 `results.csv` 比较的是本包与*它自己的*朴素实现，因此它*不是*跨语言基准。此处不引用任何计时数字：R 的那几列是归档的基准输出，而笔记本电脑上的 R 计时在不同会话之间会浮动——参考包自己的说明就警告过两倍的幅度。三个 R 脚本已在 R 4.5.3（aarch64-apple-darwin20.0.0）配合 ape 5.8.1 与 treesliceR 1.1.0 的环境下重新运行到临时文件中，把归档中位数复现到 6 % 以内（验证场景上的中位数比值 1.06，1000/2000 尖场景上 1.02-1.08），但已提交的 CSV 保持原样，好让文档引用的数字正是当初测量它们时所对照的那一组。

版本 0.1.0  ·  BSD-3-Clause  ·  Python ≥ 3.10

PhyloSlicer 是对**"切片系统发育树以推断进化时间维度上的系统发育格局"**这一范式的独立、向量化 Python 实现，该范式由 **treesliceR**（Araujo et al. 2025, *Ecography* 2025, e07364, doi:10.1111/ecog.07364；网络首发于 2024 年 10 月 28 日）提出；本包在此之上扩展了后验树不确定性量化、稳健拟合诊断、空间工作流辅助工具和命令行界面。英文文档：[USER_GUIDE.md](USER_GUIDE.md)，即本手册的记录版本/正文。

---

## 目录

1. [安装](#1-安装)
2. [五分钟快速上手](#2-五分钟快速上手)
3. [核心概念与约定](#3-核心概念与约定)
4. [输入输出：树与矩阵](#4-输入输出树与矩阵)
5. [切片](#5-切片)
6. [多样性指数（PD、PE、PB）](#6-多样性指数pdpepb)
7. [累积速率（CpD、CpE、CpB）](#7-累积速率cpdcpecpb)
8. [跨后验树的不确定性](#8-跨后验树的不确定性)
9. [空间辅助工具](#9-空间辅助工具)
10. [可视化](#10-可视化)
11. [模拟与内置数据](#11-模拟与内置数据)
12. [命令行界面](#12-命令行界面)
13. [从 treesliceR（R）迁移](#13-从-treeslicerr迁移)
14. [诊断信息、状态与报错说明](#14-诊断信息状态与报错说明)
15. [数值验证](#15-数值验证)
16. [常见问题与故障排查](#16-常见问题与故障排查)
17. [API 速查表](#17-api-速查表)

---

## 1. 安装

```bash
# core package (numpy, scipy, pandas, matplotlib)
pip install -e .

# optional extras
pip install -e ".[io]"       # + DendroPy / Biopython adapters
pip install -e ".[spatial]"  # + GeoPandas
pip install -e ".[accel]"    # + numba (batched fitting kernels)
pip install -e ".[dev]"      # + pytest, ruff, mypy, ...
```

不安装任何可选扩展，全部功能依然可用：缺失的可选依赖会抛出清晰的 `ImportError` 提示（空间网格场景则自动回退到矩形边界算法）。

## 2. 五分钟快速上手

```python
import phyloslicer as ps

# 1. read a dated (ultrametric) Newick tree
tree = ps.TreeArray.from_newick("birds.tre")

# 2. sites x species presence/absence matrix (DataFrame; see section 4)
mat = ps.io.presence_matrix(long_table, site_col="grid_cell", species_col="species")

# 3. slice the phylogeny into 100 equal-time slices and compute CpD rates
rates = ps.cpd_rate(tree, mat, n_slices=100)
print(rates.head())
#    site_id       CpD        PD       pDO  converged  r_squared  reason
# 0      s_1  1.234...  45.6...   2.42...       True     0.98...

# 4. phylogenetic beta-diversity rates over queen-adjacent neighbours
adj = ps.spatial.adjacency_from_coords(x, y, method="knn", k=8)
beta = ps.cpb_rate(tree, mat, adj, n_slices=100, component="sorensen")

# 5. uncertainty across a folder of posterior trees
post = ps.posterior_rate("posterior_trees/", mat, n_slices=100)
```

命令行等价操作：

```bash
phyloslicer rates   birds.tre sites.csv --rate cpd --n-slices 100 --out rates.csv
phyloslicer slice   birds.tre --rootward 5.0 --out crown.tre
phyloslicer posterior posterior_trees/ sites.csv --hdi 0.95 --out posterior.csv
phyloslicer validate birds.tre sites.csv --reference r_reference.csv --out report.csv
```

## 3. 核心概念与约定

### 3.1 时间、深度与年龄

PhyloSlicer 假定输入为**定年的超度量树**。内部使用两套坐标：

* **距根深度**（根为 0，枝端为 `root_age`）——区间运算的自然坐标；
* **距今年龄**（`root_age − depth`）——用户可见的坐标（Ma，或你的树所用的任何枝长时间单位）。

切片窗口一律以年龄形式报告，按 根 → 枝端 排序；速率模型使用的 `ages` 向量是**每个切片较老一侧的边界**，单位为**距今年龄（Ma）**（`root_age − depth`），而不是距根深度。非超度量树默认被拒绝：`slice_pieces()` 会抛出 `ValueError`，除非你传入 `allow_non_ultrametric=True`——此时该项检查降级为警告，并在结果表中记录 `ultrametric = False`；软件不对树做任何重标定，此时 `root_age` 取最深的枝端深度。

> 为什么默认是抛错：非超度量树照样能算出数字，但它们不再表示"时间 *t* 之前的多样性"，而一条静默告警很容易在日志里被漏掉。确实想要审计列时，请显式传入 `allow_non_ultrametric=True`。

### 3.2 贡献矩阵

每条边 `e` 在深度坐标上跨越区间 `[b_e, d_e]`，切片则是一个窗口 `[t_{j−1}, t_j]`。边 `e` 对切片 `j` 贡献的枝长为

```
C[e, j] = max(0, min(d_e, t_j) − max(b_e, t_{j−1}))
```

对所有 边 × 切片 一次向量化算出。**全程不复制任何树**：k 个切片的堆栈只是一个 `E × k` 数组，而不是 k 棵树对象（这正是与 treesliceR 的 `phylo_pieces` 的结构性差异——后者每个切片都物化一份整棵系统发育树的修改副本）。

### 3.3 存在/缺失矩阵

所有指数与速率函数都接受 **样区 × 物种** 的 `pandas.DataFrame`：行是集合区（`index` 会在所有输出中成为 `site_id`），列是物种名。取值按 `> 0` 判定。物种名会与树的枝端标签进行匹配；`align_tree_matrix()`（由速率函数自动调用）会剪掉矩阵中缺失的树枝端、丢弃树中缺失的矩阵物种、删除从未被观测到的物种——**每一次调整都发出警告**。

### 3.4 速率模型

对每个集合区，逐切片指数值自根向枝端累加、归一化到和为 1，然后用有界最小二乘拟合 `累积量 ∝ exp(−r · age)`（边界 `[1e-6, 50]`）。每次拟合都返回：

| 字段 | 含义 |
|---|---|
| `rate` | 拟合速率 r（拟合失败时为 NaN） |
| `converged` | 仅当优化器认证了最优解时为 True |
| `r_squared` | 拟合累积曲线的 1 − SSE/SST |
| `residuals` | 逐切片残差（`diagnostics=False` 时为空） |
| `start_used`、`n_iter` | 优化器溯源信息 |
| `reason` | 收敛时为空；否则给出失败的通俗语言解释 |

起源时间按 `pXO = −ln(1 − percentile) / r` 计算；默认的 `percentile=0.95` 与 treesliceR 的默认 `pDO = 5` 相对应：两者都返回 `-ln(0.05)/r`，即拟合曲线取值为 0.05 的那个年龄——参考实现直接命名剩余的那 5 %，本包命名的是累积到的 95 % 补集。

## 4. 输入输出：树与矩阵

### 4.1 树

```python
from phyloslicer import TreeArray
from phyloslicer.io import read_newick, read_nexus, parse_newick, write_newick, write_nexus

tree  = TreeArray.from_newick("tree.tre")     # file path or Newick string
tree  = read_newick("(A:1.0,B:1.0):1.0;")     # string directly
trees = TreeArray.from_nexus("posterior.nex") # list[TreeArray], multi-tree
tree  = parse_newick("((A:1,B:1)95:1,(C:1,D:1)100:1);")  # internal labels read but not stored

tree.to_newick()                # Newick string (root has no stem length)
write_nexus(trees)              # NEXUS string with all trees
tree.to_dendropy()              # requires phyloslicer[io]
tree.to_biopython()
TreeArray.from_dendropy(dnd_tree)
TreeArray.from_biopython(bp_tree)
```

说明与保证：

* 方括号注释（含 NHX 注释）在解析前会被剥离；支持引号名（`'sp one'`）。
* NEXUS `TREES` 块中的 `TREE` 语句**可以跨多行**。
* 构造时枝端标签按深度优先顺序重新编号；每个枝系的枝端构成连续的下标区间，这正是区间内核得以实现的前提。
* 内部（支持度）标签在读写往返中**不会被保留**——PhyloSlicer 的数据模型只存储枝端、边和枝长。如需支持度请保留原始文件。
* 多树 Newick 文件的第一棵树用 `read_newick` 读取；要流式读取全部树请用第 8 节的 `iter_trees()`。

`TreeArray` 的关键成员：

```python
tree.n_tips, tree.n_edges, tree.root, tree.tip_labels
tree.root_age          # deepest tip depth (== tip depth for ultrametric trees)
tree.tip_ages()        # root_age - depth, per tip
tree.is_ultrametric    # tolerance-based check
tree.check_ultrametric()  # warns when violated
tree.total_pd()        # sum of branch lengths
tree.node_depths()     # (n_nodes,) depth from root
tree.tip_intervals()   # (E, 2) [lo, hi) tip-index ranges covered by each edge
tree.subarray(["A", "B"])  # prune to a subset (MRCA-rooted)
tree.with_lengths(v)   # SlicedTree sharing the topology with new lengths
```

### 4.2 存在/缺失矩阵

```python
from phyloslicer.io import presence_matrix, align_tree_matrix, as_matrix

# wide table: sites as rows, species as columns (light validation only)
wide = presence_matrix(wide_df)

# long table: one row per (site, species[, value]) observation
wide = presence_matrix(long_df,
                       site_col="grid", species_col="species",
                       value_col="count")       # value > 0  -> presence
wide = presence_matrix(long_df, site_index_col=True)  # site ids from the index

# manual alignment with explicit reporting (rate functions do this for you)
tree2, mat2 = align_tree_matrix(tree, mat)   # both realigned + warnings
mat_df = as_matrix(numpy_array_or_sparse)    # coerce to DataFrame
```

`align_tree_matrix` 依次执行以下操作（每一项都带警告）：删除全零的物种列、删除树中不存在的矩阵物种、剪掉矩阵中不存在的树枝端、把列顺序重排为（剪枝后）树的枝端顺序。

## 5. 切片

### 5.1 单个切片

```python
from phyloslicer.slicing import slice_rootward, slice_tipward, slice_interval

crown = slice_rootward(tree, 5.0)               # keep the youngest 5 Ma
back  = slice_tipward(tree, 5.0)                # drop the youngest 5 Ma
mid   = slice_interval(tree, start=8.0, stop=5.0)   # the [5, 8] Ma band
older, younger = slice_interval(tree, 8.0, 5.0, invert=True)

pd_slice = slice_rootward(tree, 2.0, criterion="pd")  # remove root-side PD of 2
```

语义（**端点行为**与 treesliceR 一致；等 PD 窗口的摆放位置则有意不同——见 §13.1）：

| 调用 | 保留 | `criterion="time"` | `criterion="pd"` |
|---|---|---|---|
| `slice_rootward(tree, x)` | 冠部切片 | 最年轻的 `x` Ma | 整树去掉根侧 PD `x` |
| `slice_tipward(tree, x)` | 主干 | 除最年轻 `x` Ma 之外的全部 | 根侧 PD `x` |
| `slice_interval(tree, a, b)` | `[b, a]` 环带（Ma 准则下需 `a > b`） | 两个年龄之间 | 两个根侧 PD 深度之间（`a < b`） |

边界情形有明确校验：`time < 0` 报错；`time > root_age` 报错；`slice_rootward(tree, 0)` 返回**空**切片（PD = 0），`slice_tipward(tree, 0)` 返回整树——这一端点行为与参考实现一致。以 PD 度量的阈值经树的累计 PD 曲线（`depth_for_pd`）换算，该曲线统计每个深度上**真正活跃**的分支数。在节点深度互异的严格二叉树上，它复现参考实现的断点；在多叉树上则有意不复现，而且这一差异未经交叉验证（§13.1）。

返回值是 `SlicedTree`：拓扑与父树相同、枝长被截断、零长边予以保留。可用：

```python
sliced.pd              # PD stored in the slice
sliced.tree            # a TreeArray copy
sliced.collapse()      # TreeArray with zero-length internal edges removed
```

所有切片函数都接受 `collapse=True`，直接返回合并后的 `TreeArray`。

### 5.2 多个切片 — `SliceStack`

```python
from phyloslicer.slicing import slice_pieces

stack = slice_pieces(tree, n=100)                     # 100 equal-time slices
stack = slice_pieces(tree, width=0.5)                 # fixed 0.5-Ma width (rounded
                                                      # treesliceR-style to fit)
stack = slice_pieces(tree, n=50, criterion="pd")      # equal-PD slices
```

`n` 与 `width` 互斥（以此取代 R 中容易混淆的 `n` + `method` 组合）。`SliceStack` 提供：

```python
stack.n_slices
stack.windows        # (k, 2) ages [older, younger], root -> tip
stack.ages           # (k,) older boundary of each slice  <- rate-model ages
stack.widths         # (k,) slice widths
stack.contribution()             # (k,) branch length per slice (all tips)
stack.contribution(["sp_1"])     # idem restricted to a tip subset
stack.contribution_matrix()      # (E, k) full edge x slice matrix
stack.to_sliced_trees()          # materialise slices as SlicedTree objects
```

### 5.3 修剪枝端

```python
from phyloslicer.slicing import prune_tips

kept = prune_tips(tree, 0.4)                    # keep tips with terminal branch >= 0.4
kept = prune_tips(tree, 0.4, side="before")     # drop them instead (R method=2)
kept = prune_tips(tree, 0.75, quantiles=True)   # threshold = 75th percentile
lists = prune_tips(tree, np.array([0.2, 0.4]))  # one tree per threshold
```

若某个阈值会使剩下的枝端少于两个，则抛出 `ValueError`，消息中给出该阈值和可用区间。treesliceR 在此返回 `NULL`，在脚本里这个问题会推迟到对 `None` 取属性时才以属性错误暴露；向量调用中的某个元素也不会再静默产生缺失的列表项。

## 6. 多样性指数（PD、PE、PB）

> **命名提示：** `from phyloslicer.indices import pd` 会把 *PD 函数*绑定到名字 `pd` 上，遮蔽惯用的 `import pandas as pd` 别名。在同一脚本中使用 pandas 时，建议改用
> `from phyloslicer import indices` 并调用 `indices.pd(...)`，或以其他别名导入该函数。

```python
from phyloslicer.indices import pd, pe, pb, pd_per_slice, pe_per_slice, pb_per_slice

faith  = pd(tree, mat)        # (n_sites,) total PD per site
endo   = pe(tree, mat)        # (n_sites,) phylogenetic endemism (range-weighted)
beta   = pb(tree, mat, adj)   # DataFrame: site_id, PB   (full-tree beta)

stack = slice_pieces(tree, n=100)
prof_pd = pd_per_slice(tree, mat, stack)   # (n_sites, k)
prof_pe = pe_per_slice(tree, mat, stack)   # (n_sites, k)
prof_pb = pb_per_slice(tree, mat, adj, stack,
                       component="sorensen",     # "sorensen"|"turnover"|"nestedness"
                       approach="multisite",     # "multisite"|"pairwise"
                       weighted=False,           # True -> PB_RW (Laffan et al. 2018)
                       multisite_domain="mixed") # "paired" removes the drift, §6.1
# prof_pb["values"]  (n_focal, k) per-slice values, summing to the total
# prof_pb["total"]   (n_focal,) full-tree beta of the neighbourhood
# prof_pb["status"]  per-focal status string (section 14)
```

实现要点：

* 样区对边的隶属关系由枝端区间表构建，每个样区一次二分查找——没有 `(样区 × 枝端)` 中间矩阵。
* 所有逐样区、逐切片的量都是一次稀疏矩阵乘 `A @ C`。
* 特有性的范围权重为 `1 /（覆盖该边的样区数）`，与 treesliceR 及 Laffan et al. (2018) 一致。
* 每条 β 序列都**划分邻域的总 β**：对每一个有限行，逐切片值之和为 1。这是已实现的行为，不是假设——`pb_per_slice` 把每条逐切片序列除以**它自身**各份额之和，返回的 `normalisation` 列记录每一行走的是哪条路径：当切片的份额本来就把全树指数加总干净时为 `"full_tree"`（`sorensen` 属于此类，因为它的 `b + c` 分子对切片可加；以及两个分母在数值上恰好重合的邻域），当分母必须由 `sum_j series_j` 重建时为 `"profile"`。
  对 turnover 和 nestedness 而言，该性质是*新近*才成立的：`min(b, c)` 跨切片不可加（`sum_j min(b_j, c_j) <= min(sum_j b_j, sum_j c_j)`），所以若除以全树分母——A2 修复之前的代码如此，treesliceR 至今仍如此——在下面这份工作负载上，turnover 行的行和只会落在 0.96-0.98，nestedness 落在 1.05-1.40。
  支撑这一条的断言是
  `tests/test_regressions.py::TestBetaNormalisation::test_all_beta_profiles_partition_the_neighbourhood_beta`：
  用 `bd20_random` fixture（已提交的交叉验证输入 `bd20__tree.tre` + `bd20__random__mat.csv`，20 尖 / 30 样区，queen 网格邻接）以 `slice_pieces(tree, n=30, criterion="time")` 切片，覆盖全部六个 分量 × 方法 组合，用
  `np.allclose(row_sums, 1.0, atol=1e-8)` 以及更严的
  `max |row_sum - 1| < 1e-12` 校验。它替换了此前一个在 `atol=0.05` 下通过的 6 样区 / 5 切片 fixture。
  另有两条配套断言让这一说法保持诚实：
  `test_profile_and_total_share_the_same_quantities` 以 1e-12 钉住 `total` 仍等于 `pb()`
  （被归一化的是剖面，报告的指数未被重新缩放），
  `test_r_compatible_profile_keeps_the_unnormalised_curve` 则钉住
  `indices.r_compatible_profile` 仍返回 treesliceR 那条*不*作划分的曲线。


### 6.1 multisite 指数取决于邻域大小（`multisite_domain`）

treesliceR 的 multisite β 混用了两个求和域（这是继承自参考实现的，不是移植错误）：

```
a_tot  = Σ_(m sites) PD_i − PD_union        grows ~ linearly in m
bc_tot = Σ_(C(m,2) pairs) (b + c)_pair      grows ~ quadratically in m
β      = bc_tot / (2·a_tot + bc_tot)
```

共有项 `a` 按 *m* 个样区求和，而更替项 `b + c`（以及 `min(b, c)`）按 *C(m, 2)* 个数对求和，因此该比值**不是组成的不变量**：即使每一个数对都没有变化，它仍随邻域大小单调上升。在参考构造下——一根共有骨架枝，加上每个样区一条私有末端枝，长度全为 1，邻域 = 全部 *m* 个样区——Sørensen 值的解析解为 `m / (m + 2)`，默认的 `multisite_domain="mixed"` 精确复现它：

| m | `mixed`（默认，= treesliceR） | 解析值 `m/(m+2)` | 任意单个数对的比较 | `paired` |
|---|---|---|---|---|
| 2 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| 3 | 0.6000 | 0.6000 | 0.5000 | 0.5000 |
| 5 | 0.7143 | 0.7143 | 0.5000 | 0.5000 |
| 10 | 0.8333 | 0.8333 | 0.5000 | 0.5000 |
| 20 | 0.9091 | 0.9091 | 0.5000 | 0.5000 |

**因此：切勿比较在不同规模邻域上算出的 multisite 值**（k-NN 邻接中 `k` 不同，或焦点格邻居数不同的网格）。请把它们读作*按数对求和的 γ 分解*，而不是与规模无关的 β。

* `multisite_domain="mixed"`（默认；`"reference"` 是其可接受的同义词，见
  `indices.MULTISITE_DOMAIN_ALIASES`）复现参考实现的算术。它保持为默认值，是因为中途改换它会悄悄重新定义一个已发表的指标。
  请注意该警示的适用范围：它针对的是**指数**（`PB`），而不是所有下游数字。`pb_per_slice()` 把每行除以本行之和，所以对 `sorensen` 与 `turnover` 而言，multisite 剖面——以及 `cpb_rate()` 由此拟合出的 `CpB`——在两个域下*完全相同*，只有 `PB` 会变。`nestedness` 是速率会随求和域改变的唯一分量，因为它的两项分母不同。
* `multisite_domain="paired"` 把共有项也放进与 `b + c` 和 `min(b, c)` 相同的按数对求和域，从而消除漂移：上表每档 *m* 的值都回到 0.5000，并与同一邻域的 `approach="pairwise"` 聚合一致，同时 `m = 2` 不受影响（两个域在该点重合）。这是对 treesliceR 的**有意偏离**：使用 `"paired"` 之后，multisite 输出不再与 R 结果互相可比。

请注意，跟随参考算术并不等于满足某个已发表的定义。treesliceR 背后的方法学论文（Araujo et al. 2025, *Ecography*, e07364, doi:10.1111/ecog.07364，即本包所引用的条目）是本项目承载的唯一定义性出处，但参考包自带的 vignette 仍写作 "Araujo et al. (in review)"，也没有给出 multisite 公式，因此我们无法对照某个发表过的方程去核验"`a` 按样区、`b+c` 按数对"。上述内部不一致已如实如此陈述，并在该核对完成之前保持为可选项。若确需已发表的多站点指标（Baselga 2010 / Carrasco et al. 2013 一族），请用 `approach="pairwise"` 并显式聚合。

## 7. 累积速率（CpD、CpE、CpB）

```python
import phyloslicer as ps

res = ps.cpd_rate(tree, mat, n_slices=100)                 # CpD + PD + pDO
res = ps.cpe_rate(tree, mat, n_slices=100)                 # CpE + PE + pEO
res = ps.cpb_rate(tree, mat, adj,
                  component="sorensen",                    # or turnover/nestedness
                  approach="multisite",                    # or pairwise
                  weighted=False,                          # True -> CpB_RW
                  n_slices=100)                            # CpB + PB + pBO
```

每个函数返回一个 DataFrame，每行一个集合区：

```
site_id, CpD, PD, pDO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
site_id, CpE, PE, pEO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
site_id, CpB (CpB_RW), PB (PB_RW), pBO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
```

末尾四列是拟合与切片的审计线索：`at_lower_bound` / `at_upper_bound` 说明所报告的速率是边界最优而非估计值；`ultrametric` 记录被切片的树是否通过超度量检验（`False` 只有在传入 `allow_non_ultrametric=True` 时才可能出现）；`origin_within_tree` 说明拟合出的起源时刻是否落在树的范围之内。`pXO = -ln(1 - percentile) / rate` 是对拟合指数曲线的外推，因此速率偏慢时会把起源放在**根之前**，而系统发育树中并不存在对应的节点；数值原样报告，由该标记连同 `reason` 加以说明。`posterior_rate` 与 `posterior_origin` 在整个抽样池上给出 `n_ultrametric_draws` 与 `origin_within_tree`。

常用参数：

* `criterion="time"`（默认）或 `"pd"`——切片如何切分。
* `percentile=0.95`——起源时间是拟合曲线处于 `1 - percentile` 的那个时刻，即仅有 5 % 的多样性累积完成之时（一个接近根的时刻）；它**不是**已经累积到 95 % 的时刻（`pXO = −ln(1 − percentile)/r`；`0.95` ≡ R 的 `pDO = 5`）。
* `diagnostics=True`——在 `RateFit` 对象中保留残差（超大规模运行可设为 `False` 以节省内存）。
* `multistart=8`——稳健逐行重拟合回退路径的起点数。
* 并行在向量化内核内部处理；单次稀疏矩阵乘 `A @ C` 已同时求解所有组合。
  已同时求解所有集合区。

不做拟合，只要剖面：

```python
prof = ps.rate_profile(tree, mat, n_slices=100, index="PD")   # "PD"|"PE"|"PB"|"PB_RW"
prof["values"]   # (n_sites, k) relative per-slice values (PB: beta shares)
prof["ages"]     # (k,)
prof["status"]   # per-site status ("ok", "empty_site", ... )
```

切片数敏感性分析（对标 `CpR_sensitivity`，另加一条启发式建议）：

```python
sens = ps.sensitivity(tree, mat, slice_grid=[10, 25, 50, 100, 200],
                      rate="cpd",              # "cpe"|"cpb"|"cpb_rw" (need adj)
                      sample_sites=0,          # 0 = use every site (R requires > 0)
                      seed=42)                 # deterministic site subsampling
# columns: site_id, n_slices, rate, converged, at_lower_bound,
#          at_upper_bound, suggested_slices
```

`suggested_slices` 是使某样区速率在相邻网格点之间的变化小于 5% 的最小网格计数。

### 7.1 直接调用拟合函数

```python
from phyloslicer.rates import fit_rate, origin_time

fit = fit_rate(profile, stack.ages)            # profile: (k,) root -> tip
fit.rate, fit.converged, fit.r_squared
fit.residuals, fit.start_used, fit.reason
origin_time(fit.rate, percentile=0.95)         # pXO in Ma
```

`fit_rate` 依次执行对数线性粗估、几何多起点网格和 bounded least-squares 精修；任一起点失败都不会中断整次调用。`method="mle"` 切换为高斯极大似然目标，配 Nelder–Mead。批量向量化内核 `fit_rates_batch(profiles, ages)` 一次拟合整个 `(n_sites, k)` 矩阵，是 `cpd_rate` / `cpe_rate` 内部所走的路径。

## 8. 跨后验树的不确定性

```python
from phyloslicer.uncertainty import posterior_rate, posterior_origin, bootstrap_sites

# sources: a list of TreeArray, a directory of .tre files, or a
# multi-tree Newick/NEXUS file — trees are streamed, not all loaded
post = posterior_rate("posterior_trees/", mat, rate="cpd", n_slices=100, hdi=0.95)
```

`posterior_rate` 返回逐样区汇总：

```
site_id, rate_mean, rate_median, hdi_low, hdi_high, ess_proxy,
ess_proxy_caveat, n_trees, n_converged
```

* 树逐棵处理。无法**读取**的抽样会被跳过，并警告指明文件与记录号（`<file>: skipping record j of N: ...`）；能读取但无法分析的抽样按 `skipping tree i: ...` 警告；只要有任何抽样被丢弃，还会追加一条结束性警告说明丢弃了几棵，因为此时汇总所依据的池已小于磁盘上的那个池。`iter_trees` 同时接受多树 Newick **字符串**、文件或目录。
* `hdi_low/high` 是各样区速率跨树的经验最高密度区间。
* `ess_proxy` 是由滞后自相关得到的有效样本量代理；远小于 `n_trees` 说明抽样是有序的（可能自相关）。
* 每棵树独立与矩阵对齐，因此枝端集合略有不同的后验树池也能被正确处理。

`posterior_origin` 返回起源时间 pXO 的 `site_id, origin_median, hdi_low, hdi_high`。

```python
boot = bootstrap_sites(tree, mat, n_boot=1000, rate="cpd", n_slices=100, seed=0)
# DataFrame: replicate, rate   (rate refit on the pooled assemblage of each
#                              bootstrap resample of sites)
```

## 9. 空间辅助工具

```python
from phyloslicer.spatial import (adjacency_from_grid, adjacency_from_bounds,
                                 adjacency_from_coords, to_adjacency_dict,
                                 from_adjacency_dict, rate_map)

adj = adjacency_from_grid(geodataframe, method="queen")  # or "rook"
adj = adjacency_from_bounds(bounds_df, method="rook")    # (n, 4) minx/miny/maxx/maxy
adj = adjacency_from_coords(x, y, method="knn", k=8)     # symmetrised k-nearest
adj = adjacency_from_coords(x, y, method="radius", r=2.0)
d   = to_adjacency_dict(adj)          # {"0": ["1", "5", ...], ...} for interop
adj = from_adjacency_dict(d)
```

约定：

* **对角线为 1**：焦点集合区总是包含自身，与 treesliceR 一致。
* `queen` 计入任何边界接触（边或角）；`rook` 要求共享正长度的边界段。GeoPandas 路径在真实多边形上实现同样的语义；没有 GeoPandas 时，请提供 `minx/miny/maxx/maxy` 列，改用矩形运算。
* `from_adjacency_dict`/`to_adjacency_dict` 与 R 的 `AU_adj` 风格稀疏表示互通。

映射：

```python
ax = rate_map(rates, grid_df, rate="CpD")            # bounds columns ...
ax = rate_map(post, grid_gdf, rate="rate_median")    # ... or GeoDataFrame
ax = rate_map(post, grid_gdf, rate="hdi_high", quantiles=True)
# NaN rates render in light grey; ax is a matplotlib Axes
```

## 10. 可视化

```python
from phyloslicer.viz import plot_rate_line, plot_sensitivity, plot_posterior

ax = plot_rate_line(profile, fit=fit, ages=stack.ages)   # cumulative + fit + pXO line
ax = plot_sensitivity(sens_df)                           # rate vs n_slices per site
ax = plot_posterior(post)                                # ECDF with HDI shading
```

> `ps.viz.plot_rate_line(..., direction="from_root")` 是默认值：x 轴从树根朝向现生方向递增，也就是本手册示例所使用的朝向。若要 treesliceR 视图（x 为距今多少 Ma 的年龄），请传 `direction="from_present"`。

所有函数都返回 matplotlib `Axes`，默认使用静态的出版级样式，且绝不调用 `plt.show()`。

## 11. 模拟与内置数据

```python
from phyloslicer.simulate import yule_tree, birth_death_tree, make_dataset

tree = yule_tree(100, age=1.0, seed=1)             # ultrametric, exactly n tips
tree = birth_death_tree(100, age=1.0, birth_rate=1.0, death_rate=0.2, seed=1)
# grows until 100 lineages are extant, then rescales the shape onto age=1.0
mat, coords = make_dataset(n_sites=100, n_species=50,
                           structure="clustered",  # disk-shaped ranges
                           seed=1)                 # or "random"
```

内置示例（模拟规模对标 treesliceR 澳洲雀形目案例研究——网格、308 个物种、后验树池）：

```python
from phyloslicer.datasets import sim_passerines

data = sim_passerines(n_sites=208, n_species=308, n_trees=5, seed=42)
# data["trees"] (list[TreeArray], labels sp_1..), data["mat"] (sites x species,
# columns == tip labels), data["coords"], data["bounds"], data["adj"]
```

`data["mat"]` 的列携带这些树的枝端标签，因此树与矩阵可以直接配套使用：

```python
res = ps.cpd_rate(data["trees"][0], data["mat"], n_slices=50)
```

## 12. 命令行界面

```
phyloslicer [--version]
phyloslicer slice     TREE (--rootward X | --tipward X | --n-slices N)
                      [--criterion {time,pd}] --out PATH
phyloslicer rates     TREE SITES [--rate {cpd,cpe}] [--n-slices N]
                      [--criterion {time,pd}] [--percentile P] --out CSV
phyloslicer posterior TREES SITES [--rate {cpd,cpe}] [--n-slices N]
                      [--hdi H] --out CSV
phyloslicer validate  TREE SITES --reference REF_CSV [--manifest PATH]
                      [--rate-column NAME] [--n-slices N] [--tolerance T]
                      [--no-strict] --out CSV
```

* `slice` 写出单个 Newick 文件；`--n-slices` 时写出一个目录，每个切片一个文件（`slice_0001.tre`，……）。
* `--rootward` / `--tipward` / `--n-slices` 有且只有一个被接受。
* 已存在的 `--out` 会被拒绝，除非给出 `--force`。这对 `--n-slices` 影响最大：没有这项检查，把 `slice` 重复输出到同一目录会残留上一代的 `slice_*.tre` 文件，于是 `posterior` 会把两次运行的切片混作一池，而 `slice` 只报告它本次写出的数量。
* `--log FILE` 写出 JSON 运行摘要（版本、命令、`elapsed_s`、`started_at`、关键统计量）。**没有** `--seed`：CLI 不做任何随机抽样，而 `sensitivity(..., seed=)` 经由 API 提供可复现性。
* `validate` 将累积速率运行与 treesliceR 参考 CSV 比对。被比较的列是**来源清单**登记为 `rate_column` 的那一列（因此 `CpD` 与 `CpE` 参考都可用，报告会按实际比较的指标命名列）；`--reference` 若在清单中没有条目则一律拒绝，退出码 2。清单可通过 `--manifest`、`$PHYLOSLICER_REFERENCE_MANIFEST`、参考文件同目录下的 `manifest.txt`，或 `validation/reference/manifest.txt` 提供。β 类参考（`CpB`、`CpB_RW`）需要邻接矩阵，而本子命令不接受该参数——请用 `python validation/compare.py` 为它们评分。
* `--n-slices` 必须等于清单头部记录的 `n_slices`，因为两个不同的窗口数不是同一个量。
* 一致性规则：双侧均有限的按相对误差 ≤ 容差比较（默认 `1e-4`，与随包参考实际达到的约 2e-6 优化器一致水平相称）；**双侧**均为 NA 的样区计为一致；单侧 NA 计为不一致。运行时会在通过率旁打印实测的最大相对误差。
* 退出码：`0` 成功；`1` 数值不一致（可用 `--no-strict` 抑制，仅供探索性使用）；`2` 请求被拒绝或格式有误。因此 `validate` 可以用于把关 CI：等价率不足 1 即为非零退出。

## 13. 从 treesliceR（R）迁移

| treesliceR | PhyloSlicer | 说明 |
|---|---|---|
| `squeeze_root(t, x, criterion)` | `ps.slice_rootward(t, x, criterion)` | `"my"` → `"time"`，`"PD"` → `"pd"` |
| `squeeze_tips` | `ps.slice_tipward` | 同上 |
| `squeeze_int` | `ps.slice_interval` | `invert=TRUE` → `invert=True`，但主 API 返回 `(older, younger)`，而 R 返回 `(younger, older)`——要 R 的顺序请用 `ps.compat.squeeze_int` |
| `phylo_pieces(n, method)` | `ps.slice_pieces(n=` 或 `width=)` | `method=1` → `n=`，`method=2` → `width=`；`SliceStack.ages` 是递减的 Ma-before-present，而 R 的 `timeSteps` 是递增的累积深度——要 R 的向量请用 `ps.compat.phylo_pieces` |
| 等 `"PD"` 切片 | 等 `"pd"` 切片 | **数值上不可比**：见 §13.1 |
| `prune_tips(method=1/2)` | `ps.prune_tips(side="after"/"before")` | `qtl=TRUE` → `quantiles=True` |
| `r_phylo(index=)` | `ps.rate_profile(index=, normalize=False)` | R 函数对 PD/PE 返回原始逐切片值，因此须用 `normalize=False` 才能复现它；PhyloSlicer 默认返回归一化到和为 1 的剖面 |
| `CpD / CpE` | `ps.cpd_rate / ps.cpe_rate` | 增加 `converged`、`r_squared`、`reason` |
| `CpB / CpB_RW` | `ps.cpb_rate(weighted=, multisite_domain=)` | **凡参考实现能给出的 11 个 β 指标族都已交叉验证**（§15）；它返回不了的 3 个没有参考值，而 turnover/nestedness 的速率属于已记录的归一化分歧，不是不一致。在默认 `multisite_domain="mixed"` 下，multisite 序列跟随参考算术，因而会随邻域大小漂移（§6.1）；`multisite_domain="paired"` 是一个**可选的按数对求和域**（它是参数取值，不是另一种 `approach`/"paired method"），消除漂移的代价是不再与 R 可比。`approach="pairwise"` 聚合采用成对均值，因为 R 的 `CpB` pairwise 分支量纲不一致——实测：该分支会在其自身的 `nls` 调用中中止 |
| `CpR_sensitivity(_plot)` | `ps.sensitivity` / `ps.viz.plot_sensitivity` | `sample_sites=0` = 全部样区；另有建议值 |
| `CpR_graph` | `ps.spatial.rate_map` / `ps.viz.plot_rate_line` | matplotlib 后端 |
| — | `ps.posterior_rate / posterior_origin` | 新增：后验不确定性 |
| — | `ps.bootstrap_sites` | 新增：样区重抽样 |
| — | 邻接构建器、CLI、GeoPandas 互操作 | 新增 |

即插即用的兼容层保留了 R 风格的签名：

```python
from phyloslicer.compat import squeeze_root, squeeze_tips, squeeze_int, phylo_pieces, prune_tips

pieces, time_steps = phylo_pieces(tree, n=100, criterion="my", method=1, timeSteps=True)
# time_steps is R's ascending cumulative depth from the root, cumsum(rep(n, nslices)):
# time_steps[j] is the depth reached BY pieces[j], and time_steps[-1] == tree depth.
# (The main API's SliceStack.ages is the opposite convention; ask for it with
#  phylo_pieces(..., timeSteps_as_depth=False).)

complement = squeeze_int(tree, from_=7.0, to=1.0, invert=True)
# complement[0] == squeeze_root(tree, 1.0)   -> the younger/crown piece (R's tree1)
# complement[1] == squeeze_tips(tree, 7.0)   -> the older piece        (R's tree2)
```

上述两处顺序均由 `tests/test_compat.py` 断言钉住。

### 13.1 不只是改名的行为差异

* 静默的 `NA` 行被替换为 `converged=False` 加上明文 `reason`。
* 每次树/矩阵调整（删物种、剪枝端）都会警告。
* CLI、随机种子与 JSON 运行摘要使每次运行可审计。
* `n_slices` 对应 R 的 `n`（切片数量），不是宽度。
* **等 PD 切片（`criterion="PD"`）是有意的、未经交叉验证的偏离。** R 把分支数硬编码为 `2, 3, 4, …`
  （`squeeze_root.R:98`、`phylo_pieces.R:159`：
  `nBranch = 2:(length(unique(nodes$YearBegin)) + 1)`），这预设了树是严格二叉且各节点深度互异。在多叉树上，这套算术连总 PD 都不守恒：对
  `(((A:2,B:2):1,C:3,D:3,E:3):1,(F:4,G:4):0);`（T = 4，总 PD = 23），R 的断点只累积到 13 的 PD，而该树有 23。PhyloSlicer
  统计每个深度区间内*真正活跃*的分支数，因此其 `"pd"` 窗口按构造累积到全树 PD。这是一处修正而非移植，而且**未经交叉验证**——`validation/reference/`
  中的每个参考文件都以 `criterion = "my"` 生成（见 `validation/NOTES.md`）。
* `squeeze_root`/`squeeze_tips`/`squeeze_int`/`phylo_pieces` 会拒绝
  `{"my", "million years", "pd", "time"}`（大小写不敏感）之外的 `criterion`，`phylo_pieces`/`prune_tips`
  会拒绝 `{1, 2}` 之外的 `method`，而不是像 R 那样静默"不切"或返回 `NULL`。
* multisite β 的求和域可自选更改（§6.1）；默认仍保持 R 的算术。

## 14. 诊断信息、状态与报错说明

### 14.1 `pb_per_slice` 状态字符串

| status | 含义 |
|---|---|
| `ok` | 邻域正常分析 |
| `fewer_than_two_sites` | `adj` 的焦点行（含自身）不足 2 个样区 |
| `degenerate_neighbourhood` | 所有样区为空，或所有样区共享完全相同的组成（β 无定义） |
| `zero_denominator` | 某个 β 分母消失（数值上退化） |

被跳过的焦点样区会以 `converged=False` 行返回，状态写在 `reason` 中——绝不会是静默的 `c(NA, NA, NA)`。

### 14.2 `RateFit.reason` 的取值

`"fewer than two slices"`、`"non-finite input values"`、
`"profile sums to zero (empty or degenerate assemblage)"`、
`"no start converged"`、`"optimiser did not report convergence"`、
`"vectorised Newton did not satisfy criteria"`。

### 14.3 常见校验报错

* `slice_rootward(tree, t)` 且 `t > tree.root_age` → *"the threshold inputted ... is bigger than the available by the phylogeny"*。
* `slice_interval(tree, a, b)` 且 `a <= b`（time 准则）→ *"the thresholds set in arguments [start] and [stop] are incompatible"*。
* `slice_pieces(tree)` 未给 `n`/`width`，或两者都给了 → *"pass either n or width"*。
* `sensitivity(..., rate="cpb")` 未给 `adj` → *"rate=cpb requires an adjacency matrix"*。
* `cpb_rate(...)` 的 `adj` 行数与矩阵不一致 → *"adj must have one row per site of mat"*。

## 15. 数值验证

比对集是 `{PD, PE, CpD, CpE}`，外加参考实现在给定输入上能够返回的 β 多样性族，全程把 treesliceR 1.1.0 严格当作黑箱运行（R 脚本在 `validation/r/`，比对在 `validation/compare.py`，来源登记在 `validation/reference/manifest.txt`，容差依据与缺口清单在 `validation/NOTES.md`）。在这个集合内，今天真正与 R 对上的只有 `CpD` 与 `CpE` 两个文件：

* CpD 与 CpE 速率：跨全部集合区的最大相对误差 ≈ 2 × 10⁻⁶（容差 10⁻⁴；残余差距是优化器容差）；同一文件里的 PD 与 PE *总量* 相符到 ≈ 4 × 10⁻¹³。这是本项目目前能做出的唯二真正外部比对。
* 逐切片 PD 与 PE 剖面（`*__rPD.csv` / `*__rPE.csv`）：在 18000 次比较中相符到 ≈ 3.4 × 10⁻¹¹，而且这是真正的 Python 对 R 的比较。它们曾经是例外：归档副本以 slice-by-site、30 个切片列的布局存储，而这种布局是已提交的驱动（50 切片、site-by-slice）不可能产出的，于是清单把它们声明为 `status = unverified`，`validation/compare.py` 跳过它们而不是静默转置。在 R 4.5.3（aarch64-apple-darwin20.0.0）配合 ape 5.8.1 与 treesliceR 1.1.0 的环境下运行 `validation/r/generate_reference.R` 了结了这件事——其余每个参考文件都逐字节一致，而该次运行写出的剖面现已登记。被替换的旧副本不随仓库分发，证据就是上面这次重跑本身。

本比对**不**覆盖以下内容（全部细节见 `validation/NOTES.md`；不要把上面四个数字读作广义的等价证明）：

* **3 个 β 指标族**（`CpB` Sorensen/pairwise，以及 `r_phylo(index = "PB")` 的 Sorensen- 与 turnover-pairwise）——treesliceR 1.1.0 在这些输入上对它们给不出可比结果，因此没有参考值可评分。其余 11 个 β 指标族**已**比对：24 条逐切片累积 β 剖面行到 1.1e-9，1260 个整树 `PB`/`PB_RW` 总量到 7.2e-13，540 个通过校验的速率拟合在优化器容差下到 2.9e-6。
* **turnover 与 nestedness 的 β *速率***——已比对，并且相对误差最大达 196%：本包把剖面重新归一化为和 1，而参考实现拟合的是 `cumsum(profile) / PB_whole_tree`。每一这样的行都按两条曲线同时评分，而按参考约定重新拟合即可复现 R，所以这一差距是被钉住的归一化选择（见 validation/NOTES.md），不是算术；同一文件里的整树 `PB` 仍相符到 7.2e-13。
* **等 `"PD"` 切片**——每个场景都以 `criterion = "my"` 运行。`"PD"` 判据有意偏离 R（§13.1），且从未被比对过。
* **枝端修剪路径**——生成器保证每个物种至少出现一次，因此 `align_tree_matrix` 从不剪枝，该分支未被演练。
* **退化集合区的 NA 模式**——`compare.py` 现在逐元素比较 NA *模式*而不再只统计非有限值的个数，但参考生成器让每个物种至少出现一次，两侧都没有 NA 行，因此在归档集合内这项检查无事可查。此前"NA 模式在每个退化集合区上都一致"的说法缺乏支撑，应视为已移除。
* 参考树只是随包生成器产出的那几棵（20/100 尖的出生-死亡树与共合并树、单位树龄、`n_slices = 50`、30 个样区）；这里没有任何东西能验证形状或树龄不同的树。

测试套件另外还钉住：切片的解析解、指数的暴力法交叉核对、逐切片划分不变量、兼容层的返回顺序，以及 I/O 往返。它的大小有意不在这里给出：手册里的数字在任何人为仓库加一条测试的那一刻就过期了。请以树本身为准——

```bash
grep -c "^[[:space:]]*def test_" tests/test_*.py                               # functions
python -m pytest --collect-only -q -o addopts="" 2>/dev/null | grep -c "::"     # cases
```

* 第一条统计 `def test_` 函数，第二条统计 pytest 实际*收集*到的东西，两者不同，因为参数化的函数会按每个参数组展开成一个用例。`-o addopts=""` 是必需的：`pyproject.toml` 里设了 `addopts = "-q"`，因此裸跑 `pytest --collect-only -q` 相当于 `-qq`，只按文件打印计数而不输出节点 ID（此时 `grep` 返回 0）。
运行全部内容：

```bash
python -m pytest tests/
python validation/compare.py validation/reference
```

## 16. 常见问题与故障排查

**Q：矩阵列名与树枝端标签对不上。**
A：`align_tree_matrix`（由速率函数自动调用）会修好重叠部分，并对每一次调整发出警告。要显式控制对齐——或消除本可避免的警告——请先自行子集化矩阵。

**Q：我可以用非超度量树吗？**
A：默认不行：`slice_pieces()` 会抛出 `ValueError`，因为切片语义要求所有枝端终于
同一年龄。传入 `allow_non_ultrametric=True` 可得到"告警并继续"的行为，它会在结果
表中记录 `ultrametric = False`，但切片年龄依旧失去"时间 *t* 之前的多样性"这一
含义。请先重标度或做超度量化。

**Q：有些样区返回 `converged=False`，原因是 "profile sums to zero"。**
A：该集合区为空。这属于被报告的情况，不是致命错误。按你的问题决定是填充还是剔除该样区。

**Q：某个样区的 CpD 大得离谱 / pDO 不合理。**
A：查看 `r_squared` 并用 `plot_rate_line(profile, fit, ages)` 检查。近似平坦或由单个切片主导的剖面会把速率推到边界（50）；这类边界解会通过诊断字段暴露出来，而不是被藏起来。

**Q：我该用多少个切片？**
A：在一个网格上跑 `sensitivity`；逐样区采用 `suggested_slices`，并在方法部分报告所用的网格。

**Q：树的支持度值去哪了？**
A：`TreeArray` 只存枝端、边和枝长；内部标签在读入时被丢弃。需要支持度请保留原始 Newick/NEXUS 文件。

**Q：`percentile=0.95` 与 R 的 `pDO = 5` 是一回事吗？**
A：是——`pXO = −ln(1 − 0.95)/r = −ln(0.05)/r`，即 R 的 `−log(pDO/100)/r` 取 `pDO = 5`。

**Q：有内置并行吗？**
A：并行在向量化内核内部完成——单次稀疏矩阵乘 `A @ C` 已同时计算所有样区；如需跨
树或跨数据块并行，请使用工作流引擎。

## 17. API 速查表

```python
# core
ps.TreeArray.from_newick(x) / .from_nexus(p) / .from_dendropy(o) / .from_biopython(o)
tree.root_age / .n_tips / .tip_labels / .total_pd() / .is_ultrametric
tree.subarray(labels) / .to_newick() / .to_dendropy() / .to_biopython()

# slicing
ps.slice_rootward(tree, time, criterion="time", collapse=False)
ps.slice_tipward(tree, time, criterion="time", collapse=False)
ps.slice_interval(tree, start, stop, invert=False, criterion="time", collapse=False)
ps.slice_pieces(tree, n=None, width=None, criterion="time")  -> SliceStack
ps.prune_tips(tree, threshold, quantiles=False, side="after")
SliceStack: .windows .ages .widths .contribution(subset) .contribution_matrix()
            .to_sliced_trees(collapse=False)

# indices
ps.indices.pd(tree, mat) / .pe(tree, mat) / .pb(tree, mat, adj, component, approach, weighted)
ps.indices.pd_per_slice / .pe_per_slice / .pb_per_slice
ps.indices.site_edge_membership / .edge_range_sizes

# rates
ps.fit_rate(profile, ages, method="lsq", multistart=8, bounds=(1e-6, 50)) -> RateFit
ps.origin_time(rate, percentile=0.95)
ps.cpd_rate(tree, mat, n_slices=100, criterion="time", percentile=0.95, ...)
ps.cpe_rate(...);  ps.cpb_rate(tree, mat, adj, component=..., approach=..., weighted=..., multisite_domain=...)
ps.rate_profile(tree, mat, adj=None, n_slices=100, index="PD")
ps.sensitivity(tree, mat, adj=None, slice_grid=(...), rate="cpd", sample_sites=0, seed=None)

# uncertainty
ps.posterior_rate(trees, mat, rate="cpd", n_slices=100, hdi=0.95)
ps.posterior_origin(trees, mat, rate="cpd", percentile=0.95, hdi=0.95)
ps.uncertainty.bootstrap_sites(tree, mat, n_boot=1000, rate="cpd", seed=None)
ps.uncertainty.iter_trees(source)

# spatial
ps.spatial.adjacency_from_grid(gdf, method="queen")
ps.spatial.adjacency_from_bounds(bounds, method="rook")
ps.spatial.adjacency_from_coords(x, y, method="knn", k=8, r=None)
ps.spatial.to_adjacency_dict(adj) / .from_adjacency_dict(d)
ps.spatial.rate_map(rates_df, grid, rate="CpD", quantiles=False)

# viz / simulate / datasets / compat
ps.viz.plot_rate_line / .plot_sensitivity / .plot_posterior
ps.simulate.yule_tree / .birth_death_tree / .make_dataset
ps.datasets.sim_passerines()
ps.compat.squeeze_root / .squeeze_tips / .squeeze_int / .phylo_pieces / .prune_tips

# io
ps.io.read_newick / .parse_newick / .read_nexus / .write_newick / .write_nexus
ps.io.presence_matrix / .align_tree_matrix / .as_matrix
```

## 引用

如果你使用 PhyloSlicer，请同时引用本软件包与方法学来源：

> Araujo, M.L., Ferreira, L.G.S.S., Nakamura, G., Coelho, M.T.P., Rangel, T.F. (2025)
> 'treesliceR': a package for slicing phylogenies and inferring phylogenetic patterns
> over evolutionary time. *Ecography* 2025, e07364. https://doi.org/10.1111/ecog.07364
>
> **关于年份。** 该文于 2024 年 10 月 28 日网络首发，作为文章 e07364 归入 *Ecography* 2025 卷；此处按卷号年份引用。已于 2026 年 9 月 26 日依据 doi:10.1111/ecog.07364 的 Crossref 记录核对（journal-article，Wiley，volume 2025，article-number e07364，issued 2024-10-28，print 2025-01）；该记录按顺序列出的五位作者，与生成验证参考所用的 treesliceR 1.1.0 构建的 DESCRIPTION 中的五位一致。该包自带的 vignette 仍保留旧的 "in review" 占位写法，期刊记录已取代它。若你的引用规范以网络首发日期计年，请改引 2024——DOI 不变。

## 许可证

代码：BSD-3-Clause。文档：CC-BY 4.0。
