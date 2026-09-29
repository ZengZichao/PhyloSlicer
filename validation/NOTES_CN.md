# 交叉验证说明（验证协议的第 4 步）

本文件是 `NOTES.md` 的中文译文，正文以 [NOTES.md](NOTES.md) 为准（无后缀的英文文件是记录版本）。

## 容差依据（tolerance rationale）

| 量 | 容差（tolerance） | 理由 |
|---|---|---|
| 逐切片 PD / PE 剖面（`r_phylo` 原始切片） | rtol 1e-8 | 完全确定性的算术；残余差异只来自浮点求和顺序（稀疏乘积 vs 逐边求和）以及树的 12 位 Newick 往返读写 |
| CpD / CpE 速率 | rtol 1e-4 | *速率*是一个非线性优化的解：treesliceR 用 `stats::nls`（单起点 r = 0.2，Gauss-Newton）拟合，而 phyloslicer 用有界多起点最小二乘。两者收敛到同一目标的同一个最优值，但只到优化器容差的精度（约 1e-5 相对值），因此 1e-4 才是合适的等价带 |
| pDO / pEO，pBO | 跟随速率 | 由 `-ln(1 - percentile) / rate` 导出 |
| 逐切片 PB 剖面（`r_phylo(index = "PB")`） | rtol 1e-8 | 与 PD/PE 相同的确定性算术，驱动它的 queen 邻接是由生成器作为输入导出，而不是从矩阵重建 |
| CpB / CpB_RW 速率 | rtol 1e-4 | 与 CpD/CpE 相同的优化器差异。在两个实现对拟合曲线归一化方式不同的地方（见下文*β 多样性分歧*），该行会被打分两次：一次基于随包发布的剖面，一次基于参考自己的曲线，且只允许第二次去弥合差距 |
| PB / PB_RW 整树总量 | rtol 1e-8 | 参考把它与自己拟合的结果并列报告，所以它是一个确定性量，不是优化器的输出 |

## 退化情形

treesliceR 对空群落组合（assemblage）返回 `c(NA, NA, NA)`，并让 `nls`
的错误中止整个调用。phyloslicer 返回 `converged=False`，附带一个明确的
`reason` 和 `NaN` 值。比较器把两侧的 NA 都视为一致，并把 NA
模式的差异报告在 `note` 列中。

## β 多样性分歧（已记录的设计决策）

treesliceR 的 `CpB` pairwise 分支构造一个 k x P 的逐对相对值矩阵，然后把它
的 `cumsum()` 与一个长度为 k 的 `age` 向量一起喂给 `nls` —— 维度上不一致。
`CpB_RW` 的 pairwise 分支则在每个切片内部对逐对比值取平均。PhyloSlicer
对每一个组分／方法（component/approach）组合都采用逐对平均值（单一代码路径，
记录在 `indices.api._pb_series` 的 docstring 中）。

对源码的这一解读如今已经是一次测量：在 2026-09-21 于 R 4.5.3 / ape 5.8.1 /
treesliceR 1.1.0 下的那次运行中，Sorensen 的 pairwise 分支在它自己的优化器
内部失败，报 `task 1 failed - "'qr' and 'y' must have the same number of rows"`
（`CpB.R:479` 处的逐对向量对 `CpB.R:525` 处的切片向量做了循环补齐），而
Sorensen 与 turnover 的 pairwise `r_phylo(index = "PB")` 输出返回的是一个展平的
k x P 块，其元素长度为 300 / 750 / 1800，而不是每个群落组合一条 50 切片的曲线。
生成器中的 `write_beta_rates()` / `write_beta_profile()` 拒绝这些形状，而不是把
它们强行塞成一份参考，因此这三个族在 `validation/reference/` 中是*缺失*的，
报告把它们记为 `no reference` 行——一个被点名的缺口，绝不是一次沉默的通过。

第二处分歧是一个约定，而它正是"某一行 β 速率读到 `documented divergence`、
周围的一切却都已核验通过"的原因：参考拟合的是
`cumsum(profile) / PB_whole_tree`，而 PhyloSlicer 在拟合之前把剖面重新归一化到
和为 1（一个被钉住的设计决策）。在周转（turnover）与嵌套（nestedness）、两种
方法（approach）上，这会让拟合速率最多偏离 1.96 倍。因此比较器会在参考自己的
曲线上为这些行重新拟合（`fit_reference_convention`）；重拟合复现了参考速率，
而同一份参考文件里的整树 `PB` / `PB_RW` 列与我们的结果相符至 7.2e-13。
那一行记录的是归一化的选择，不是算术上的差异。

> **仍然暂定的是什么。** 多点位 β *指数*随邻域大小增长，
> 因为参考约定把共享项分量跨样区求和，
> 而周转分量是跨配对求和。PhyloSlicer 有意跟随该约定；
> 来自不同邻域大小的数值
> 不得彼此比较。

## 结论

在 `Rscript validation/r/generate_reference.R` 之后运行
`python validation/compare.py`，把打印出的摘要连同所用的场景清单一并追加到这里。

### 当前状态

- **PD/PE 逐切片剖面**（rPD/rPE）：**已验证。** 2026-09-21 在
  `micromamba env r-4.5.3` 下运行：R 4.5.3 (aarch64-apple-darwin20.0.0), ape 5.8.1, treesliceR 1.1.0。那次运行
  写出的剖面已登记在
  `validation/reference/manifest.txt`（`site-by-slice`，30 行 x 50 个切片列，sha256 已钉住），并由 `compare.py` 打分：相对于 1e-8 的确定性容差，
  **18000 次比较中最大相对误差 3.396e-11**，通过率 100%，其中 **157**
  个单元格逐位相同——即报告的 `identical_cells` 列，其计数与最大值取自同一批有限
  单元格。在那次运行之前，没有
  办法确立存档文件的来源（第 4 项），这两个族被声明为 `status = unverified`
  并被跳过；这段历史保留在清单的 `limitations` 块里，而不是被删改抹去。
- **CpD / CpE 速率**：真正做过比较，最大相对误差 ~2e-6（容差
  1e-4）；同一批文件里的 PD 与 PE *总量*相符至 ~4e-13。这是仅有的两个
  *α* 速率族；β 速率是 `CpB` 与 `CpB_RW`，完整族列表在下面一条里。
- **CpB、CpB_RW 和 `r_phylo(index = "PB")`**：**在参考能给出对应物的地方
  完成了交叉验证。** 14 个 β 族里有 11 个有已登记的参考并被打分：24 行
  逐切片累积 β 剖面相符至 **1.08e-09（36000 个数值）**（容差 1e-8；其中 10244 个单元格
  逐位相同，多数是数值上为零的最老切片），42
  行整树 `PB`/`PB_RW` 相符至 **7.2e-13（1260 个数值）**，18 行已验证的速率相符至
  **2.88e-06（540 个群落组合拟合）**（容差 1e-4）。另有 24
  行速率——周转与嵌套，两种方法——由于上文的归一化原因被记为
  `documented divergence`，18 行（3 个族 x 6 个场景）因 treesliceR 1.1.0
  没有返回可比对象而被记为 `no reference`。这三类都不计入通过，
  `compare.py --strict` 对它们全部失败。
- **`r_phylo` 原始 PD/PE/PB 剖面**：这些*就是*等价性报告中的 `:PD` / `:PE` 行，
  而在 2026-09-21 于 R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 下的那次运行之后，它们是
  真正的 Python 对 R 的比较（最大相对误差 3.396e-11），不是
  "Python 对一份来源不明的 CSV 文件"——第 4 项记录了发生了什么变化。
  `r_phylo(index = "PB")` 同样被比较，只要参考为每个群落组合返回一条曲线
  （见上面 β 那一条与第 2 项）。
- **`criterion="PD"`（等 PD 切片）**：**仍未交叉验证**，而现在的理由是一个
  被实证的原因，不是缺少工具。
  生成器中受保护的 `criterion = "PD"` 代码块*确实*运行过（2026-09-21 于 R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 下的那次运行），并在全部六个场景上失败（18 条警告，未写出任何 `*__PD.csv`）。失败信息是 `Error: object 'cutted_tree' not found`，而不是
  这些说明的早前一个版本记录的 `length of 'dimnames' [2] not equal to array extent`：`phylo_pieces.R:42` 和 `:152` 在小写的 `criterion ==
  "pd"` 上分支，而 `CpD`/`CpE`/`r_phylo` 文档写明并向下传递 `"PD"`，于是
  任何切片都不会被计算。因此被文档写明的拼写根本到不了任何比较。（是否*应该*
  使用那个未被文档记载的小写拼写，是另一个决定，第 3 项的结构性分歧对它同样适用。） 从这些参考出发，关于等 PD 切片不能在任一方向上做出任何数值论断；第 3 项
  描述的结构性分歧仍只是对源码的解读，不是测量。
- **剪枝路径**（矩阵缺失整个演化支（clade））：不在当前验证集的
  覆盖范围内；参考生成器保证每个物种至少出现一次，因此 `align_tree_matrix`
  从不触发剪枝。
- **比较集是 `{PD, PE, CpD, CpE}` 加上 β 集**（`CpB`、
  `CpB_RW`、`r_phylo(index = "PB")`），**只覆盖 treesliceR 能产出的那些族。**
  超出该范围的一切，等 PD 切片包括在内，都没有外部参考，
  无论本文件的标题数字看起来如何。


## 本次交叉验证未覆盖的内容（在引用任何通过率之前请先读）

1. **三个 β 族没有外部参考**：`CpB`
   `sorensen`/`pairwise`，以及 `sorensen`/`pairwise` 和
   `turnover`/`pairwise` 的 `r_phylo(index = "PB")`。它们不是测试框架里的缺口——treesliceR 1.1.0 不为它们
   返回可比对象（已测量，见*β 多样性分歧*），而生成器也拒绝写出一份它不得不
   凭空造出的形状。关于这三个族的一切论断仍只是内部一致性的结果。
2. **累积-β 参考现已存在，覆盖参考能产出的那些族。** `generate_reference.R` 对每一个组分×方法组合调用 `CpB`、`CpB_RW` 和
   `r_phylo(index = "PB")`，把它所使用的 queen 邻接作为输入导出（`*__adj.csv`，因此 Python 侧是在参考实际
   消费的那些邻域集合上被打分的），并登记返回的内容：报告 126 行 β 行中有 84 行已验证。第 1 项
   点名的三个族完全不产出文件，而 A2 归一化约定的四个速率族被记为分歧而不是通过。`r_phylo` 的 PD/PE
   输出*确实*被写出，并回到报告的 `:PD` / `:PE` 行上，是一次真正的 Python 对 R 的比较——该结论如何在 2026-09-21 落定见第 4 项。
3. **`criterion="PD"`（等 PD 切片）从未被比较。** 全部六个存档
   场景都运行 `criterion="my"`。这一偏离是有意为之并有记录：
   PhyloSlicer 统计真实的活动分支，因此即使多叉树也守恒 PD 总量，
   而参考对每一层深度硬编码分支数
   `2, 3, 4, ...`（`squeeze_root.R:98`、`phylo_pieces.R:159`、
   `nBranch = 2:(length(unique(nodes$YearBegin)) + 1)`），因此在这样的树上它甚至
   不守恒 PD 总量 —— 例如在
   `(((A:2,B:2):1,C:3,D:3,E:3):1,(F:4,G:4):0);`（T = 4，PD 总量 = 23）上，参考
   的切点累加到 13，而 PhyloSlicer 累加到 23。
   `validation/r/generate_reference.R` *尝试*在 `criterion = "PD"` 下生成
   `CpD`/`CpE`/`r_phylo` 参考（后缀为 `__PD.csv` 的文件）。**它们无法被产出**：该代码块被运行过（2026-09-21 于 R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 下的那次运行），每一次调用都在
   treesliceR 内部以 `object 'cutted_tree' not found` 失败（*结论*一节描述的
   小写 `criterion == "pd"` 分支问题），因此 `validation/reference/` 中不存在任何
   `*__PD.csv`，也不得把任何一个加入通过率表。即使未来的 treesliceR
   版本能产出它们，脚本中的警告代码块仍然适用：存档场景的树都是二分的，
   因此那里的一致并不能验证两个实现产生分歧的多叉树情形。
4. **`*__rPD.csv` 与 `*__rPE.csv` 的来源：已于 2026-09-21 通过运行生成器解决。** 2026-09-21 在 `micromamba env r-4.5.3` 下运行：R 4.5.3 (aarch64-apple-darwin20.0.0), ape 5.8.1, treesliceR 1.1.0。`Rscript
   validation/r/generate_reference.R <dir>` 逐字节复现了它与 `validation/reference/` 共有的 **21 个**时间准则文件中的 **21 个**——三棵 Newick 树、
   六个存在矩阵以及全部十二份 `*__CpD.csv` / `*__CpE.csv`
   参考——这就钉住了存档集背后的随机种子、RNG 流与包版本。唯一不同的输出是那十二个
   剖面文件：生成器写出 **30 行 × 51 列**（`site_id` +
   `slice_1..slice_50`），也就是**site-by-slice、50 个切片列**，与其源码所述完全一致，
   而提交的文件是 **50 行 × 31 列**、**slice-by-site、30 个列**，
   `site_id` 先跑 1..30 再跑 1..20（把一个 `seq_len(nrow(mat)) = 30` 的
   向量粘到 50 行的块上，这是 R 的循环补齐签名）。所以提交的那批剖面出自本脚本更早的一个版本，从来不是已登记输入的配套物；重新生成的那批才是，而 `manifest.txt` 现在把它们登记为
   `status = registered`，这正是把报告的 `:PD` / `:PE` 行从"Python 与这些 CSV 文件一致"变成"Python 与
   treesliceR 一致"的原因。清单仍然做不到的：它的哈希能发现事后的改动，但哈希本身并不能重新推导出某个 CSV 就是
   被提交的脚本所写出的东西——只有上面那次运行能做到，而它记录在每份剖面的 `note` 里。被取代的文件不属于仓库：它们只存在于作者的工作历史中，所以
   这一矛盾是由上面的重跑落定的，而不是由出示旧的字节落定的。
   保留在案的原始测量：每份提交的文件是 **50 行 × 31 列**
   （`site_id` + `slice_1..slice_30`），即按 **slice-by-site 存储、30 个切片列**，因为存档场景集使用 30 个样区。而提交在案的
   `generate_reference.R` 却设置 `n_slices <- 50L` 并写出
   `matrix(unlist(rpd_ref), nrow = length(rpd_ref), byrow = TRUE)`，对一个
   30 样区的矩阵这会产出 **30 行 × 51 列** — **site-by-slice、50 个切片列**。两种朝向与两种列数都不同，所以
   要么是另一个驱动脚本产出了存档文件，要么它们是从 PhyloSlicer 自己的输出转置而来的——而 2026-09-21 的那次运行支持前者。`validation/compare.py` 的早前一个版本会自动转置形状不符的参考（"legacy column-bound layout"），正是这一点让该
   不匹配在那么长的时间里未被察觉。`manifest.txt` 为每个随包发布的文件钉住 sha256、行数、列数、朝向以及速率/总量/来源列。
5. **剪枝路径有单元测试，但从未在 R 比较中被真正走到**，因为
   每个 fixture 矩阵都包含树的全部枝端。
6. **退化群落组合上的 NA 模式一致性**在 `validation/compare.py` 中现在是对
   模式的真实比较（而不是非有限值的个数）（`na_pattern_report`，逐元素），但六个参考场景的构造方式保证没有
   任何一个群落组合为空 —— 每一份存档的 `*__CpD.csv` / `*__CpE.csv` 都完全有限（实测：`CpD`/`PD`/`pDO` 和 `CpE`/`PE`/`pEO` 全程 0 个 NaN）—— 所以
   这项检查在存档参考集内没有可咬住的东西。
7. **R 与包版本现在已被记录，但未被校验。**
   `validation/reference/manifest.txt` 带有 `written_in`、`r_available` 和
   `generator_sha256`，每份重新生成的剖面的 `note` 点名
   `R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0`；内置副本的 `DESCRIPTION`
   写着 `Version: 1.1.0`。运行时没有任何东西*断言*这些版本，所以
   "针对 treesliceR 1.1.0 验证过"是参考集的一个文档属性，而不是比较器被强制满足的
   前置条件——未来在一个不同的 treesliceR 下的一次运行可能与钉住的字节
   悄悄不一致。sha256 钉值能抓住这一点（不同的输出会让 `check_registration` 失败），这是较弱但自动的保证。

因此，本目录里 100% 的通过率意味着：*phyloslicer 在四个 α 族 `{PD, PE, CpD, CpE}` 和 treesliceR 能产出的 11 个 β 族之内，按所述容差复现了存档的 `criterion = "my"` 参考 CSV 中的数值* ——
逐切片 PD/PE 剖面也包括在内，因为生成器被运行过、清单也登记了它们（第 4 项）。它不延伸到等 PD 切片、剪枝路径、
退化群落组合的 NA 模式，也不延伸到第 1 项的三个 β 族，而 A2 归一化约定的四个
速率族报告为分歧而不是通过。

## 报告自身的可复现性

`validation/reference/equivalence_report.csv` 是参考目录加上每个场景两份
Newick 树／存在矩阵输入的字节函数，别无其他：

- 重新运行 `Rscript validation/r/generate_reference.R <dir>`（2026-09-21 于 R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 下的那次运行）
  逐字节重写 `validation/reference/` 的全部 **105** 个数据文件——三棵树、六个矩阵、六个邻接、十二份 CpD/CpE 速率、十二份 rPD/rPE
  剖面，以及整个 β 块（CpB、CpB_RW、rPB）。正是这一点使 β 比较成为与 treesliceR 的比较，而不是与这些 CSV 的比较。
- `compare.py` 钉住它的 BLAS 线程数（`OMP_NUM_THREADS` 及同类变量，在 numpy 加载之前设置），因为不加钉住，报告并*不*是字节可复现的：
  连续多次运行实测，162 行中有 7 行在最后一位有效数字上变动（一个剖面误差由 9.161e-14 -> 9.141e-14，一个拟合速率由
  1.47158405e-06 -> 1.47158404e-06），因为多线程归约按与调度相关的顺序求和。状态没有变化，但文档引用了这些
  数字，所以重点不在容差问题——在于可复现性。
  `tests/test_validation_chain.py::TestReportDeterminism` 在两个新进程里各跑一遍比较器并比对字节。
- 字节一致性成立在*同一个*解释器内部，不跨解释器。以 numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3 环境为基准实测：全部 162 个状态结果一致，7 行在末位数字上变动（9.161e-14 -> 9.141e-14），所以
  随包报告所支持的论断是"每一行状态相同、误差在有意义的那些位上相同"，而上面的 `written_in` 点名了产出这些存档数字所用的环境。
