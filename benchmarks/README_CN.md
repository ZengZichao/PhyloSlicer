# 本目录中的基准测试表

中文译文，正文以 [README.md](README.md) 为准。

这里放着四套彼此独立的协议。它们**不可**互换，有两个列名在不同文件里指的是不同的
东西，而现在每一行都会写明它是在哪个输入上测得的。文档中给出的每一条加速比都能追溯
到下面某一张表，并且 `recompute_ratios.py` 会把这条追溯关系打印出来。

## 协议

| 文件 | 计时对象 | 统计量 | 重复次数 | 跨语言？ |
|---|---|---|---|---|
| `results.csv` | `bench_core.py`：向量化内核 vs 进程内朴素基线 | **最小值** | 3（朴素基线：1） | **否** - Python vs Python |
| `results_py.csv`, `results_r.csv` | 验证场景输入上的切片内核 + CpD/CpE 拟合 | 中位数 | 5 | 是 |
| `results_scaling_py.csv`, `results_scaling_r.csv` | **仅切片内核**（`slice_pieces(...).contribution_matrix()` vs `phylo_pieces(...)`），每个树规模在 k = 25/50/100 的每一档上 | 中位数 | 15（+1 次被丢弃的热身运行） | 是 |
| `results_py_large.csv`, `results_r_large.csv` | 250 个组合、1000/2000 tips、k = 100 上的端到端 CpD | 中位数 | 5 | 是 |
| `results_memory_py.csv` | `benchmark_memory.py`：`E x k` 数组与 `k` 份数组级副本的 tracemalloc 峰值 | 单次运行，无重复 | 1 | **否** - 仅 Python 侧的模拟 |

`slice_s` 在 `results_py.csv`、`results_scaling_py.csv` 和 `results_py_large.csv`
里指的是 **E x k 贡献矩阵**，但在 `testing_suite/benchmarks/benchmark_results.csv`
里只指**窗口构建本身**。不要把这两个目录里的数字互相比较。

每张计时表都带一列 `reps`，写明该行是几次计时运行的中位数；`recompute_ratios.py`
遇到 `reps` 缺失、或与上面这套协议文档声明的次数不一致的表会直接拒绝通过。

## 每一行是什么——以及它不是什么

`results_py.csv` 与 `results_r.csv` 用 **12 行覆盖 4 个 (tips, n_sites, n_slices)
单元格**：网格是 3 棵树（`bd100`、`bd20`、`coal100`）× 2 种矩阵结构
（`clustered`、`random`）× k = 25/50，而 `bd100` 和 `coal100` 都是 100 tips。这些多
出来的行**并不是那 5 次重复**。每一行本身已经是它自己那 5 次计时运行的中位数；行与
行之所以不同，是因为它们测自**不同的树和不同的矩阵**。因此并不存在一个诚实的
“(100, 30, 25) 四行取中位数”——那会把两棵出生-死亡树和两次溯祖会话平均成一个没有
人测过的统计量——配对是按
`(tips, n_sites, n_slices, source, structure)` 一一对应的。本目录更早的版本用
`(tips, n_sites, n_slices)` 做索引并让**最后一行**胜出，这就是归档的网格比值读到的
是 `coal100__random` 而不是某个已声明场景的原因。

`results_scaling_r.csv` 用 **21 行覆盖 18 个（树规模, k）单元格**，成因相同但形式不
同：`tips = 100` 确实被测量了两次，一次在 `bd100`、一次在 `coal100`，而 Python 表里
只有 `bd100`。

## `source` 列的数据来源

两个计时驱动都会写出标明输入标签的 `source` 列，而且两边的拼写现在一致
（`benchmark_scaling.py` 过去写 `bd100__tree`，R 驱动写 `bd100`；两边都写 `bd100`）。
归档表早于这一列，所以此处记录的值是**按驱动的枚举顺序追溯归属**得到的，而这个顺序
可以精确复现归档表的行序：

| tips | `results_scaling_r.csv` 的 source | `results_scaling_py.csv` 的 source | 证据 |
|---|---|---|---|
| 1000 | `large1000`（第 1-3 行） | `large1000` | `benchmarks/large/large1000__tree.tre`，由 `make_large_inputs.py` 写出（`yule_tree(1000, age=1, seed=1100)`） |
| 2000 | `large2000`（第 4-6 行） | `large2000` | `benchmarks/large/large2000__tree.tre`，`seed=2100` |
| 200 | `large200`（第 7-9 行） | `large200` | `benchmarks/scaling/large200__tree.tre`，由 `make_scaling_inputs.py` 写出（`seed=300`） |
| 500 | `large500`（第 10-12 行） | `large500` | `benchmarks/scaling/large500__tree.tre`，`seed=600` |
| 100 | **`bd100`（第 13-15 行）与 `coal100`（第 19-21 行）** | `bd100`（第 4-6 行） | `validation/reference/bd100__tree.tre` 与 `coal100__tree.tre`，两者都是 100 tips / 198 条边 |
| 20 | `bd20`（第 16-18 行） | `bd20`（第 1-3 行） | `validation/reference/bd20__tree.tre` |

* R 驱动遍历 `sort(list.files(dirs, "__tree\\.tre$"))`，跨三个目录的完整路径顺序是
  `benchmarks/large/large1000`、`benchmarks/large/large2000`、
  `benchmarks/scaling/large200`、`benchmarks/scaling/large500`、
  `validation/reference/bd100`、`validation/reference/bd20`、
  `validation/reference/coal100`——归档表里 `tips` 序列 1000, 2000, 200, 500,
  **100, 20, 100** 正是这个顺序。
* Python 驱动按 `(边数, 文件名)` 排序，并且每个 (tips, k) 只保留一行。`bd100` 与
  `coal100` 同为 198 条边，于是按文件名的平局判定把 `bd100` 排在前面；归档表的 `tips`
  序列 20, 100, 200, 500, 1000, 2000 正是去掉重复档之后的该枚举。这个平局判定现在已经
  在 `benchmark_scaling.py` 里显式写成 `PUBLISHED_DUPLICATE`，不再是
  `drop_duplicates` 的副作用。

这里没有任何一档需要写成 `unknown`：结合驱动与磁盘上的树，每一档都可以归属确定。
`recompute_ratios.py` 仍然把 `unknown`（或缺失标签）当作硬错误，所以未来的表无法再
依赖这种宽容。

## 每条对外的数字来自哪张表

| 对外引用的结论 | 表配对 | 统计量 | 数值 |
|---|---|---|---|
| “76-556x”标题级加速比 | `results_scaling_py.csv` / `results_scaling_r.csv` | `slice_s`，全部 18 个单元格，100 tips 处 `source` = `bd100` | 75.6 - 556.0 |
| “100 切片时 148-302x”的固定 k 缩放 | 同一对表，`n_slices = 100`（6 个单元格） | `slice_s` | 148.2 - 302.2 |
| 最宽端点的四分位距系数 4.3 | `results_scaling_r.csv`，单元格 `source=large1000`、tips=1000、k=25 | `slice_q75 / slice_q25` | 4.29 |
| “6.7-8.1x”端到端速率拟合 | `results_py_large.csv` / `results_r_large.csv` | `cpd_s`，tips=1000/2000 | 6.72 - 8.06 |
| 峰值内存 1.93 / 3.49 / 6.69 MB 对 2.81 / 5.62 / 11.23 MB | `results_memory_py.csv` | `array_peak_bytes`、`copies_peak_bytes`，**MB = 10^6 字节** | 见下文“单位” |
| 不支撑任何对外结论 | `results_py.csv` / `results_r.csv` | `slice_s`，12 个一一对应的单元格（45-161x） | 仅报告，不引用 |
| 不支撑任何对外结论 | `results.csv` | Python 与自身基线比较 | 不是跨语言数字 |

`recompute_ratios.py` 会在每条结论前打印 `source: <Python 表> / <R 表>`，并在每一个
比值行后面再次点名两张表，因此上面这张映射表是被机器强制的，而不是一句口头承诺。

## 单位：`results_memory_py.csv`

`array_peak_bytes` 与 `copies_peak_bytes` 的单位是**字节**，即 `tracemalloc` 报告的
原始值。文档按 **MB = 10^6 字节**打印它们，这也正是 `benchmark_memory.py` 自己打印的
口径：

| tips | k | 数组，字节 | 数组 MB (10^6) | 数组 MiB (2^20) | 副本，字节 | 副本 MB (10^6) | 副本 MiB (2^20) |
|---|---|---|---|---|---|---|---|
| 2000 | 25 | 1933888 | 1.93 | 1.84 | 2808112 | 2.81 | 2.68 |
| 2000 | 50 | 3492618 | 3.49 | 3.33 | 5616224 | 5.62 | 5.36 |
| 2000 | 100 | 6691010 | 6.69 | 6.38 | 11232672 | 11.23 | 10.71 |

若按 MiB 口径读，同样三行是 1.84 / 3.33 / 6.38 与 2.68 / 5.36 / 10.71；对外引用的
1.93 / 3.49 / 6.69 与 2.81 / 5.62 / 11.23 是 10^6 口径的读法。`copies_over_array`
（1.45 / 1.61 / 1.68）与单位无关。

**`copies_peak_bytes` 是一份模拟，不是对 R 的测量。** 它是在本包的 Python 进程内追踪
参考设计的结构性代价——把树的边数组与长度数组复制 `k` 份。本目录里没有 R 的内存列，
也没有做过这样的测量；`benchmark_memory.py` 的文档字符串已经写明。这些行支持的结论
是“同一个 Python 进程里，一个 (E, k) 数组对 k 份副本”，即结构性比较，而不是跨语言
测量。

## 复现某一张表

```sh
python make_scaling_inputs.py && python benchmark_scaling.py
Rscript benchmark_scaling_r.R validation/reference benchmarks/large benchmarks/scaling \
  benchmarks/results_scaling_r.csv
python benchmark_phyloslicer.py                     # 验证网格，Python 侧
Rscript benchmark_r.R validation/reference benchmarks/results_r.csv
python make_large_inputs.py    && python benchmark_large.py
Rscript benchmark_large_r.R benchmarks/large benchmarks/results_r_large.csv
python benchmark_memory.py                          # Python，峰值内存
python recompute_ratios.py                          # 核对文档引用的每一个比值
```

`recompute_ratios.py` 既不需要 R，也不需要重新计时：它读取归档的表，重算文档所引用的
每一个区间（76-556、100 切片时的 148-302、速率拟合的 6.7-8.1，以及最宽那个端点的四分
位距），打印支撑每个单元格的表、场景与统计量；一旦某个被引用的区间不再能从表中推出、
某行的场景标签缺失或为 `unknown`、某张表的 `reps` 不再是文档声明的次数，或两张表对
“存在哪些单元格”的意见不一致，它就返回非零退出码。

## 必须与这些数字同行的限制

1. **R 的列是归档输出。** 它们是在更早的一次会话中、在 R 4.5.3 下用 treesliceR 1.1.0
   产出的；`results_scaling_r.csv` 在每个中位数旁边都存着 R 的四分位距。同一台机器上，
   R 侧重跑在不同会话之间最多能差到两倍，这就是每条对外引用的比值最多只保留两位有效
   数字的原因。
2. **会话漂移影响的是两个引擎，不只是 treesliceR。** 同一个 `E x k` 贡献矩阵调用，在
   2000 tips、k = 100 处被归档了两次：**Python 0.95158 ms**
   （`results_scaling_py.csv`）与 **1.4513 ms**（`results_py_large.csv`），对应的 R 分别
   是 140.98 ms 与 156.17 ms。于是缩放那一对给出 148.2x，端到端那一对给出 107.6x
   （1000 tips 处是 176.6x 对 104.7x）。两个数字都是同一操作在不同会话里的真实测量，
   任何一张表都不是另一张表的更正；这正是每条对外引用的比值都要点名其来源表的原因。
3. **`tips = 100` 有两棵输入树**（`bd100`，一棵出生-死亡树，以及 `coal100`，一棵溯祖
   树），因此 `results_scaling_r.csv` 对 18 个（树规模, k）单元格带着 21 行。比值是相对
   `bd100` 那一行计算的，与 Python 侧一致，并且这个选择是按 `source` 标签做的，不再按行
   序。该选择不改变被引用的区间——两种配对都给出总体 76-556、k = 100 时 148-302——但会
   改变这一档本身：100-tip 的比值在 `bd100` 下是 93.8 / 148.9 / 210.8x，而若把 Python 的
   `bd100` 行与 R 的 `coal100` 行相配则是 92.1 / 144.9 / 209.4x；k = 50 处比值对树规模的
   Spearman 相关分别是 +0.43 与 +0.54。`recompute_ratios.py` 会把两种读法都打印出来，只
   使用已声明的那一种，并在表失去让二者可区分的标签时直接失败。
4. **没有 BLAS/线程绑定。** 只有 `validation/compare.py` 会绑定
   `OMP/OPENBLAS/MKL/VECLIB` 线程；这些基准脚本不会，所以亚毫秒级的中位数会在最后一位
   数字上随运行而变化。
5. **峰值内存只在 Python 侧测量**，而且其中“副本”一侧是对参考设计的模拟，不是对 R 进程
   的测量（见上文“单位”）。
6. **原始重复向量只归档了一半。** `benchmark_scaling.py`、`benchmark_scaling_r.R`、
   `benchmark_phyloslicer.py` 与 `benchmark_r.R` 现在会写出每一次单独的计时运行
   （`repeat_times_s`，或当一行包含多个统计量时的 `repeat_slice_s` / `repeat_cpd_s` /
   `repeat_cpe_s`），因此中位数可以从它的输入重新算出。归档表早于这些列：它们每个单元
   格存的是中位数、`q25`、`q75`、最小值和 `reps`，没有重复向量。所以一个有据可查的比值
   只能审计到它的四分位数，审计不到它的逐次重复；这就是端点论断上残留的限制，
   `recompute_ratios.py` 在自己的输出里也这么写明。这里没有任何中位数是从未归档的重复
   里重新推导出来的。
7. **`results_r_large.csv` 里 `cpe_s = NA`**，而 `results_py_large.csv` 没有 `cpe_s` 列，
   因此在大树规模上不存在 CpE 的比较。

## 机器环境

归档的表是在一台 Apple silicon 笔记本上产出的。随包验证产物
（`validation/reference/manifest.txt`）所记录的 Python 环境是 Python 3.14.6，
配 numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3；R 一侧用的是 R 4.5.3
（aarch64-apple-darwin20.0.0），配 ape 5.8-1 与 treesliceR 1.1.0。
