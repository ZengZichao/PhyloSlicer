# 本目录中的基准测试表

中文译文，正文以 [README.md](README.md) 为准。

这里放着四套彼此独立的协议。它们**不可**互换，而且有两个列名在不同文件里指的是
不同的东西。文档中给出的每一条加速比都能追溯到下面四张表中的一张。

## 协议

| 文件 | 计时对象 | 统计量 | 重复次数 | 跨语言？ |
|---|---|---|---|---|
| `results.csv` | `bench_core.py`：向量化内核 vs 进程内朴素基线 | **最小值** | 3（朴素基线：1） | **否** - Python vs Python |
| `results_py.csv`, `results_r.csv` | 验证场景输入上的切片内核 + CpD/CpE 拟合 | 中位数 | 5 | 是 |
| `results_scaling_py.csv`, `results_scaling_r.csv` | **仅切片内核**（`slice_pieces(...).contribution_matrix()` vs `phylo_pieces(...)`），每个树规模在 k = 25/50/100 的每一档上 | 中位数 | 15（+1 次被丢弃的热身运行） | 是 |
| `results_py_large.csv`, `results_r_large.csv` | 250 个组合、1000/2000 tips、k = 100 上的端到端 CpD | 中位数 | 5 | 是 |

`slice_s` 在 `results_py.csv`、`results_scaling_py.csv` 和 `results_py_large.csv`
里指的是 **E x k 贡献矩阵**，但在 `testing_suite/benchmarks/benchmark_results.csv`
里只指**窗口构建本身**。不要把这两个目录里的数字互相比较。

## 复现某一张表

```sh
python make_scaling_inputs.py && python benchmark_scaling.py       # Python，缩放
Rscript benchmark_scaling_r.R                                      # R，缩放
python make_large_inputs.py    && python benchmark_large.py        # Python，端到端
Rscript benchmark_large_r.R                                        # R，端到端
python benchmark_memory.py                                         # Python，峰值内存
python recompute_ratios.py                                         # 核对文档引用的每一个比值
```

`recompute_ratios.py` 既不需要 R，也不需要重新计时：它读取归档的表，重算文档所
引用的每一个区间（76–556、100 切片时的 148–302、速率拟合的 6.7–8.1，以及最宽
那个端点的四分位距），打印支撑它们的单元格；一旦某个被引用的区间不再能从表中
推出，它就返回非零退出码。

两侧读取的都是同一批 `*__tree.tre` 输入，两个写出程序都会输出一列 `source`，
写明输入标签，因此任何一行都可以追溯到它是在哪一棵具体的树上测得的。

**归档的 `results_scaling_*.csv` 表早于这一列。** 它们是更早的输出，所以其中重复
出现的 100-tip 那一档其实是同一个网格点上两次没有标签的运行。没有任何一条对外
给出的论断依赖这两行谁是谁，也正是出于这个原因，每个比值最多只引用到两位有效
数字。要恢复数据来源，请重新运行：

```bash
python make_scaling_inputs.py && python benchmark_scaling.py
Rscript benchmark_scaling_r.R
```

## 必须与这些数字同行的限制

1. **R 的列是归档输出。** 它们是在更早的一次会话中、在 R 4.5.3 下用
   treesliceR 1.1.0 产出的；`results_scaling_r.csv` 在每个中位数旁边都存着 R 的
   四分位距。同一台机器上，R 侧重跑在不同会话之间最多能差到两倍，这就是每条
   对外引用的比值最多只保留两位有效数字的原因。
2. **`tips = 100` 有两棵输入树**（`bd100`，一棵出生-死亡树，以及 `coal100`，
   一棵溯祖树）。因此 `results_scaling_r.csv` 对 18 个（树规模, k）单元格带着 21
   行。比值是相对 `bd100` 那一行计算的，与 Python 侧一致。
3. **没有 BLAS/线程绑定。** 只有 `validation/compare.py` 会绑定
   `OMP/OPENBLAS/MKL/VECLIB` 线程；这些基准脚本不会，所以亚毫秒级的中位数会在
   最后一位数字上随运行而变化。
4. **峰值内存只在 Python 侧测量。** `benchmark_memory.py` 在
   `results_memory_py.csv` 里按树规模记录 `E x k` 贡献矩阵的 tracemalloc 峰值，
   并与同样数目的序列化 `TreeArray` 副本相对照。参考实现一侧不做测量，因为
   "每切片复制一次树"是已发布的 R 实现本身的属性，而不是要拿它来计时的对象；
   这里给出的比较是本包内部的"数组 vs 副本"比较。
5. **原始重复向量只归档了一半。** `benchmark_scaling.py` 和
   `benchmark_phyloslicer.py` 现在会写一列 `repeat_times_s`，保存每一次单独的
   运行，因此中位数可以从它的输入重新算出。R 驱动每个单元格只写中位数、`q25`、
   `q75` 和最小值，不写重复向量，而归档的 Python 表早于这个新列。所以一个有据
   可查的比值只能审计到它的四分位数，审计不到它的逐次重复；这就是端点论断上
   残留的限制，`recompute_ratios.py` 在自己的输出里也这么写明。
6. **`results_r_large.csv` 里 `cpe_s = NA`**，而 `results_py_large.csv` 没有
   `cpe_s` 列，因此在大树规模上不存在 CpE 的比较。

## 机器环境

归档的表是在一台 Apple silicon 笔记本上产出的。随包验证产物
（`validation/reference/manifest.txt`）所记录的 Python 环境是 Python 3.14.6，
配 numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3；R 一侧用的是 R 4.5.3
（aarch64-apple-darwin20.0.0），配 ape 5.8.1 与 treesliceR 1.1.0。
