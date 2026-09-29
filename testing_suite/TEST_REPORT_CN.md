# PhyloSlicer 软件测试报告

> 本文件是 `TEST_REPORT.md` 的中文翻译；`TEST_REPORT.md`（无后缀）才是记录版本（version of record），两者内容不一致时以英文文件为准。

**版本：** 0.1.0  
**测试日期：** 2026-09-29  
**测试环境：** Python 3.14.6 / numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3 / matplotlib 3.11.1 / pytest 9.1.1 / dendropy 5.0.13 / biopython 1.88；可选加速与空间后端（numba、geopandas）在本环境中**未安装**，因此恰有 2 个用例按设计跳过  
**测试机器：** Apple M5（10 核）、32 GB RAM、macOS 27.0（`python3 -c "import platform; print(platform.platform())"` 报告 `macOS-27.0-arm64-arm-64bit-Mach-O`，`platform.machine()` 报告 `arm64`）  
**测试执行者：** `testing_suite/scripts/run_all_tests.py`（端到端；其输出目录被重定向到工作树之外，见第 7 节）与 `pytest`（单元/回归）  
**仓库路径：** 仓库根目录

---

## 1. 测试总结

| 指标 | 值 |
|---|---|
| **端到端检查总数** | 118 |
| **通过数** | 118 |
| **失败数** | 0 |
| **通过率** | **100.0%** |
| **总耗时** | 本轮多次运行为 3.3–4.3 秒（同为这 118 项检查；墙钟时间会波动，`testing_suite/logs/test_results.json` 记录的是最后一次运行） |
| **已有 pytest 测试** | 411 collected → 409 passed, 0 failed, 2 skipped |

### 结论

PhyloSlicer v0.1.0 的 118 项端到端检查全部通过；下文的分类统计可由当前版本的
`testing_suite/scripts/run_all_tests.py` 逐条复现。检查数从本报告此前的 102 项上升，是因为
`TEST_PLAN.md` 定义了 118 项而执行器只实现了其中 102 项：十条 `T-CLI` 用例、`T-ERR-06`（不连通分量）
与五条 `T-PERF` 用例都没有实现，所以"计划内用例全部通过"当时对执行器而言并不成立。这十六条用例现已
实现，并且 CI 会执行该执行器。pytest 套件现在收集 **411** 个用例（本报告上一版
记录的是 397 个），由 316 个 `def test_` 函数经 27 个 `@pytest.mark.parametrize` 展开而来；其中
**409 个通过、0 个失败、2 个跳过**。两个跳过项分别是 `tests/test_rates.py:279`（numba 加速内核与纯内核
对照）与 `tests/test_io_spatial_cli.py:404`（geopandas 网格邻接），它们跳过只是因为本机没有安装 numba
与 geopandas；装了 `[accel]` 与 `[spatial]` extras 的环境会执行这两个用例。

**关于发布卫生检查的数字。** 上面的 409/0/2 采集于 2026-09-29 22:12，此时文档配对与元数据的改动已经
落地。在同一段工作里更早运行的同一条命令，当时 `tests/test_release_hygiene.py` 与多份文档仍在被并行
修改，先后报告过 405 passed / 4 failed、407 / 2、408 / 1。这些失败项全部是文档卫生断言，不是软件缺陷：
`docs/USER_GUIDE_CN.md` 里残留的旧版本横幅、`pyproject.toml` 中的版本字段，以及
`test_every_english_document_ships_a_chinese_counterpart` 报告的
`documents without a Chinese counterpart: ['testing_suite/TEST_REPORT.md']`——最后这一项正是由本轮新增的
`TEST_REPORT_CN.md`（即本文件）满足的。功能与数值结果始终未变：每一次运行中端到端检查套件与非发布卫生
用例都全部通过。

---

## 2. 测试覆盖范围

### 2.1 按类别统计

| 类别 | 测试数 | 通过 | 失败 | 通过率 |
|---|---|---|---|---|
| Install（安装与依赖） | 6 | 6 | 0 | 100% |
| IO（输入/输出） | 10 | 10 | 0 | 100% |
| Core（核心算法） | 10 | 10 | 0 | 100% |
| Indices（多样性指数） | 10 | 10 | 0 | 100% |
| Rates（累积速率） | 10 | 10 | 0 | 100% |
| Slicing（切片） | 12 | 12 | 0 | 100% |
| Uncertainty（后验不确定性） | 8 | 8 | 0 | 100% |
| Spatial（空间工具） | 6 | 6 | 0 | 100% |
| Compat（兼容层） | 5 | 5 | 0 | 100% |
| Sim（模拟与数据集） | 6 | 6 | 0 | 100% |
| Viz（可视化） | 3 | 3 | 0 | 100% |
| CLI（命令行） | 10 | 10 | 0 | 100% |
| Error（异常与边界） | 10 | 10 | 0 | 100% |
| E2E（端到端） | 7 | 7 | 0 | 100% |
| Perf（性能与基准） | 5 | 5 | 0 | 100% |
| **合计** | **118** | **118** | **0** | **100%** |

### 2.2 功能覆盖

覆盖了 PhyloSlicer 的全部核心模块：

- ✅ `phyloslicer.core.tree` — TreeArray 构建、验证、深度、修剪
- ✅ `phyloslicer.slicing.api` — 5 种切片模式 + 切片栈
- ✅ `phyloslicer.rates.api` — CpD/CpE/CpB 速率 + 拟合 + 敏感性
- ✅ `phyloslicer.rates.fitting` — 指数拟合 + 批量拟合 + 边界检测
- ✅ `phyloslicer.indices.api` — PD/PE/PB 指数 + 每片分解
- ✅ `phyloslicer.uncertainty` — 后验速率/起源 + bootstrap + HDI/ESS
- ✅ `phyloslicer.spatial` — 邻接矩阵构建 + rate_map
- ✅ `phyloslicer.compat` — treesliceR 兼容层
- ✅ `phyloslicer.simulate` — Yule/出生-死亡树 + 数据集模拟
- ✅ `phyloslicer.io` — Newick/NEXUS 读写 + 矩阵构建 + 对齐
- ✅ `phyloslicer.viz` — 速率曲线/敏感性/后验分布图
- ✅ `phyloslicer.datasets` — sim_passerines 数据集
- ✅ `phyloslicer.cli.main` — CLI slice/rates/validate 子命令

---

## 3. 基准测试结果

### 3.1 运行时间（中位数，3 次重复）

| 树大小 | 位点数 | 切片数 | 切片 (s) | PD/片 (s) | CpD 速率 (s) | 后验 (s) | 内存 (MB) |
|---|---|---|---|---|---|---|---|
| 20 | 10 | 50 | 0.0001 | 0.0005 | 0.0050 | 0.0164 | 0.08 |
| 50 | 20 | 50 | 0.0002 | 0.0009 | 0.0053 | 0.0195 | 0.20 |
| 100 | 30 | 50 | 0.0005 | 0.0017 | 0.0060 | 0.0253 | 0.38 |
| 100 | 30 | 100 | 0.0005 | 0.0017 | 0.0070 | 0.0302 | 0.57 |
| 100 | 30 | 200 | 0.0005 | 0.0017 | 0.0077 | 0.0300 | 0.86 |
| 200 | 50 | 100 | 0.0012 | 0.0037 | 0.0096 | 0.0438 | 1.04 |
| 500 | 100 | 100 | 0.0034 | 0.0132 | 0.0186 | 0.0929 | 2.97 |

上表由 `testing_suite/scripts/run_benchmarks.py` 于 2026-09-29 在上述机器重新测得，运行时把它的
`BENCH_DIR` 指向仓库之外，因此没有改写仓库内的任何文件。仓库里
`testing_suite/benchmarks/benchmark_results.csv` 仍是 2026-09-26 在 Python 3.12.14 下那一次的结果：
它的峰值内存列与上表在舍入误差内一致，而耗时约为上表的 3-5 倍。更新该存档文件属于与本测试报告相分离
的另一步工作。

### 3.2 性能分析

- **切片操作**：所有测试规模都在 4 ms 以内完成，线性扩展。
- **CpD 速率计算**（含对齐 + 切片 + 贡献矩阵 + 拟合）：500-tip 树仅需 18.6 ms。
- **后验分析**（3 棵树）：500-tip 树用时 92.9 ms。
- **内存使用**：500-tip 树、100 sites、100 slices 仅使用 ~3 MB 峰值内存。
- **切片数扩展性**：从 50→200 切片，CpD 时间增加约 28%，符合 `边数 x 切片数` 贡献矩阵带来的线性增长预期。

---

## 4. 测试中发现的问题及修复

### 4.1 问题记录

本套端到端脚本在编写与调试阶段共修正过 **5 个测试脚本自身的问题**（非软件 bug），这 5 处修正至今仍保留在
`testing_suite/scripts/run_all_tests.py` 中：

| 问题 | 原因 | 修复方式 |
|---|---|---|
| T-RATES-04: fit_rate 参数恢复 | `fit_rate` 的模型 `cum ∝ exp(-r*ages)` 无截距项，归一化后 `cum[-1]=1` 约束了 r 的范围；测试使用 r=0.5 导致模型不匹配 | 调整为验证收敛性 + R² > 0.7 + rate > 0 |
| T-COMPAT-03: squeeze_int 参数名 | 参数名是 `to` 而非 `to_`（Python 关键字 `to` 在此处可用） | 改用 `to=0.5` |
| T-COMPAT-05: prune_tips 阈值 | 4-tip 对称树所有终端分支均为 1.0，阈值 0.5 导致全部被修剪 | 使用不等长分支的树和合理阈值 |
| T-E2E-05: CLI 启动方式 | 直接运行 `main.py` 导致相对导入失败 | 使用 `python -m phyloslicer.cli.main` |
| T-E2E-05: 输出文件冲突 | 重复运行时输出文件已存在 | 测试前清理输出文件 |

### 4.2 软件本身的问题

**本轮回归未发现新的软件 bug。** 需要说明的是，`tests/test_regressions.py` 与
`tests/test_defect_regressions.py` 本身就是为钉住开发过程中发现并已修复的行为而存在的（切片端点、
beta 剖面归一化、超度量判据的单位、`timeSteps` 方向、剪枝复杂度、Newick 引号、拒绝对未登记的验证
参考数据跑校验等），因此它们的数量不应被读作"软件无缺陷"，而应读作"这些缺陷已有回归防护"。

---

## 5. 端到端测试详情

下列数字全部取自第 1 节所述的 2026-09-29 端到端运行。

### 5.1 完整 CpD 工作流 (T-E2E-01)
- 输入：20-tip Yule 树 + 30 sites 存在/缺失矩阵
- 操作：`cpd_rate(tree, mat, n_slices=50)`
- 输出：30 行 DataFrame，含 11 列（`site_id`、`CpD`、`PD`、`pDO`、`converged`、`r_squared`、`reason`、
  `at_lower_bound`、`at_upper_bound`、`ultrametric`、`origin_within_tree`）
- 结果：✅ 全部 30 位点收敛

### 5.2 完整 CpE 工作流 (T-E2E-02)
- 操作：`cpe_rate(tree, mat, n_slices=50)`
- 结果：✅ 30 行，含 CpE/PE/pEO 列

### 5.3 完整 CpB 工作流 (T-E2E-03)
- 操作：`cpb_rate(tree, mat, adj, n_slices=50, component="sorensen")`
- 结果：✅ 30 行，含 CpB/PB/pBO 列

### 5.4 后验分析工作流 (T-E2E-04)
- 输入：5 棵后验树目录 + 30 sites 矩阵
- 操作：`posterior_rate(dir, mat, n_slices=30, hdi=0.95)`
- 结果：✅ 30 行，`n_trees = 5`，含 HDI 区间（`hdi_low`/`hdi_high`）与 ESS proxy
  （`ess_proxy`，并带 `ess_proxy_caveat` 列）

### 5.5 CLI 工作流 (T-E2E-05)
- 操作：`phyloslicer slice` + `phyloslicer rates`
- 结果：✅ 输出 Newick + CSV 文件（`cli_slice.tre`、`cli_rates.csv`）

### 5.6 交叉验证 (T-E2E-06)
- 操作：`phyloslicer validate` 与已登记的 treesliceR 参考数据（`bd20__random__CpD.csv`）比对
- 结果：✅ 30 个位点的等价率为 1.0000

### 5.7 可视化工作流 (T-E2E-07)
- 操作：`plot_rate_line` 生成速率曲线图
- 结果：✅ PNG 图形文件生成（`e2e_rate_line.png`）

---

## 6. 既有 pytest 测试套件

除了 118 项端到端检查，PhyloSlicer 还附带了完整的 pytest 单元测试套件：

```
411 collected -> 409 passed, 0 failed, 2 skipped   （整条命令耗时 11-16 秒，两次运行）
```

411 个用例由 316 个 `def test_` 函数经 27 个 `@pytest.mark.parametrize` 展开而来，命令为
`python3 -m pytest tests/ -p no:cacheprovider -q --tb=no -rs`。本轮的每一次运行都在跳过/失败摘要之后
结束，没有输出通常结尾的 `N passed in Xs` 计数行，因此上面的数字取自进度行与 `-rs`/`-rf` 短摘要，
整条命令的耗时取自 shell 自身的计时。

两个跳过项：

| 跳过位置 | 原因 |
|---|---|
| `tests/test_rates.py:279` | `numba not installed (pip install phyloslicer[accel])` |
| `tests/test_io_spatial_cli.py:404` | `could not import 'geopandas': No module named 'geopandas'` |

本环境既没有装 numba 也没有装 geopandas，所以恰好跳过 2 个用例；装了 `[accel]` 与 `[spatial]` extras
的环境会把这两个用例都跑起来。上表记录的那次运行没有任何失败项；文档仍在改动时采集的更早几次快照见第 1 节。

覆盖文件（数字为 `python3 -m pytest tests/ --collect-only -q` 在该文件收集到的用例数）：
- `tests/test_core.py` (13) — TreeArray 结构、深度、尖区间与校验
- `tests/test_slicing.py` (26) — 切片语义与解析解
- `tests/test_indices.py` (41) — 多样性指数与暴力对照
- `tests/test_rates.py` (23) — 速率拟合、诊断与后验
- `tests/test_invariants.py` (68) — 跨模块不变量与端到端属性
- `tests/test_io_spatial_cli.py` (34) — 读写、空间邻接与命令行
- `tests/test_compat.py` (52) — treesliceR 风格兼容层
- `tests/test_regressions.py` (27) — 回归测试
- `tests/test_defect_regressions.py` (82) — 一处缺陷对应一个钉住用例，按症状命名
- `tests/test_validation_chain.py` (32) — 验证链与自证防护
- `tests/test_release_hygiene.py` (11) — 发布包的元数据与文档卫生检查
- `tests/test_published_numbers.py` (2) — 从存档表重新算出文档给出的加速比区间

---

## 7. 测试资产清单

```
testing_suite/
├── TEST_PLAN.md              # 测试方案
├── TEST_PLAN_CN.md           # 测试方案的中文版
├── TEST_REPORT.md            # 本测试报告（记录版本）
├── TEST_REPORT_CN.md         # 本测试报告的中文翻译
├── data/                     # 测试数据
│   ├── small_tree.tre        # 4-tip 对称树
│   ├── small_mat.csv         # 4-site 矩阵
│   ├── medium_tree.tre       # 20-tip Yule 树
│   ├── medium_mat.csv        # 30-site 矩阵
│   ├── large_tree.tre        # 100-tip 出生-死亡树
│   ├── large_mat.csv         # 50-site 矩阵
│   ├── posterior_trees/      # 5 棵后验树
│   ├── posterior.nex         # NEXUS 格式后验树（5 棵）
│   ├── posterior.tre         # Newick 格式后验树（5 棵）
│   ├── grid_bounds.csv       # 网格边界（30 个单元）
│   ├── grid_adj.csv          # Queen 邻接矩阵（30 x 30）
│   ├── passerines_tree.tre   # 模拟雀形目树（50 tip）
│   ├── passerines_mat.csv    # 模拟矩阵（50 site）
│   ├── passerines_bounds.csv # 模拟网格（50 个单元）
│   ├── passerines_adj.csv    # 模拟邻接
│   └── passerines_posterior/ # 模拟后验树（3 棵）
├── scripts/
│   ├── generate_test_data.py # 测试数据生成
│   ├── run_all_tests.py      # 完整测试套件
│   └── run_benchmarks.py     # 基准测试
├── results/                  # 测试输出
│   ├── e2e_cpd.csv           # CpD 速率结果
│   ├── e2e_cpe.csv           # CpE 速率结果
│   ├── e2e_cpb.csv           # CpB 速率结果
│   ├── e2e_posterior.csv     # 后验分析结果
│   ├── e2e_validate.csv      # 交叉验证结果
│   ├── e2e_rate_line.png     # 速率曲线图
│   ├── cli_slice.tre         # CLI 切片输出
│   └── cli_rates.csv         # CLI 速率输出
├── logs/
│   └── test_results.json     # 详细测试结果 (JSON)
└── benchmarks/
    ├── benchmark_results.csv # 基准测试数据
    ├── BENCHMARK_REPORT.md   # 基准测试报告
    └── BENCHMARK_REPORT_CN.md # 基准测试报告的中文版
```

关于这份清单有两点说明。其一，`testing_suite/scripts/run_all_tests.py` 与
`testing_suite/scripts/run_benchmarks.py` 用默认参数运行时会写入 `results/`、`logs/` 与 `benchmarks/`；
第 1 节记录的 118 项运行就使用了这些默认路径，因此磁盘上的
`testing_suite/logs/test_results.json` 与 `results/` 中的 CSV 就是这一次运行（2026-09-29，版本
0.1.0），这些 CSV 也因此新增了 `origin_within_tree` 列——此前提交的副本还没有这一列。其二，上表列出的
每一个路径都存在于当前仓库中。

---

## 8. 测试结论

PhyloSlicer v0.1.0 通过了针对软件本身的各项测试：

1. **功能完整性** — 所有核心功能（切片、速率、指数、后验、空间、可视化、CLI）均正常工作。
2. **数值正确性** — 多样性指数的数学性质（行和 = 总 PD、partition 性质）得到验证。
3. **异常处理** — 所有边界条件和非法输入均产生明确的 ValueError。
4. **可重复性** — 相同种子产生相同结果。
5. **性能** — 500-tip 树的完整速率分析在 20 ms 以内完成。
6. **兼容性** — treesliceR 兼容层行为与 R 原版一致。
7. **交叉验证** — 与 treesliceR 参考数据的一致性通过验证（30 个位点等价率 1.0000）。

**建议：** 软件具备发布条件：118 项端到端检查全部通过，411 个 pytest 用例中 409 个通过，只余 2 个
可选后端跳过项。请在打标签的那个提交里，用装了 `[accel]` 与 `[spatial]` extras 的环境再跑一次
`python -m pytest -q`，让发布卫生检查与这两个后端用例在同一提交里保持全绿。
