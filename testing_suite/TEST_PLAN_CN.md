# PhyloSlicer 软件测试方案

本文档是 [`TEST_PLAN.md`](TEST_PLAN.md) 的中文译文；无后缀的 `TEST_PLAN.md` 是记录版本（version of record），以英文正文为准。

**版本：** 0.1.0  
**日期：** 2026-09-26  
**环境：** conda `phyloslicer` (Python 3.12, numpy/scipy/pandas/matplotlib/numba)  
**工作路径：** `testing_suite/`（相对仓库根目录）

---

## 1. 测试目标

对 PhyloSlicer 进行完整、系统的端到端测试，覆盖：

1. **安装与依赖** — 确认包可安装、导入，所有核心/可选依赖可用。
2. **命令行接口 (CLI)** — `phyloslicer slice | rates | posterior | validate` 四个子命令的端到端测试。
3. **输入解析** — Newick/NEXUS 读取、存在/缺失矩阵构建、树-矩阵对齐。
4. **核心算法** — 贡献矩阵、切片窗口、PD/PE/PB 指数、速率拟合。
5. **输出结果** — CSV/Newick 文件正确性、字段完整性、数值一致性。
6. **异常处理** — 非法输入、边界条件、退化场景的错误信息。
7. **边界条件** — 2-tip 树、空位点、退化邻域、NaN 处理。
8. **真实数据全流程** — 使用模拟数据集完成从树读取到速率计算、后验分析、可视化的完整工作流。
9. **性能/基准测试** — 不同规模数据集的运行时间和内存。
10. **可重复性验证** — 相同随机种子结果一致性、数值稳定性。

---

## 2. 测试环境

| 项目 | 版本 |
|---|---|
| Python | 3.12 |
| numpy | >= 1.24 |
| scipy | >= 1.11 |
| pandas | >= 2.0 |
| matplotlib | >= 3.7 |
| numba | >= 0.59 (可选加速) |
| dendropy | >= 5.0 (I/O 适配器) |
| pytest | >= 8.0 |

---

## 3. 测试用例分类

### 3.1 安装与依赖测试 (T-INSTALL)

| ID | 描述 | 预期 |
|---|---|---|
| T-INSTALL-01 | `import phyloslicer` 成功 | 无 ImportError |
| T-INSTALL-02 | `phyloslicer.__version__` 返回字符串 | 非空字符串 |
| T-INSTALL-03 | 核心依赖 numpy/scipy/pandas/matplotlib 可导入 | 无异常 |
| T-INSTALL-04 | 可选依赖 numba 可导入 | 无异常 |
| T-INSTALL-05 | 可选依赖 dendropy 可导入 | 无异常 |
| T-INSTALL-06 | CLI 入口 `phyloslicer --version` 正常输出 | exit 0 |

### 3.2 命令行接口测试 (T-CLI)

| ID | 描述 | 预期 |
|---|---|---|
| T-CLI-01 | `phyloslicer slice tree.tre --rootward 5.0 --out out.tre` | exit 0, 输出 Newick 文件 |
| T-CLI-02 | `phyloslicer slice tree.tre --tipward 5.0 --out out.tre` | exit 0 |
| T-CLI-03 | `phyloslicer slice tree.tre --n-slices 10 --out dir/` | exit 0, 10 个 .tre 文件 |
| T-CLI-04 | `phyloslicer rates tree.tre sites.csv --rate cpd --out r.csv` | exit 0, CSV 输出 |
| T-CLI-05 | `phyloslicer rates tree.tre sites.csv --rate cpe --out r.csv` | exit 0 |
| T-CLI-06 | `phyloslicer posterior trees/ sites.csv --out p.csv` | exit 0 |
| T-CLI-07 | `phyloslicer validate tree.tre sites.csv --reference ref.csv --out v.csv` | exit 0 (等价率=1) |
| T-CLI-08 | `phyloslicer slice` 缺少必需参数 | exit 2 |
| T-CLI-09 | `phyloslicer slice tree.tre --rootward 5.0 --tipward 3.0 --out o.tre` | exit 2 (互斥参数) |
| T-CLI-10 | `--log` 生成 JSON 日志 | 日志文件包含 version/elapsed_s |

### 3.3 输入解析测试 (T-IO)

| ID | 描述 | 预期 |
|---|---|---|
| T-IO-01 | Newick 字符串解析 | 正确的 tip_labels/edges/lengths |
| T-IO-02 | Newick 文件读取 | 与字符串解析一致 |
| T-IO-03 | NEXUS 多树文件读取 | 返回所有树 |
| T-IO-04 | 带方括号注释的 Newick | 注释被正确剥离 |
| T-IO-05 | 带引号标签的 Newick | 引号标签正确解析 |
| T-IO-06 | Newick 往返 (parse → write → parse) | 结构一致 |
| T-IO-07 | 长格式表转存在/缺失矩阵 | 正确的 sites×species DataFrame |
| T-IO-08 | 宽格式表转存在/缺失矩阵 | 正确转换 |
| T-IO-09 | align_tree_matrix 树-矩阵对齐 | 正确修剪/丢弃 + 警告 |
| T-IO-10 | 树和矩阵无公共物种 | ValueError |

### 3.4 核心算法测试 (T-CORE)

| ID | 描述 | 预期 |
|---|---|---|
| T-CORE-01 | TreeArray.from_newick 构建正确 | n_tips/n_edges/root_age 正确 |
| T-CORE-02 | 节点深度计算 | 每个节点深度 = 父深度 + 分支长度 |
| T-CORE-03 | tip_intervals 正确性 | 每条边的 [lo, hi) 覆盖正确的 tip 范围 |
| T-CORE-04 | 贡献矩阵 C[e,j] | 等于 max(0, min(d_e, t_j) - max(b_e, t_{j-1})) |
| T-CORE-05 | 切片窗口 (equal-time) | 窗口等宽，覆盖 [0, T] |
| T-CORE-06 | 切片窗口 (equal-PD) | 每片 PD 相等 |
| T-CORE-07 | subarray 修剪 | 保留 MRCA 子树 |
| T-CORE-08 | ultrametric 检测 | 超度量树返回 True |
| T-CORE-09 | 非超度量树检测 | 返回 False + 警告/错误 |
| T-CORE-10 | total_pd 计算 | 等于所有分支长度之和 |

### 3.5 多样性指数测试 (T-INDICES)

| ID | 描述 | 预期 |
|---|---|---|
| T-INDICES-01 | pd() 总 PD 计算正确 | 与手算一致 |
| T-INDICES-02 | pe() 系统发育特有性计算正确 | 与手算一致 |
| T-INDICES-03 | pb() 全树 beta 多样性 | 值在 [0, 1] |
| T-INDICES-04 | pd_per_slice 每片 PD | 行和等于总 PD |
| T-INDICES-05 | pe_per_slice 每片 PE | 行和等于总 PE |
| T-INDICES-06 | pb_per_slice 行和为 1 (partition) | np.allclose(row_sums, 1.0) |
| T-INDICES-07 | pb_per_slice sorensen/turnover/nestedness | 所有 component 正确 |
| T-INDICES-08 | pb_per_slice multisite/pairwise | 两种 approach 正确 |
| T-INDICES-09 | multisite_domain mixed vs paired | 值符合预期 |
| T-INDICES-10 | 退化邻域处理 | converged=False + reason |

### 3.6 累积速率测试 (T-RATES)

| ID | 描述 | 预期 |
|---|---|---|
| T-RATES-01 | cpd_rate 返回正确 DataFrame | 含 site_id/CpD/PD/pDO/converged 等列 |
| T-RATES-02 | cpe_rate 返回正确 DataFrame | 含 CpE/PE/pEO |
| T-RATES-03 | cpb_rate 返回正确 DataFrame | 含 CpB/PB/pBO |
| T-RATES-04 | fit_rate 对已知指数曲线恢复参数 | r 误差 < 1% |
| T-RATES-05 | origin_time 计算 | pXO = -ln(1-p)/r |
| T-RATES-06 | 空位点处理 | converged=False, reason 非空 |
| T-RATES-07 | 边界拟合标记 | at_lower_bound/at_upper_bound |
| T-RATES-08 | sensitivity 切片数敏感性 | 返回 suggested_slices |
| T-RATES-09 | rate_profile (normalize=True) | 行和为 1 |
| T-RATES-10 | rate_profile (normalize=False) | 行和等于总 PD |

### 3.7 切片测试 (T-SLICING)

| ID | 描述 | 预期 |
|---|---|---|
| T-SLICING-01 | slice_rootward(time=0) → 空切片 | PD=0 |
| T-SLICING-02 | slice_rootward(time=T) → 全树 | PD=total_pd |
| T-SLICING-03 | slice_tipward(time=0) → 全树 | PD=total_pd |
| T-SLICING-04 | slice_tipward(time=T) → 空切片 | PD=0 |
| T-SLICING-05 | slice_interval(start, stop) | 正确区间 |
| T-SLICING-06 | slice_interval invert=True | 返回 (older, younger) 元组 |
| T-SLICING-07 | slice_pieces(n=10) → 10 片 | n_slices=10 |
| T-SLICING-08 | slice_pieces(width=0.5) | 正确宽度 |
| T-SLICING-09 | slice_pieces criterion="pd" | 等PD切片 |
| T-SLICING-10 | prune_tips 正确修剪 | 保留/丢弃正确 tips |
| T-SLICING-11 | 时间超过 root_age → ValueError | 异常信息明确 |
| T-SLICING-12 | 负时间 → ValueError | 异常信息明确 |

### 3.8 后验不确定性测试 (T-UNCERTAINTY)

| ID | 描述 | 预期 |
|---|---|---|
| T-UNC-01 | posterior_rate 目录输入 | 正确的 per-site 汇总 |
| T-UNC-02 | posterior_rate Newick 字符串输入 | 同上 |
| T-UNC-03 | posterior_origin 计算 | 含 origin_median/hdi_low/high |
| T-UNC-04 | bootstrap_sites | 返回 n_boot 行 |
| T-UNC-05 | iter_trees 多格式输入 | 正确迭代所有树 |
| T-UNC-06 | 坏树跳过 + 警告 | 不中断, n_trees 减少 |
| T-UNC-07 | HDI 区间计算 | hdi_low <= rate_median <= hdi_high |
| T-UNC-08 | ESS proxy 计算 | ess_proxy <= n_trees |

### 3.9 空间工具测试 (T-SPATIAL)

| ID | 描述 | 预期 |
|---|---|---|
| T-SPATIAL-01 | adjacency_from_bounds queen | 对角线=1, 角接触=1 |
| T-SPATIAL-02 | adjacency_from_bounds rook | 角接触=0 |
| T-SPATIAL-03 | adjacency_from_coords knn | 对称, 对角线=1 |
| T-SPATIAL-04 | adjacency_from_coords radius | 对称 |
| T-SPATIAL-05 | to/from_adjacency_dict 往返 | 一致 |
| T-SPATIAL-06 | rate_map 生成图形 | 返回 matplotlib Axes |

### 3.10 兼容层测试 (T-COMPAT)

| ID | 描述 | 预期 |
|---|---|---|
| T-COMPAT-01 | squeeze_root == slice_rootward | 结果一致 |
| T-COMPAT-02 | squeeze_tips == slice_tipward | 结果一致 |
| T-COMPAT-03 | squeeze_int invert 返回顺序 | (younger, older) R 风格 |
| T-COMPAT-04 | phylo_pieces timeSteps | 递增深度向量 |
| T-COMPAT-05 | prune_tips method=1/2 | side after/before |

### 3.11 模拟与数据集测试 (T-SIM)

| ID | 描述 | 预期 |
|---|---|---|
| T-SIM-01 | yule_tree 生成正确 | n_tips 正确, ultrametric |
| T-SIM-02 | birth_death_tree 生成正确 | n_tips 正确, ultrametric |
| T-SIM-03 | make_dataset random | 正确维度矩阵 |
| T-SIM-04 | make_dataset clustered | 空间自相关 |
| T-SIM-05 | sim_passerines 数据集 | 含 trees/mat/coords/bounds/adj |
| T-SIM-06 | 种子可重复性 | 相同种子相同结果 |

### 3.12 可视化测试 (T-VIZ)

| ID | 描述 | 预期 |
|---|---|---|
| T-VIZ-01 | plot_rate_line | 返回 Axes |
| T-VIZ-02 | plot_sensitivity | 返回 Axes |
| T-VIZ-03 | plot_posterior | 返回 Axes |

### 3.13 真实数据端到端测试 (T-E2E)

| ID | 描述 | 预期 |
|---|---|---|
| T-E2E-01 | 完整 CpD 工作流 (API) | 有意义的速率结果 |
| T-E2E-02 | 完整 CpE 工作流 (API) | 有意义的速率结果 |
| T-E2E-03 | 完整 CpB 工作流 (API) | 有意义的速率结果 |
| T-E2E-04 | 完整后验分析工作流 (API) | 有意义的 HDI 区间 |
| T-E2E-05 | 完整 CLI 工作流 | 所有子命令正常运行 |
| T-E2E-06 | 交叉验证 (validate) | 等价率 = 1.0 |
| T-E2E-07 | 可视化全流程 | 图形文件生成 |

### 3.14 异常与边界测试 (T-ERROR)

| ID | 描述 | 预期 |
|---|---|---|
| T-ERR-01 | 1-tip 树 → ValueError | 明确错误信息 |
| T-ERR-02 | NaN 分支长度 → ValueError | 明确错误信息 |
| T-ERR-03 | 负分支长度 → ValueError | 明确错误信息 |
| T-ERR-04 | 重复 tip 标签 → ValueError | 明确错误信息 |
| T-ERR-05 | 环 → ValueError | 明确错误信息 |
| T-ERR-06 | 不连通树 → ValueError | 明确错误信息 |
| T-ERR-07 | n_slices < 1 → ValueError | 明确错误信息 |
| T-ERR-08 | percentile 超出 (0,1) → ValueError | 明确错误信息 |
| T-ERR-09 | hdi 超出 (0,1] → ValueError | 明确错误信息 |
| T-ERR-10 | 无效 criterion → ValueError | 明确错误信息 |

### 3.15 性能与基准测试 (T-PERF)

| ID | 描述 | 预期 |
|---|---|---|
| T-PERF-01 | 100-tip 树 CpD 速率 | < 5s |
| T-PERF-02 | 1000-tip 树 CpD 速率 | < 60s |
| T-PERF-03 | 切片数扩展性 (10→200) | 线性时间 |
| T-PERF-04 | 内存使用合理 | 无内存泄漏 |
| T-PERF-05 | 与内置基准 (`benchmarks/`) 对比 | 结果可复现 |

---

## 4. 测试数据

### 4.1 模拟数据

- **小树**: 4-tip 对称树 `((A:1,B:1):1,(C:1,D:1):1);`
- **中等树**: 20-tip Yule 树 (seed=42, age=10.0)
- **大树**: 100-tip 出生-死亡树 (seed=42, age=10.0)
- **后验树集**: 5 棵 50-tip Yule 树 (不同 seed)
- **存在/缺失矩阵**: 30 sites × 20 species (seed=7)
- **邻接矩阵**: queen 邻接 (6×6 网格，取前 30 个单元)

### 4.2 已有验证参考数据

使用 `validation/reference/` 中的 treesliceR 参考文件进行交叉验证。

---

## 5. 测试执行流程

```
1. 生成测试数据 (testing_suite/scripts/generate_test_data.py)
2. 运行单元测试 (pytest，tests/)
3. 运行端到端测试 (testing_suite/scripts/run_all_tests.py)
4. 运行基准测试 (testing_suite/scripts/run_benchmarks.py)
5. 汇总结果 → 测试报告 (testing_suite/logs/test_results.json → testing_suite/TEST_REPORT.md)
```

---

## 6. 通过/失败标准

- **通过**: 所有断言成立，无异常，输出文件正确。
- **失败**: 任何断言失败、异常未捕获、输出缺失/错误。
- **整体**: 100% 通过方可进入仓库整理阶段。
