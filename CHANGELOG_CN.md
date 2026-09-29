# 变更日志

本文档记录 PhyloSlicer 的所有重要改动。版本记录以 `phyloslicer.__version__`
为准；下面的日期是该条目被打标签的日期。

## 0.1.0 — 2026-09-29

首个实现版本。

### 新增

- 切片表示为单张 `(E, k)` 区间相交贡献矩阵，因此不会为任何切片复制整棵树
  （`slicing/`）。
- 通过单次稀疏矩阵乘得到逐样区逐切片的 PD、PE 与 PB(+RW)（`indices/`）。
- 有界多起点最小二乘的"可认证"累积速率拟合，每行附带 `converged`、`r_squared`、
  `residuals` 与明确的 `reason`（`rates/`）。
- `origin_within_tree`：标记拟合速率推出的多样性起源时间早于树根的情形，并在
  `reason` 中记录 `origin extrapolated beyond root`（`rates/`）。
- 后验树不确定性：`posterior_rate`、`posterior_origin`、`bootstrap_sites`
  （`uncertainty/`）。
- 空间工作流助手：rook/queen/knn/radius 邻接构建器，以及接受 GeoDataFrame 的速率
  地图（`spatial/`、`viz/`）。
- `treesliceR` 风格的兼容层（`compat/`）与内置模拟数据集（`datasets/`）。
- 命令行界面：`phyloslicer slice | rates | posterior | validate`，附带 JSON 运行
  摘要（`cli/`）。
- `examples/quickstart.py`：一条命令即可在内置模拟数据上跑通完整流程；CI 会执行它，
  以免文档示例腐坏。
- 与 treesliceR v1.1.0 的交叉验证，覆盖 20/100 尖的生灭与凝聚场景，并带 SHA-256
  来源清单（`validation/`）。
- 基准测试阶梯：进程内朴素对照、跨语言切片内核（固定 `k`）、跨树规模的缩放曲线、
  1000/2000 尖端到端 CpD；另有 `benchmarks/benchmark_memory.py` 记录 Python 侧峰值
  内存，以及 `benchmarks/recompute_ratios.py` 从存档表格重算文档引用的每个加速比区间。
- 测试套件（`tests/`）与端到端测试方案/报告（`testing_suite/`）。
- `.pre-commit-config.yaml`（ruff check、ruff format），与 `dev` 可选依赖相匹配。
