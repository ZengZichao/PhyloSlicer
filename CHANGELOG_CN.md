# 变更日志

本文档记录 PhyloSlicer 的所有重要改动。版本记录以 `phyloslicer.__version__`
为准；下面的日期是该版本公开的日期，更早的源码定稿日期写在条目正文里。

## 0.1.1 — 2026-10-18

首次发布到 PyPI。0.1.0 谱系按 2026-10-08 审阅意见加固：修复了唯一的库级缺陷，
扩充验证套件以钉住修正后的行为，并让发布管线获得 PyPI 发布能力。源码于
2026-10-18 定稿并公开。

### 新增

- 通过 GitHub Actions Trusted Publishing 发布到 PyPI
  （`.github/workflows/publish.yml`）：推送 `v*` 标签即构建 sdist 与 wheel、
  校验元数据并以 OIDC 上传，仓库中不保存任何 API token。
- Python 3.14 进入 CI 矩阵与 classifiers。
- 钉住修正后行为的回归测试：混合超度量池下的后验审计列、公布数字、输出模式、
  CI 镜像与审阅包；β 剖面划分不变式的覆盖从 6 个组合补到 24 个。
- 审阅包守卫：DOI 形状与 AST 片段两级校验、显式的 `--allow-no-git` /
  `--allow-dirty` 并记录脏路径、身份字面量拆分至 `tools/identity_literals.py`。

### 修复

- `posterior_rate` / `posterior_origin` 的审计列算在了错误的轴上（`np.vstack`
  把抽样拉平后才按站点索引）；现改为 `np.stack` 成 `(n_draws, n_sites, 2)` 后
  按该形状取列。
- `validation/compare.py` 默认拒绝覆写受版本控制的 `equivalence_report.csv`
  （新增 `--in-place`），且报告新增 `n_finite` 列，使正文引用的
  "35 950 / 36 000" 可直接从交付表重算。
- CI 日常维护加固：防双跑、CodeQL、构建校验与 ruff 钉版保护；Dependabot 覆盖
  GitHub Actions 生态并附 major 屏蔽规则。
- `LICENSE` 补齐第二作者版权人；`NOTICE` 与 README 声明英文文档为权威版、
  中文为同结构翻译。

## 0.1.0 — 2026-09-30

首次实现。0.1.0 源码于 2026-09-29 定稿；公开的 GitHub Release 与 Zenodo
存档均于 2026-09-30 建立，这一天才是该版本的发布日期，也是 Software Notes
规则所要求记录的日期。

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
- 手册之外的发布文档集：`docs/CI_WORKFLOWS*.md`（CI 工作流定义的逐字镜像——因为
  `.github/` 与 `tools/` 无法进入审稿代码包）、`docs/WORD_COUNT*.md`（让稿件篇幅统计
  可复现的计数口径），以及 `docs/ZENODO_DEPOSIT_CHECKLIST*.md`（存档记录的逐字段更正
  清单）。每一份英文文档都配一份中文对应件，该规则由
  `tests/test_release_hygiene.py` 强制校验。
