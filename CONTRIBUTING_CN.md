# 贡献指南

本文件是 [`CONTRIBUTING.md`](CONTRIBUTING.md) 的中文译文，正文以英文文件为准。

感谢参与 PhyloSlicer。这是一个科研用软件包，因此一项变更的合格标准是它保持
**可核验**：文档中每个数字都能从仓库内的文件推出，每条"与参考实现一致"的论断
都能复现。

版本 0.1.0 · BSD-3-Clause · Python ≥ 3.10

## 搭建开发环境

```bash
git clone https://github.com/ZengZichao/PhyloSlicer.git
cd PhyloSlicer
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,io,spatial,accel]"
pre-commit install
```

`[accel]`（numba）与 `[spatial]`（geopandas）是可选后端。缺少后端时相关用例会
跳过，因此完整本地运行应两者都装；CI 故意让加速内核与纯 Python 内核互相对照。

## 变更必须通过的检查

```bash
ruff check src tests benchmarks examples testing_suite validation
ruff format --check src tests benchmarks examples testing_suite validation
python -m pytest -q
python testing_suite/scripts/run_all_tests.py
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python benchmarks/recompute_ratios.py
```

`.pre-commit-config.yaml` 在提交时执行两条 ruff 检查；CI 在 Ubuntu 与 macOS 上、
Python 3.10-3.13 全部版本上执行上述七项。

## 约定

- **版本号只有一个真源。** `src/phyloslicer/__init__.py` 中的 `__version__`
  即版本记录，`pyproject.toml` 动态读取它。不要再加第二处字面量。任何横幅或
  引用字段与它不一致，`tests/test_release_hygiene.py` 就会失败。
- **数值改动必须有参考基准。** 新增或改动的指数/速率路径必须与
  `validation/reference/manifest.txt` 中登记的 treesliceR 输出比对（SHA-256、
  形状、列名模式、存储朝向）。未登记的基准会被拒绝而不是被信任。
- **不要引用未经重新生成的基准数字。** `benchmarks/` 的四套协议彼此不可互换；
  先读 `benchmarks/README_CN.md`，如实填写 `source` 列，跨语言比较一律用中位数
  （不是最小值）。
- **诊断信息是 API 的一部分。** 拟合失败时返回 `converged=False` 且带可读
  `reason` 的行；不要把读者必须能看到的退化情形用异常吞掉。
- **文档成对交付。** 无后缀的英文文件是记录版本；每个 `X.md` 都有小节结构相同的
  `X_CN.md` 中文版。改动英文正文时，同一个 PR 里一并改中文。
- **不得出现个人路径。** 树中任何位置的绝对本地路径都会让
  `tests/test_release_hygiene.py::test_no_personal_paths_or_names_anywhere_in_the_tree`
  失败。

## 提出变更

凡涉及公开行为、指数定义或 `validation/` 中容差的改动，请先开 issue。PR 需要写
明你执行了哪些检查输出、关联的 issue，并附一个"没有该修复就会失败"的回归测试
（房内风格见 `tests/test_defect_regressions.py`：一个缺陷一条钉住测试，以症状
命名）。

## 发布

见 [RELEASE_CN.md](RELEASE_CN.md) · 英文：[RELEASE.md](RELEASE.md)。
