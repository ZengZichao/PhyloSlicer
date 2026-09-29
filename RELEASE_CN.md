# 发布流程

本文件是 [`RELEASE.md`](RELEASE.md) 的中文译文，正文以英文文件为准。

PhyloSlicer 版本 0.1.0 · BSD-3-Clause

这里的"一次发布"指文档、引用元数据与存档验证基准在同一时刻描述的是同一个版本。
请按下面的顺序切版。

## 1. 固定版本记录

`src/phyloslicer/__init__.py` 是唯一允许出现版本字面量的位置，`pyproject.toml`
动态读取它。所有标注版本的文件都要同步更新：

- `README.md` 与 `README_CN.md` 的版本横幅
- `docs/USER_GUIDE.md` 与 `docs/USER_GUIDE_CN.md` 的版本横幅
- `CONTRIBUTING.md` 与 `CONTRIBUTING_CN.md` 的版本横幅
- `CHANGELOG.md` 与 `CHANGELOG_CN.md`——新增一条带日期的 `## X.Y.Z — YYYY-MM-DD` 条目
- `CITATION.cff`——`version:` 与 `date-released:`
- `.zenodo.json`——`version` 与 `publication_date`

上述任何一处发生漂移，`python -m pytest tests/test_release_hygiene.py` 就会失败。

## 2. 跑一遍 CI 的全部检查

```bash
ruff check src tests benchmarks examples testing_suite validation
ruff format --check src tests benchmarks examples testing_suite validation
python -m pytest -q
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python - <<'PY'
import hashlib, re
from pathlib import Path
man = Path("validation/reference/manifest.txt").read_text(encoding="utf-8")
entries = re.findall(r"\[file ([^\]]+)\]\nsha256 = ([0-9a-f]{64})", man)
bad = [n for n, h in entries
       if hashlib.sha256((Path("validation/reference") / n).read_bytes()).hexdigest() != h]
print(f"{len(entries)} registered references, {len(bad)} hash mismatches")
assert not bad and entries
PY
python benchmarks/recompute_ratios.py
```

`validation/compare.py` 会重写 `validation/reference/equivalence_report.csv`。
只有当这份报告确实应随本次发布移动时才提交它；而被钉住的基准文件本身，在不重跑
`validation/r/` 中 R 生成器的前提下绝不允许改动。

## 3. 打标签并发布

```bash
git tag -a "v0.1.0" -m "PhyloSlicer 0.1.0"
git push origin main --follow-tags
```

然后为该标签创建 GitHub Release，发布说明直接取 `CHANGELOG.md` 中对应小节。

## 4. 生成并记录 DOI

先在 <https://zenodo.org/account/settings/github/> 一次性开启 Zenodo 的 GitHub
集成，并打开本仓库。此后 Zenodo 会抓取打了标签的发布：`.zenodo.json` 用于预填
草稿，其中的 `version` 与 `publication_date` 会被标签元数据覆盖，这是正常现象。
**不要**在 `.zenodo.json` 里写 `conceptdoi` 字段——DOI 由 Zenodo 生成，提前写入会
造成冲突。

返回的 DOI 有两个：**版本 DOI**（每个发布各一个）与**与版本无关的概念 DOI**
（长期不变，引用它即引用所有版本）。要记录的是概念 DOI：

- 在 `CITATION.cff` 中新增 `identifiers:` 条目记录概念 DOI，并在
  `README.md` / `README_CN.md` 里与仓库地址并列引用
- 绝不要把占位 DOI 提交进仓库：一个解析不到任何东西的保留值比没有 DOI 更糟，
  因为会有人本着信任去引用它
- 若要挂版本徽章，用版本 DOI：`https://zenodo.org/badge/DOI/<version-doi>.svg`
- 这个占位符存在，只是为了让元数据在首个标签出现之前就能写全；带着它发布是不允许的

## 5. 在干净环境里确认

```bash
python -m venv /tmp/phyloslicer-check && source /tmp/phyloslicer-check/bin/activate
pip install "phyloslicer @ git+https://github.com/ZengZichao/PhyloSlicer@v0.1.0"
python -c "import phyloslicer; print(phyloslicer.__version__)"
phyloslicer --version
```

安装后的版本、标签、变更日志条目、引用元数据与 DOI 指向同一次构建时，本次发布才算完成。
