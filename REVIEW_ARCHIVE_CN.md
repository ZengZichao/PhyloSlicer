# 匿名化同行评审存档包

英文原版：[REVIEW_ARCHIVE.md](REVIEW_ARCHIVE.md)（本节内容与英文版逐节对应）。

Ecography《Software Notes》栏目允许投稿以补充压缩包的形式送审，而不必借助匿名仓库
镜像：“It can be made available to reviewers either via a service such as Anonymous
Github … or by submitting the full source-code package and data as a supplementary
zip file.” 本文说明生成该压缩包的打包工具、它移除了什么、刻意保留了什么，以及如何
验证结果。

投稿清单中由此得到支撑的陈述只有一句、且可检验：**CI 通过打包工具构建评审存档，
只要残留任何个人路径或作者身份，该任务即失败。** 该陈述由
`tools/make_review_archive.py`（构建后重新读取并扫描自身产物）与
`tests/test_review_archive.py`（解包并实际运行）共同保证，并作为具名步骤存在于
`.github/workflows/tests.yml` 中。

## 构建

```bash
python tools/make_review_archive.py                      # -> dist/review_archive/PhyloSlicer_review_<version>.zip
python tools/make_review_archive.py --format tar.gz
python tools/make_review_archive.py --out /tmp/review --verbose
python tools/make_review_archive.py --allow-dirty        # 仅当允许清单登记未提交路径时
python tools/make_review_archive.py --allow-no-git       # 仅在确实没有 .git 处；记录 commit = null
```

只用标准库：无需安装、无需 `PYTHONPATH`——版本号是从
`src/phyloslicer/__init__.py` 正则读取而非 import，因此工具不引入任何运行时依赖。
输出目录已被 git 忽略（`.gitignore` 同时列出 `dist/review_archive/` 与
`review_archive/`）。

退出码：守护检查通过为 `0`；残留身份或构建产物为 `1`；存档无法作出诚实的来源声明
（读不到 git 元数据，或工作树有未提交改动）而又未给出对应的 `--allow-*` 开关时，同样
以 `PROVENANCE GATE` 信息退出 `1`；指定的存档文件不存在为 `2`。

## 验证

```bash
# 不重新构建，直接复核已有压缩包
python tools/make_review_archive.py --check-only dist/review_archive/PhyloSlicer_review_0.1.0.zip

# 完整契约，其中包含“反向用例”
python -m pytest -q tests/test_review_archive.py
```

`--check-only` 会重新读取每个条目、扫描身份模式、把每个条目的 SHA-256 与来源清单
逐条比对，并拒绝任何构建/VCS 产物。测试中有两个用例故意在正常压缩包内植入作者姓名
和一个 `__pycache__` 文件，并断言守护检查**必须失败**——一个不可能失败的检查毫无价值。

确定性：条目排序、所有时间戳统一取 HEAD 的提交日期（可用 `SOURCE_DATE_EPOCH`
覆盖）、归档元数据归一化（权限 0644、uid/gid 为 0、不写所有者名）、gzip 以
`mtime=0` 写出。连续两次构建产物逐字节相同，`test_two_builds_are_byte_identical`
对此作了断言，实测哈希为：

| 存档 | SHA-256 |
| --- | --- |
| `PhyloSlicer_review_0.1.0.zip` | `d3750aa2e2abd553cf19fc3343a9440d270b953d01dea7dfa8e22a61b81d1e79` |
| `PhyloSlicer_review_0.1.0.tar.gz` | `0af47241b70ed8069f895fc774870495d3806a2181a6fd1004f3b7e1198dd59a` |

（均基于提交 `0cba60f`、工作树干净时用 CPython 3.11 构建三次，三份产物逐字节相同；存档含 263 个条目。
仍欢迎自行重建后比对——印出哈希的意义就在于它必须能被复现。）`SOURCE_DATE_EPOCH` 会被采用，也是唯一
改变时间戳的开关。）

解包后即是一份可直接工作的检出：

```bash
unzip PhyloSlicer_review_0.1.0.zip && cd PhyloSlicer_review_0.1.0
python -m pytest -q
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python -m pip install -e .
```

## 收录范围与排除清单

收录：`pyproject.toml`、`MANIFEST.in`、`.gitignore`、`LICENSE`、`CITATION.cff`、
`.zenodo.json`、`src/`、`tests/`、`examples/`、`docs/`、`benchmarks/`、`validation/`
（含参考 CSV、`manifest.txt` 与随包发布的 `equivalence_report.csv`）、
`testing_suite/`，以及全部根目录文档——`README`、`CHANGELOG`、`CONTRIBUTING`、
`RELEASE`、`CODE_OF_CONDUCT`、`SECURITY`、`NOTICE`。

这份清单由两条规则而非主观偏好决定：

* 评审者必须能够安装、运行并独立复核该软件包；
* 软件包自身的卫生测试必须在存档内通过，因此这些测试读取的文件必须在包内——
  `test_version_licence_and_release_date_are_stated` 会读取 `CITATION.cff`、
  `.zenodo.json` 与十份手册/变更日志，
  `test_gitignore_excludes_artifacts_that_carry_absolute_paths` 会读取 `.gitignore` 本身。

排除（原因同时写入来源清单，并在 `--verbose` 下逐条打印）：`.git/`、`.github/`
（其中 `CODEOWNERS` 把每个文件的所有者声明为 `@ZengZichao`）、
`.pre-commit-config.yaml`、`__pycache__/`、`.pytest_cache/`、`.ruff_cache/`、
`.mypy_cache/`、`*.egg-info/`、`build/`、`dist/`、`.DS_Store`、`*.pyc`、输出目录本身，
以及关键的 `tools/` 与 `REVIEW_ARCHIVE.md`/`REVIEW_ARCHIVE_CN.md`。打包工具的替换表
和本文档必然包含它们所删除的字面量，收录进来等于自行拆掉匿名化，而运行软件包并不需要
它们。替代做法是：包内由打包工具生成
`PEER_REVIEW_NOTICE.md` / `PEER_REVIEW_NOTICE_CN.md`，向评审者解释这些占位符的含义。

所有被收录的文档都保留其中英文对照版本，因为项目自身的卫生规则如此要求：生成的说明
以成对形式写入；排除 `.github/` 是整目录排除，不会拆散任何一对文档。

## 脱敏内容与做法

动手之前先把审阅方给出的线索清单逐条与仓库比对，命中位置与处理方式如下：

| 位置 | 内容 | 处理 |
| --- | --- | --- |
| `.github/CODEOWNERS:3`、`.github/ISSUE_TEMPLATE/config.yml:4` | 所有者账号、仓库 URL | 整个 `.github/` 排除 |
| `LICENSE:3` | `Copyright (c) 2026, Zichao Zeng and Yinzhao Wang` | 仅以占位符替换两位持有人姓名，许可全文保留（见下文） |
| `pyproject.toml:12,56-60` | 作者姓名与邮箱、五条仓库 URL | 占位符替换 |
| `CITATION.cff:14-27,31,35` | 两位作者的姓名、邮箱、ORCID、机构、Zenodo DOI、仓库 URL | 占位符替换 |
| `.zenodo.json:11-19,31,35` | 同上，位于 `creators`/`related_identifiers` | 占位符替换 |
| `README.md:3-6,12,128,143-145` / `README_CN.md:3-6,12,117,132-134` | CI 与 Zenodo 徽章、仓库 URL、DOI、含机构与 ORCID 的作者栏 | 占位符替换 |
| `SECURITY.md:19` / `SECURITY_CN.md:18` | 安全报告邮箱 | 占位符替换 |
| `CODE_OF_CONDUCT.md:65,137` / `CODE_OF_CONDUCT_CN.md:50` | 举报邮箱 | 占位符替换 |
| `CONTRIBUTING.md:13` / `CONTRIBUTING_CN.md:14`、`RELEASE.md:85` / `RELEASE_CN.md:80` | clone/install URL、Zenodo GitHub 设置页 URL | 所有者段替换 |
| `tests/test_release_hygiene.py:191,194,202` | 卫生测试本身断言了真实姓名、邮箱与 URL | 用同一张映射表一并替换（理由见下） |
| `docs/USER_GUIDE.md`、`docs/USER_GUIDE_CN.md`、`validation/*.py`、`examples/*.py`、`benchmarks/*`、`testing_suite/*` | 经检查：**未发现**作者字符串、仓库 URL 或本机绝对路径 | 原样收录 |

未发现任何家目录绝对路径（macOS 或 Windows 的用户目录前缀，也没有本检出自身的路径），
未发现资助号 `42422209`、`42272354`、`2023YFC3108600`：这些信息在稿件里，不在仓库里。
它们仍留在守护检查的拒绝清单中，一旦将来被提交，构建就会失败而不是把它们漏出去。
（本文档刻意不书写任何绝对路径字面量，正因为 `test_release_hygiene.py` 也会扫描它：
它是源树内一个受扫描的 `.md` 文件，写下被禁模式反而会弄巧成拙地弄挂测试。）

脱敏采用**一次性、统一的字面量到占位符替换，作用于所有被收录的文本文件**——文档、
代码、元数据*以及测试*。这是整个设计的关键取舍：`tests/test_release_hygiene.py` 要求
`pyproject.toml`、`CITATION.cff`、`.zenodo.json`、`README*.md` 与 `LICENSE` 指向同一
作者与同一仓库，因此只替换其中一部分会让评审者拿到一个“自身一致性测试失败”的包。
全量替换则让测试读取的字面量与元数据携带的字面量始终相同，套件在存档内原样通过。
所用占位符：`Author One Withheld`、`AuthorTwo`、`Withheld`、
`AUTHOR_WITHHELD_FOR_PEER_REVIEW`、`author.one.withheld@example.invalid` /
`author.two.withheld@example.invalid`（RFC 2606 保留域 `.invalid`，不会被误投递）、
`WITHHELD_FOR_PEER_REVIEW`（ORCID 字段）、`Institution Withheld For Peer Review`、
`School Withheld For Peer Review`、`10.5281/zenodo.WITHHELD_FOR_PEER_REVIEW`、
`GRANT_WITHHELD_FOR_PEER_REVIEW`。另有两条规则是为替换自身的盲区而设：其一是家目录与
本检出的绝对路径（构建时由 `Path.home()` 与账户名推导，不写死）；其二是 `README.md:143`
把 `Shanghai⏎Jiao Tong University` 拆在了两行，整串机构名规则无法命中——若没有
`Jiao Tong University` 与 `Shanghai` 两条规则，机构信息就会原样留在包里。

源文件从未被修改：脱敏后的内容只写进压缩包，`test_source_files_were_not_modified`
钉住了工作树仍然书写真实作者。

**ORCID 字段用占位符，绝不用伪造的编号。** 早先的构建把两位作者的 iD 换成
`0000-0000-0000-0001` 与 `-0002`：它们长得像已注册的编号，却解析不到任何记录，元数据
采集方会把它读作存档缺陷而不是匿名化。（这两个串的校验和其实合法：全零主体恰好使 ISO 7064 校验位为 `1`，
这正说明单靠校验和不足以充当安全测试。）
现在的脱敏在 `CITATION.cff` 与 `.zenodo.json` 中写入 `WITHHELD_FOR_PEER_REVIEW`，保留键名
以免元数据不可解析，并由守护检查以两种方式执法：合成的全零形式是 scan-only 字面量；
`orcid_violations()` 对包内**每一个**形如 ORCID 的串按 ISO 7064 Mod 11-2 重算校验位，
不合格即拒绝构建。因此第三方的真实 iD 照旧通过，植入的假 iD 通不过。

### 刻意保留的姓名

`NOTICE.md`、`NOTICE_CN.md` 以及 `CITATION.cff` 的 `references:` 条目记录了 treesliceR
的作者（Araujo、Ferreira、Nakamura、Coelho、Rangel）及其 ORCID。他们是第三方，不是本
投稿的作者；删除他们等于删掉必需的署名致谢，也不会带来任何匿名效果。因此守护检查是按
字面量拒绝*本投稿两位作者*的 ORCID，而不是拒绝任意 `orcid.org/…` 字符串，也不把
`Araujo` 视为身份信息。来源清单中的 `third_party_names_kept` 记录了这一点。

## `LICENSE`：唯一无法两全的要求

BSD-3-Clause 是*由具名持有人授予*的许可；其第 1 项条件要求保留版权声明，而项目自身的
卫生测试断言 `assert "Zichao Zeng" in license_text, "the licence must name the copyright
holder"`。于是评审包要么含作者姓名，要么就没有许可——不存在既诚实又合法的第 3 种做法，
而“剥掉许可的压缩包”是不可接受的：那会歪曲代码所提供的授权条款。

**决定已经作出并写入代码：送审存档就是 BSD-3-Clause 许可的 `placeholder` 构建——许可
全文照发，仅把两位持有人姓名替换为占位符**——
`Copyright (c) 2026, Author One Withheld and Author Two Withheld`。三项条件与免责声明逐字未改，评审者拿到的
授权范围完全不变；真实持有人在公开发布版本中恢复，而那才是法律上必须载名的副本。
替代方案是已实现的、不只是文字描述：

```bash
python tools/make_review_archive.py --license preserved
```

`preserved` 模式原样收录 `LICENSE`（含真实姓名），在来源清单的 `declared_exceptions`
中记录原因，守护检查随后只容忍这一个文件中的这一处命中，其他任何命中仍然导致失败。
若你更希望许可文本一字不动、并依赖“该包仅私密递交编辑部”这一事实，请选 `preserved`；
若你希望送审副本中一个名字都不出现，请用默认模式。注意 `preserved` 会让包内出现一行
真实身份，因此 Cover letter 中必须如实说明；并且它会牺牲存档自身的套件：实测
`preserved` 构建的包内 `tests/test_release_hygiene.py` 为 10 通过 / 1 失败，失败项正是
`assert "Author One Withheld" in license_text`——被脱敏的测试遇上了逐字保留的许可。
这就是该取舍本身的如实陈述，而不是把检查调弱以把它藏起来。

这件事不留任何待作者确认的悬项：模式已在打包工具中钉死为
`RESOLVED_LICENSE_MODE = "placeholder"`，每次构建都把结果写入来源清单的
`anonymisation.licence_decision`（`resolved: true`、既定模式、本次构建所用模式、上面那段
BSD-3-Clause 陈述，以及在刻意偏离时附上的 `deviation` 说明），守护检查的收尾一行会打印
该既定值，而 `check()` 会拒绝任何来源清单未把许可模式判定为已解决、或包内文件仍带有
“未完成决定”标记（`AUTHOR ACTION`、`confirm which mode`、`TODO:`）的存档。

## 来源清单（provenance sidecar）

包内写入 `REVIEW_ARCHIVE_PROVENANCE.json`，记录：从
`src/phyloslicer/__init__.py` 读到的包版本；工作树提交、分支、`git describe`、发布标签
`v0.1.0` 及其提交，以及 HEAD 与发布准备提交 `b151683` 的关系（**HEAD 是发布之后的维护
提交，并非发布标签**——Zenodo 存档记录切自 `b151683`，因此把包与 DOI 对照的编辑会落到
另一个提交上，必须提前说明）；构建所用的 Python 版本与平台；统一使用的时间戳及其来源；
许可模式，**以及已解决的许可决定（本次构建所用模式与对既定模式的任何偏离）**；ORCID
策略；语言分布（`language_mix`：包内有多少文档含中文、其中多少是 `*_CN.md` 对应件、
中文字符总数与最长的 12 个文件——这就是说明中“英文用户指南 + 中文用户手册”那句话的
证据）；实际触发的脱敏规则、逐文件的脱敏记录与命中次数；每一条排除路径及其原因；
以及包内每个文件的 SHA-256 与字节数。

来源清单刻意不复述*被删除的*字面量——只记录规则名（例如“first author e-mail”）与其
替换结果。若把 `Zichao Zeng → Author One Withheld` 写进审计记录，等于在 documenting
匿名化的同时泄露匿名对象，所以这张映射表只存在于本文档与打包工具中，而两者都不进包。
`test_provenance_never_quotes_a_removed_literal` 负责把关。清单无法记录自身的摘要，
因此它对自身除外；`--check-only` 知道这一个例外，其余任何未登记、缺失或字节不符的条目
一律拒绝。

构建时读取工作树带来一个后果：存档只能声称工作树能证明的事，因此打包工具现在拒绝
过度声称。两种状态会直接停止构建：

* **读不到 git 元数据**（没有 `.git`，或 PATH 上没有 `git`）。此前它把字面量 `"unknown"`
  写进 `source_control.commit`，于是包内套件以 `assert 7 == 40` 失败——丢掉 `.git` 的评审者
  既遇到红灯，又遇到一个伪装成哈希的假话。现在构建会以 `PROVENANCE GATE` 信息被拒绝并
  说明缘由；加 `--allow-no-git` 才继续，并把 `commit: null` 与 `git_unavailable: true`
  写入清单，即完全不声明任何提交；`tests/test_review_archive.py` 读到该标记后跳过提交断言
  而不是失败。
* **工作树存在未提交的改动。** 在几十个路径与 HEAD 不一致时只写 HEAD，是假的来源声明，
  因此除非显式加 `--allow-dirty` 否则拒绝构建；加了之后，清单会记录 `git_dirty: true`、
  `git_dirty_count` 与完整的 `git_dirty_paths`，并在 `source_control.note` 里明确声明 HEAD
  不足以描述包内字节。`check()` 会拒绝“脏构建却未如实记录”，也会在 git 可读时拒绝任何不是
  40 位小写十六进制的提交。

诚实的工作流因此是：先由主编提交，再从干净检出构建（无需任何开关），并以 CI 中位于测试
之后的那一步生成最终递交文件。`--allow-dirty` 的存档之所以合法，正是因为它自我声明为
脏构建产物。

## 评审者应当能够运行什么

以下每条都在解包后的存档内、仅凭存档本身、无需安装步骤实测通过：

| 命令 | 结果（CPython 3.11 实测） |
| --- | --- |
| `python -m pytest -q` | 全套通过：测试台内约 12 s，首次冷解包约 27 s，含发布卫生测试 |
| `python examples/quickstart.py --sites 24 --species 60 --trees 3` | 约 0.6 s，正常输出拟合表 |
| `python validation/compare.py` | 约 2.5 s，`verdict: PASS`；无需 `PYTHONPATH`，无需安装 |
| `python benchmarks/recompute_ratios.py` | 重新推导出文中引用的全部加速比区间 |

上面每一条都能照原样跑通：不需要安装，也不需要设 `PYTHONPATH`。
`examples/quickstart.py` 与基准脚本自行把 `src/` 加入 `sys.path`，pytest 由
`pyproject.toml` 的 `pythonpath = ["src"]` 获得，而 `validation/compare.py`
现在也在 import 处做了同样的事（在它的 `phyloslicer` 导入之前
`sys.path.insert(0, <repo>/src)`），所以 README 里那句不加环境变量的
`python validation/compare.py` 字面上就是成立的。2026 年 10 月 8 日在未设该
变量的环境下实测：`python validation/compare.py` 退出码 0，打印
`verdict: PASS`。本文件早先的版本曾声称比较器是唯一仍需要 `PYTHONPATH=src`
或安装的脚本，那是照旧快照写的，此说法撤回。

其中 `validation/compare.py` 对稿件的“哈希登记参考值”主张最关键：它读取
`validation/reference/manifest.txt`，清单条目只记录裸文件名（比较器相对 `reference/`
目录解析它们，不依赖任何被记录的路径），在比对前先校验每个 SHA-256，并拒绝未登记的
参考文件。这条链上没有任何一处引用作者的开发机，因此它能穿过匿名化、在存档内重新跑
通——论文里“评审者可自行核验哈希”那句话所需要的正是这一点。

存档内执行类测试是 `tests/test_review_archive.py` 中最慢的部分：实测该模块 19 个测试
共 26 s，其中约 15 s 来自解包后的套件、示例与比较器。这远在 3 分钟预算之内，因此
默认不做跳过；若某台机器连这段耗时也不希望承担，可设
`PHYLOSLICER_SKIP_ARCHIVE_RUN=1`，三条执行测试会跳过，而构建、脱敏、来源清单、
可复现性与守护检查测试照常运行。CI 中**不**跳过。

## 致编辑

作者身份**仅**在评审期间、且仅在压缩包内被隐去。稿件标题页与仓库中的记录版本仍然载有
身份：`CITATION.cff`（姓名、邮箱、ORCID、机构）、`LICENSE`
（`Copyright (c) 2026, Zichao Zeng and Yinzhao Wang`）、`.zenodo.json`，以及公开版本中 `v0.1.0` 标签提交
上的 README 作者栏。匿名化包里没有任何一处是永久改动：打包工具只把脱敏后的字节写进
压缩包，从不编辑源文件，因此 Zenodo 与 GitHub 在录用后承载的存档仍然写明作者；包内元
数据字段读作 `Withheld`，而不是“尚未完成”的占位符。来源清单记录了构建所用的确切提交，
便于两相核对；包内的 `PEER_REVIEW_NOTICE.md` / `PEER_REVIEW_NOTICE_CN.md` 会向评审者
解释这些占位符的含义，但不会解释它们是谁。
