# Zenodo 存档核对清单 —— 只有作者能执行的外部步骤

**适用版本：** PhyloSlicer 0.1.0 · **撰写日期：** 2026 年 10 月 8 日
**读者：** 两位作者。本文件的每一步都无法在代码仓库内部完成，也不应当在仓库内部
完成：它记录的是软件包、稿件与投稿信*对外部登记元数据*的主张，以及如何让公开记录
与主张一致。
**本文件与 `ZENODO_DEPOSIT_CHECKLIST.md` 互为对应件**，英文件为准，本件为同结构译稿。

**为什么需要它。** 投稿材料里有四份文件对"存档元数据"作了事实陈述——说它登记了两位
作者、说它记录的是 9 月 30 日发布。线上记录并不是这样。另有一个字段被本文件的第一轮
读反了方向，而这个误读随后被三份稿件文件抄走：许可证在线上记录上**是**登记了的
（DataCite kernel-4 的 `rightsList`，SPDX `bsd-3-clause`），所以"存档不登记许可证"是
错的、不该再印——第 1 节解释为什么测量会显示 `None`。本地的 `CITATION.cff` 与
`.zenodo.json` 是对的；但登记库不是文件，改树里的文本无法弥合剩余的缺口。本文件把这
个缺口逐字段写清，并给出测量它的命令。

---

## 1. 操作前先验证，操作后再验证（一条命令，无需凭据）

```bash
curl -s https://api.datacite.org/dois/10.5281/zenodo.23058460 > record.json
python - <<'PY'
import json
record = json.load(open("record.json"))["data"]["attributes"]
for key in ("titles", "creators", "contributors", "dates", "publicationYear",
            "rights", "rightsList", "fundingReferences", "relatedIdentifiers",
            "version"):
    print(f"{key}: {record.get(key)}")
PY
```

请读 `rightsList`，不要读 `rights`。Zenodo 写出的 kernel-4 payload 里根本没有 `rights`
这个键，所以即使许可证登记得完全正确，`record.get("rights")` 也会打印 `None`——本文件
第一轮、以及抄用它的三份稿件文件，就是这样把记录说成"没有许可证"的。许可证是登记了的。

上面的 version DOI 是稿件所引用的可cite对象；concept DOI 解析到该系列的最新版本，
因此它**不**钉住 0.1.0。请对 concept DOI 跑同一条命令——它当前的元数据与 0.1.0 记录
一致，这正是"引用 concept DOI 就继承了下面每一条缺陷"的原因。

**本轮（2026-10-08）亲自复核到的线上取值**（DataCite kernel-4 JSON，`state: findable`，
`created` 与 `registered` 均为 `2026-09-30T12:13:27Z`，`publisher: Zenodo`，`language: en`，
`types.resourceTypeGeneral: Software`，`types.ris: COMP`）：

| JSON 路径 | 当前值 |
|---|---|
| `data.attributes.titles` | `[{"title": "PhyloSlicer"}]` |
| `data.attributes.creators` | 只有 1 条——通讯作者，含 ORCID iD 与完整单位 |
| `data.attributes.contributors` | 1 条——第二作者被登记为 `contributorType: "Supervisor"`（导师），并带其 ORCID iD |
| `data.attributes.dates` | `[{"date": "2026-09-29", "dateType": "Issued"}]` |
| `data.attributes.publicationYear` | `2026` |
| `data.attributes.rightsList` | 1 条——`BSD 3-Clause "New" or "Revised" License`，`rightsIdentifier: bsd-3-clause`，SPDX 词表。许可证确实挂在可cite对象上（`data.attributes.rights` 在 payload 里是"不存在"，不是"空值"） |
| `data.attributes.fundingReferences` | `[]` |
| `data.attributes.version` | `0.1.0` |
| `data.attributes.relatedIdentifiers` | `IsSupplementedBy` 仓库地址、`IsDerivedFrom` doi:10.1111/ecog.07364、`IsVersionOf` concept DOI |

## 2. 逐字段：线上当前值 → 目标值

| 字段 | 线上当前值 | 目标值 | 如何验证 |
|---|---|---|---|
| `titles[0].title` | `PhyloSlicer` | `CITATION.cff` 里的完整标题串，即 `PhyloSlicer: temporal slicing of phylogenies and cumulative-diversity rate analysis for macroecology in Python` | `titles` |
| `creators` | 1 位 creator | 2 位 creator，按论文作者顺序，每人带 ORCID iD **与**单位；逐字照抄 `CITATION.cff` 的 `authors:` 块（或 `.zenodo.json` 的 `creators` 数组），两处本地文件已经写全 | `creators` 长度与各自 `nameIdentifiers` |
| `contributors` | 第二作者被放在这里，类型 `Supervisor` | 删除该条：合著者是 creator，不是导师。`contributors` 只留给真正的非作者（本作品没有） | `contributors` |
| `dates[dateType=Issued]` | `2026-09-29` | 稿件、`CITATION.cff`（`date-released`）、`.zenodo.json`（`publication_date`）三处一致的发布日：`2026-09-30`。DataCite 的 `Issued` 由 Zenodo 的 *Publication date* 字段导出，所以改那个字段，而不是改 DOI | `dates` |
| `publicationYear` | `2026` | 不变 | `publicationYear` |
| `rightsList` | 上面那条 BSD-3-Clause，已经登记 | 保持原样——无需新增。Software Notes 要求"许可证在软件与说明中都要写清"，可cite对象上就靠这一字段兑现，而它已经达标；不要重开这个字段，也不要再印"存档不登记许可证" | `rightsList` |
| `fundingReferences` | `[]` | 若作者决定登记，则登记 3 项：两项国家自然科学基金、一项国家重点研发计划。只使用标题页上出现的编号——不要凭记忆转录、不要虚构；仓库内部的匿名化对照表（`tools/make_review_archive.py`，按设计不进审稿包）保存着同样三条串，可作内部交叉核对 | `fundingReferences` |
| `version` | `0.1.0` | 保持 `0.1.0`，除非编辑部的版本裁定要求重新打 tag——见第 3 节 | `version` |
| `description` | 存档里的摘要 | 用 `.zenodo.json` 的 `description` 重新生成，使记录与该文件逐字一致 | `descriptions` |
| `relatedIdentifiers` | 见上表 | 保持不变；若发新版本，则两条 version DOI 之间由 Zenodo 在"new version"流程里自动写入 `IsPreviousVersionOf`/`IsNewVersionOf` | `relatedIdentifiers` |

## 3. 关于 0.1.0 的版本裁定会怎样影响 DOI

投稿信请编辑部裁定：`0.1.0` 是否满足 Software Notes 对"稳定可cite版本"的期待。元数据
路线取决于这条裁定，而两种结果并不对称：

* **若 `0.1.0` 就是版本记录。** 已被引用的 Zenodo 版本不可变：它的 version DOI 与其
  DataCite 元数据无法就地更正。于是只有两条路：(i) 那四份断言"存档元数据登记了两位
  作者、9 月 30 日发布"的文件在印刷上就是错的，可用性一节必须补上
  披露句——线上记录使用短标题、单一 creator、9 月 29 日的 Issued 日期，而它的许可证
  （BSD-3-Clause）确实登记了；作者名单以 `CITATION.cff` 为准；或 (ii) 作者既写披露句，又另外发布一份"代码
  相同、元数据已修正"的新版本：稿件仍引 `0.1.0`，但面向后来者的可cite对象被修好了。
  只做 (i) 是诚实的下限；(ii) 是更好的记录。
* **若编辑部裁定 `0.1.0` 不够。** 那就需要重新打 tag（`1.0.0`），这必然在**同一个
  concept DOI 之下铸造一个新的 version DOI**。此时本地文件要与存档按如下顺序一起改：
  `.zenodo.json`（`version`、`publication_date`）→ 打 tag 与发布 → `CITATION.cff`
  （`version`、`date-released`，以及 `identifiers:` 中新 version DOI 那一条；concept DOI
  那条保留）→ `CHANGELOG.md` / `CHANGELOG_CN.md` → `README.md` / `README_CN.md` 的版本横幅
  → 最后才改稿件里的版本句与日期句及其软件自引条目。不要删除或撤下 `0.1.0` 记录：
  已有引用要靠它解析，新版本的 `IsPreviousVersionOf` 链接才是让谱系保持可读的东西。

无论哪一种，**都不要**把 concept DOI 当作钉死的引用：它永远解析到系列的最新版本。
需要钉住"论文随附的那份代码"的说明，必须引用 version DOI；`CITATION.cff` 自己在
`identifiers:` 的 description 里就写明了这一点。

## 4. Zenodo 上的操作顺序

1. 重跑第 1 节命令并保存输出——这就是"操作前"的证据。
2. 打开 `0.1.0` 记录，使用 **New version**（保留谱系并铸造新的 version DOI；只有在
   准备一份尚未发布的草稿时，"就地改草稿"这条路才存在）。确认新草稿继承了文件，
   然后按第 2 节逐字段修正。许可证只能通过 Zenodo 自带的许可证选择器设置（它写进
   `rightsList`）；把许可证名字写进描述文本并不构成登记权利。`0.1.0` 记录已经由选择器
   带上了 `bsd-3-clause`，所以"new version"只需要把它留住。
3. 用期刊 Funder Registry 的条目填写 **Funding**，并逐字挂上标题页上的项目编号。登记库里
   写错一个编号，比留空更糟。
4. 发布，然后对**新的** version DOI 与 concept DOI 重跑第 1 节，确认：2 位 creator、
   没有 `Supervisor` 条目、完整标题、`rightsList` 里仍是 BSD-3-Clause、`Issued` = 发布日、3 条
   `fundingReferences`。
5. 只有第 4 步与仓库一致之后，那四份文件关于"存档元数据"的句子才为真。若作者选择
   保持 `0.1.0` 记录原样，请把第 3 节的披露句在同一次修订里写进可用性一节与投稿信。

## 5. 本文件刻意不做的事

它不改动登记库，也不在这里复述作者姓名、ORCID iD 或项目编号：审稿包会收录本文件，
而包的匿名守卫会拒绝任何审稿可见文件里残留的作者身份字面量。这些字面量存在于
`CITATION.cff`、`.zenodo.json`、`LICENSE` 与 `NOTICE.md`，打包器在构建审稿包前会把它们
替换掉——这也正是"本地文件必须成为唯一事实来源"、而本清单只指向它们的原因。

## 6. 参见

* `CITATION.cff` / `.zenodo.json`——已修正的本地元数据。
* `CHANGELOG.md`——版本记录与其发布日期。
* `docs/CI_WORKFLOWS.md` / `docs/CI_WORKFLOWS_CN.md`——哪些由 CI 自动核验，哪些永远不可能。
* `docs/WORD_COUNT.md` / `docs/WORD_COUNT_CN.md`——同一套纪律的另一份说明：公布口径，而不只是数字。
