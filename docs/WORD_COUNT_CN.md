# 稿件字数统计口径（可复现）

**适用版本：** PhyloSlicer 0.1.0 · **撰写日期：** 2026 年 10 月 8 日
**读者：** 编辑、审稿人，以及任何需要证明投稿检查清单里的文体统计量不是"口头主张"的人。
**本文件是 [`WORD_COUNT.md`](WORD_COUNT.md) 的中文译文**，正文以英文文件为准，两者保持同一章节结构。

投稿检查清单给出了篇幅统计量，并承诺"读者可以从交付的文件里重新算出每一个数字"。
本文件就是让这句承诺对字数统计真正有意义的口径：给出逐字的计数规则、规则作用的
段落区间、各区间最近一次运行时所针对的那份稿件快照测得的数字，以及一段只用标准库
的脚本——对它指向的*任何一份*稿件 `.docx`，它都能重新算出那一份的字数；不需要 Word、
不需要 pandoc、也不需要 `python-docx`。引用本文件里任何一个数字之前，请先读 §2
最后一段：可迁移的是规则，不是数字。

---

## 1. 计数规则

1. 把稿件 `.docx` 当作 zip 打开，读取 `word/document.xml`。
2. 按文档顺序以 `<w:p …> … </w:p>` 切分段落。
3. 段落的*样式*取其 `w:pStyle`/`w:val` 属性值；没有 `w:pStyle` 的段落样式记为空。
4. 段落的*文本*是该段落内所有 `<w:t>` 元素按文档顺序的拼接。超链接里的文本
   （例如 DOI）计入。
5. *词*指一段尽可能长的连续非空白字符，且其中至少含一个字母或数字；"字母或数字"按 Python
   `str.isalnum()` 判定（即 Unicode 字母与数字都算，所以 `α`、`CpD`、`e07364`、
   `2025` 各算一个词）。纯标点不算词：`—`、`;`、`%`、`(1)`、`·`、单独的减号都是零贡献。
6. 一个区间的字数 = 该区间内（剔除被排除样式后）所有段落的词数之和。
7. 恒定排除的样式：`FigureCaption`（两条图注）与 `Equation`（四个行间公式）。
   `FigureFrame` 是承载插图的空段落，没有 `<w:t>`，因此无论计入还是剔除都是 0。
   章节标题（`Heading1`、`Heading2`）与正文中的代码/输出块**计入**——这是该期刊
   "main text" 的含义，下面引用的每一个数字都用同一约定。

## 2. 三个区间

三个区间的起点相同：文本恰为 `Introduction` 的 `Heading1` 段落。它们只差终点。

| 编号 | 终止边界 | 通俗名称 | 包含内容 |
|---|---|---|---|
| **D** | 以 `Statement on the use` 开头的 `Heading1` 段落 | "从引言到讨论结束"——作者须知页面自己的措辞 | 在 AI 声明之前结束；含 `Concluding remarks` |
| **M** | `References` 的 `Heading1` 段落 | 检查清单所采用的头条数字 | 区间 **D** 再加上 AI 声明的标题与正文 |
| **W** | 文件最后一段 | 最宽口径 | 区间 **M** 再加上整份参考文献表 |

摘要用同一套词规则测量，区间是 `Abstract` 标签段与 `Key-words:` 标签段之间的那
一段正文；标题段与两个标签段不计入。

**实测数字。**

| 区间 | 2026-10-08 测得（修复稿件元数据那一轮） | 2026-10-08 22:2x 在同一台机器上复测 |
|---|---|---|
| **D** 引言 → 讨论结束（不含 AI 声明） | **3665** | **3991** |
| **M** 引言 → 含 AI 声明（不含参考文献） | **3829** | **4141** |
| **W** 引言 → 文件结尾（含参考文献） | 4532 / 4552 / 4556 | 5030 |
| 摘要正文段 | **244** | 278 |

表里刻意同时给出两个头条值，因为有一轮检查清单把 **3665** 印成了头条数字，
而它写下的规则却是"截止到 `References` 标题"，按那条规则算出来是 **3829**。
这两个数从来不是同一件事的两次测量，而是 **D** 与 **M** 两个不同区间；两者的差
恰好是 AI 声明正文加它自己的小节标题（在测得 3665/3829 的那一轮是 157 + 7 个词，
复测这一轮是 150 个词）。在那一轮里，检查清单把 **M = 3829** 作为它背书的数字，
并把 **D** 说明为作者须知页面所命名的更窄口径。无论按哪一个，在期刊所测量的那个区间
上，文章都在 Software Note 的 4000 词上限之内。

**W** 那一行看起来"不一致"——同一轮里出现 4532、4552、4556 三个值，语言修订把若干
段落写长之后又测出 5030。这不是规则坏了。**W** 这个区间包含的文字最多，因而在稿件当天被
反复编辑时最敏感：三个数字是同一个文件在同一天的三个不同状态，各自都是当时如实
测出的结果。任何人对*自己收到的那一份*文件重跑脚本，得到的就是那一份的数。这正
是本文件公布规则而不是只公布数字的理由。

**这些数字各属于哪一份文件。** 上表里每一个数字都是稿件在某一列所写日期那一份*快照*
的测量结果，而它们之中的任何一个都无法从今天项目文件夹里那份稿件 `.docx` 重新算出。
对 `PhyloSlicer-手稿正文.docx`（SHA-256 以 `5c90250d` 开头、以 `9de92598` 结尾，
2026 年 11 月 13 日复测）运行 §3 的脚本，得到的是 **D 2714、M 2821、W 3694、摘要 197**
——上表中没有一个值。这不是同一篇文字的不同读法，而是不同的文件：该文件全部段落、
全部样式的词数合计只有 4288，因此对它而言 **W = 5030** 与 **M = 4141** 在算术上不可
能成立。交付的检查清单又给出第三套数字：它把 **3969** 词作为自己背书的数字，
`References` 之前那一区间是 **3990**，更窄那一区间是 **3850**，摘要是 **276** 词，
而且它的规则在 `FigureFrame`、`FigureCaption`、`Equation` 之外还剔除 `Reference`。
三套数字，三份文件。
上文发表的数字一律照原样保留，并没有换成 2714/2821/3694/197，因为稿件字数是关于*所投
那篇文章*的编辑主张，而本仓库既不收稿件源文件，也不收组装这些文档的脚本，因此仓库里
不存在可以重新推出这三套数字中任何一套的凭据；把今天的数换进去，只会把一个有日期的
测量变成一个没日期的测量，并且会与真正投出去的那份检查清单悄悄对不上。本文件能担保的是
*规则*：把 §3 指向任何一份稿件 `.docx`，它都返回那一份的数——标题里"可复现"能表达的
也就是这个意思。三套快照中哪一套才是准的，属于作者的决定，不是本文件能够裁定的事。

## 3. 脚本

只用标准库（`zipfile`、`re`、`sys`）。存成 `word_count.py`（放在仓库之外的任意
位置），或直接粘进 `python -`，参数是稿件 `.docx` 的路径：

```python
"""Re-derive the manuscript word counts. Pure standard library."""
import re
import sys
import zipfile

EXCLUDED_STYLES = {"FigureCaption", "Equation"}

xml = zipfile.ZipFile(sys.argv[1]).read("word/document.xml").decode("utf-8")
paras = re.findall(r"<w:p\b[^>]*>.*?</w:p>", xml, re.S)


def style(par: str) -> str:
    found = re.search(r'<w:pStyle w:val="([^"]+)"', par)
    return found.group(1) if found else ""


def text(par: str) -> str:
    return "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", par, re.S))


def words(par: str) -> int:
    return len([t for t in text(par).split() if any(c.isalnum() for c in t)])


def heading(name: str) -> int:
    for i, par in enumerate(paras):
        if style(par) == "Heading1" and text(par).strip().startswith(name):
            return i
    raise SystemExit(f"no Heading1 paragraph starting {name!r}")


start = heading("Introduction")
ai = heading("Statement on the use")
refs = heading("References")
spans = {
    "D  Intro -> end of Discussion (excl. AI statement)": (start, ai),
    "M  Intro -> through the AI statement (the headline)": (start, refs),
    "W  Intro -> end of file (incl. references)": (start, len(paras)),
}
for label, (first, last) in spans.items():
    total = sum(words(p) for p in paras[first:last] if style(p) not in EXCLUDED_STYLES)
    print(f"{label}: {total}")

body = [p for p in paras[:start] if text(p).strip() not in {"", "Abstract"}
        and not text(p).strip().startswith("Key-words:")]
print("Abstract paragraph:", sum(words(p) for p in body if style(p) != "ManuscriptTitle"))
```

2026 年 10 月 8 日 22:3x 对刚被语言修订重写的那份稿件 `.docx`——即当时放在项目
文件夹中本仓库旁的那一份快照，*不是*今天放在那里的那一份——执行
`python word_count.py <manuscript>.docx`，实际输出：

```
D  Intro -> end of Discussion (excl. AI statement): 3991
M  Intro -> through the AI statement (the headline): 4141
W  Intro -> end of file (incl. references): 5030
Abstract paragraph: 278
```

同一段脚本今天原样跑在项目文件夹所交付的那份稿件 `.docx` 上（2026 年 11 月 13 日），
打印出来的是一篇不同的文章——差额由 §2 最后一段解释，而它跑的那份文件正是 §5 所
描述的那份：

```
D  Intro -> end of Discussion (excl. AI statement): 2714
M  Intro -> through the AI statement (the headline): 2821
W  Intro -> end of file (incl. references): 3694
Abstract paragraph: 197
```

还值得跑一次样式普查，因为整条规则就是按段落样式取数的：

```bash
# the paragraph styles the rule keys on, and how many of each the file holds
python -c "import zipfile,re,collections,sys;c=collections.Counter(re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p).group(1) if re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p) else '' for p in re.findall(r'<w:p\b[^>]*>.*?</w:p>',zipfile.ZipFile(sys.argv[1]).read('word/document.xml').decode('utf-8'),re.S));print(c)" MANUSCRIPT.docx
```

## 4. 两点说明：本规则不是什么

* **不是 `python-docx` 的读法。** `python-docx` 通过 `p.runs` 给出段落文本，会跳过
  嵌在超链接里的 run；本规则拼接段落内每一个 `<w:t>`。在上面引用的 10 月 8 日那一轮里
  两种读法结果完全一致（两遍都测得 3991 / 4141 / 5030），因为当时正文里的超链接只出现
  在参考文献表，而那里的 DOI 在同一段落里也以纯文本写了一遍。今天项目文件夹交付的这份
  稿件上两者依然完全一致，但成立的理由换成了另一个：该文件**一个超链接都没有**
  （`word/document.xml` 里 `<w:hyperlink>` 元素为 0 个，用 `python-docx` 1.2.0 遍历
  同样是 0 个），而它的 26 个 `Reference` 段落已经把 23 个标识符以纯文本写了出来
  （形如 `doi:10.…`；这几个字母另外唯一一次出现，是某条注释里的单词 `DOI`），因此两种
  读法没有可分歧之处，两边都返回 2714 / 2821 / 3694（2026 年 11 月 13 日复测）。若今后
  某次修订用超链接*取代*了 DOI 文本，那么标准库这一侧才是诚实的计数。
* **不能替代编辑部自己的计数。** Software Note 的 4000 词上限是按期刊所命名的那个
  区间来统计的；如果编辑部另有一套数法，以它为准，本文件解释差异而不是争辩差异。

## 5. 稿件文件自身的来源说明

五份稿件 `.docx` 是程序生成的：构建脚本用 `python-docx` 写出它们，而这一点在每个文件的
`docProps/core.xml` 里体现为两个字段，而不是一个。`dc:description` 带的是生成器自己的
印记 `generated by python-docx`。`dc:creator` 里没有任何人名，但它装的也不是
`python-docx` 这个字符串：那里写的是构建过程填进去的出处说明句 `PhyloSlicer submission
files, assembled programmatically (no individual author metadata)`。因此 `python-docx`
这个字面量出现在 description 里而不在 creator 里，只 grep `<dc:creator>` 的检查看到的
是那句出处说明而不是生成器名。这里记录这一点，是因为稿件中的 AI 使用声明必须诚实描述
生产路径——这些文档是由脚本从作者撰写并反复修改的文本生成的，而生成脚本本身是可复现性
叙事的一部分，不是需要遮掩的东西。一条命令即可核对这两个字段：

```bash
python -c "import zipfile,sys;print(zipfile.ZipFile(sys.argv[1]).read('docProps/core.xml').decode('utf-8'))" MANUSCRIPT.docx | grep -o "<dc:creator>[^<]*</dc:creator>\|<dc:description>[^<]*</dc:description>"
```

对五份文件它都会按上述顺序打印出同样的两行：

```
<dc:creator>PhyloSlicer submission files, assembled programmatically (no individual author metadata)</dc:creator>
<dc:description>generated by python-docx</dc:description>
```

同一批文件的 `docProps/app.xml` 里 `Application` 字段写的是某桌面字处理软件；该字符串
由生成器模板带出，并不表示文件在 Word 里编辑过。在 Word 或 LibreOffice 中把这五份文件
各另存一次，即可对齐属性并刷新 PAGE/NUMPAGES 缓存字段，代价是脚本构建不再是这些文件
的事实来源。

## 6. 参见

* `docs/CI_WORKFLOWS.md` / `docs/CI_WORKFLOWS_CN.md`——CI 各条流水线分别测什么。
* `docs/ZENODO_DEPOSIT_CHECKLIST.md` / `_CN.md`——外部元数据记录（Zenodo 记录）。
* `testing_suite/TEST_REPORT.md` / `_CN.md`——测试计数，并附上重新推出这些数字的命令；
  正是本文件对字数所采用的那一套纪律。
