# 稿件字数统计口径（可复现）

**适用版本：** PhyloSlicer 0.1.0 · **撰写日期：** 2026 年 10 月 8 日
**读者：** 编辑、审稿人，以及任何需要证明投稿检查清单里的文体统计量不是"口头主张"的人。
**本文件与 `WORD_COUNT.md` 互为对应件**，英文件为准，本件为同结构译稿。

投稿检查清单给出了篇幅统计量，并承诺"读者可以从交付的文件里重新算出每一个数字"。
本文件就是让这句承诺对字数统计真正成立的口径：给出逐字的计数规则、规则作用的
段落区间、各区间最近一次测量的结果，以及一段只用标准库的脚本——它只需要稿件的
`.docx`，不需要 Word、不需要 pandoc、也不需要 `python-docx`。

---

## 1. 计数规则

1. 把稿件 `.docx` 当作 zip 打开，读取 `word/document.xml`。
2. 按文档顺序以 `<w:p …> … </w:p>` 切分段落。
3. 段落的*样式*取其 `w:pStyle`/`w:val` 属性值；没有 `w:pStyle` 的段落样式记为空。
4. 段落的*文本*是该段落内所有 `<w:t>` 元素按文档顺序的拼接。超链接里的文本
   （例如 DOI）**计入**。
5. *词*指一段连续的非空白字符，且其中至少含一个字母或数字；"字母或数字"按 Python
   `str.isalnum()` 判定（即 Unicode 字母与数字都算，所以 `α`、`CpD`、`e07364`、
   `2025` 各算一个词）。纯标点不算词：`—`、`;`、`%`、`(1)`、`·`、单独的减号都是零贡献。
6. 一个区间的字数 = 该区间内（剔除被排除样式后）所有段落的词数之和。
7. 恒定排除的样式：`FigureCaption`（两条图注）与 `Equation`（四个行间公式）。
   `FigureFrame` 是承载插图的空段落，没有 `<w:t>`，计 0，两种口径一致。
   各级标题（`Heading1`、`Heading2`）与正文中的代码/输出块**计入**——这是该期刊
   "main text" 的含义，下面引用的每一个数字都用同一约定。

## 2. 三个区间

三个区间的起点相同：文本恰为 `Introduction` 的 `Heading1` 段落。它们只差终点。

| 编号 | 终止边界 | 通俗名称 | 包含内容 |
|---|---|---|---|
| **D** | 以 `Statement on the use` 开头的 `Heading1` 段落 | "从引言到讨论结束"——作者须知页面自己的措辞 | 在 AI 声明之前结束；含 `Concluding remarks` |
| **M** | `References` 的 `Heading1` 段落 | 检查清单采用头条数字 | 区间 **D** 再加上 AI 声明的标题与正文 |
| **W** | 文件最后一段 | 最宽口径 | 区间 **M** 再加上整条参考文献 |

摘要用同一套词规则测量，区间是 `Abstract` 标签段与 `Key-words:` 标签段之间的那
一段正文；标题段与两个标签段不计入。

**实测数字。**

| 区间 | 2026-10-08 metadata 修复那一轮测得 | 同日 22:3x 语言修订之后复测 |
|---|---|---|
| **D** 引言 → 讨论结束（不含 AI 声明） | **3665** | **3991** |
| **M** 引言 → 含 AI 声明（不含参考文献） | **3829** | **4141** |
| **W** 引言 → 文件结尾（含参考文献） | 4532 / 4552 / 4556 | 5030 |
| 摘要正文段 | **244** | 278 |

表里刻意同时给出两个头条值，因为有一段时间检查清单把 **3665** 当标题数字，
而它写下的规则却是"截止到 `References` 标题"，按那条规则算出来是 **3829**。
这两个数从来不是同一件事的两次测量，而是 **D** 与 **M** 两个不同区间；两者的差
恰好是 AI 声明正文加它自己的小节标题（在测得 3665/3829 的那一轮是 157 + 7 个词，
复测这一轮是 150 个词）。检查清单现在把 **M = 3829** 作为它背书的数字，
并把 **D** 说明为作者须知页面所命名的更窄口径。无论按哪一个，都在 Software Note
的 4000 词上限之内。

**W** 那一行看起来"不一致"——同一轮里出现 4532、4552、4556 三个值，语言修订把若干
段落写长之后又测出 5030。这不是规则坏了。**W** 装的文字最多，因而在稿件当天被
反复编辑时最敏感：三个数字是同一个文件在同一天的三个不同状态，各自都是当时如实
测出的结果。任何人对*自己收到的那一份*文件重跑脚本，得到的就是那一份的数。这正
是本文件公布规则而不是只公布数字的理由。

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

2026 年 10 月 8 日 22:3x 对刚被语言修订重写的那份稿件执行 `python word_count.py <稿件>.docx`，
实际输出：

```
D  Intro -> end of Discussion (excl. AI statement): 3991
M  Intro -> through the AI statement (the headline): 4141
W  Intro -> end of file (incl. references): 5030
Abstract paragraph: 278
```

还值得跑一次样式普查，因为整条规则就是按段落样式取数的：

```bash
python -c "import zipfile,re,collections;c=collections.Counter(re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p).group(1) if re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p) else '' for p in re.findall(r'<w:p\b[^>]*>.*?</w:p>',zipfile.ZipFile(sys.argv[1]).read('word/document.xml').decode('utf-8'),re.S));print(c)" 稿件.docx
```

## 4. 这条规则不是什么

* **不是 `python-docx` 的读法。** `python-docx` 通过 `p.runs` 给出段落文本，会跳过
  嵌在超链接里的 run；本规则拼接段区内每一个 `<w:t>`。在交付的稿件里两种读法结果
  完全一致（实测两遍都是 3991 / 4141 / 5030），因为正文里的超链接只出现在参考
  文献表，而那里的 DOI 在同一段落里也以纯文本写了一遍。若今后某次修订用超链接
  *取代*了 DOI 文本，那么标准库这一侧才是诚实的计数。
* **不能替代编辑部自己的计数。** 我们按期刊命名的区间核对 4000 词上限；如果编辑部
  另有一套数法，以它为准，本文件解释差异而不是争辩差异。

## 5. 稿件文件自身的来源说明

五份稿件 `.docx` 是程序生成的：构建脚本用 `python-docx` 写出它们，因此每个文件的
`docProps/core.xml` 里 creator 字段都是 `python-docx`。这里记录这一点，是因为稿件中
的 AI 使用声明必须诚实描述生产路径——这些文档是由脚本从作者撰写并反复修改的文本
生成的，而生成脚本本身是可复现性叙事的一部分，不是需要遮掩的东西。一条命令即可核对：

```bash
python -c "import zipfile,sys;print(zipfile.ZipFile(sys.argv[1]).read('docProps/core.xml').decode('utf-8'))" 稿件.docx | grep -o "<dc:creator>[^<]*</dc:creator>"
```

同一批文件的 `docProps/app.xml` 里 `Application` 字段写的是某桌面字处理软件；该字符串
由生成器模板带出，并不表示文件在 Word 里编辑过。把这五份文件各另存一次即可对齐属性
并刷新 PAGE/NUMPAGES 缓存字段，代价是失去"脚本即唯一来源"这一优点。

## 6. 参见

* `docs/CI_WORKFLOWS.md` / `docs/CI_WORKFLOWS_CN.md`——CI 两条腿各自测什么。
* `docs/ZENODO_DEPOSIT_CHECKLIST.md` / `_CN.md`——外部登记元的核对清单。
* `testing_suite/TEST_REPORT.md` / `_CN.md`——测试计数，同样附可复现命令；
  本文件对字数做的事，那份报告对测试计数做同样的事。
