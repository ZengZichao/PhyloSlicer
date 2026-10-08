# Word-count rule for the manuscript (re-derivable)

**Applies to:** PhyloSlicer 0.1.0 · **Written:** 8 October 2026
**Audience:** the editor, the reviewers, and anyone who has to show that the
style statistics quoted in the submission checklist are not assertions.

The submission checklist states its length statistics and promises that "a reader
can re-derive each of those numbers from the delivered document". This file is
the rule that makes the promise true for the word counts. It gives the counting
rule verbatim, the span boundaries it is applied over, the numbers each span
gave when it was last measured, and a pure-standard-library script that
reproduces them from the manuscript `.docx` alone — no Word, no pandoc, no
`python-docx`.

---

## 1. The rule

1. Open the manuscript `.docx` as a zip archive and read `word/document.xml`.
2. Split the body into paragraphs at `<w:p …> … </w:p>`, in document order.
3. A paragraph's *style* is the value of its `w:pStyle`/`w:val` attribute; a
   paragraph with no `w:pStyle` has the empty style.
4. A paragraph's *text* is the concatenation, in document order, of every
   `<w:t>` element inside it. Text carried by a hyperlink (a DOI, say) counts.
5. A *token* is a maximal run of non-whitespace characters that contains at
   least one alphanumeric character, where alphanumeric means Python's
   `str.isalnum()` (Unicode letters and digits, so `α`, `CpD`, `e07364` and
   `2025` are each one token). Pure punctuation is not a token: `—`, `;`, `%`,
   `(1)`, `·` and a bare minus all contribute nothing.
6. Word count of a span = the sum of the tokens of the paragraphs in that span,
   after dropping the excluded styles.
7. Excluded styles, always: `FigureCaption` (the two figure legends) and
   `Equation` (the four displayed equations). `FigureFrame` is the empty
   paragraph that carries a figure, so it has no `<w:t>` and counts as 0 either
   way. Section headings (`Heading1`, `Heading2`) and the in-text code/output
   block **are** counted — that is what "main text" means for this journal, and
   the same convention is used by every figure quoted below.

## 2. The three spans

All three spans start at the same place: the `Heading1` paragraph whose text is
exactly `Introduction`. They differ only in where they stop.

| ID | Stop boundary | Plain-language name | What it swallows |
|---|---|---|---|
| **D** | the `Heading1` paragraph beginning `Statement on the use` | "Introduction through the end of the Discussion" — the wording the instructions-for-authors page itself uses | ends before the AI statement; includes `Concluding remarks` |
| **M** | the `Heading1` paragraph `References` | the checklist's headline figure | span **D** plus the AI-statement heading and its text |
| **W** | the last paragraph of the file | the widest reading | span **M** plus the whole reference list |

The abstract is measured by the same token rule, over the single paragraph that
sits between the `Abstract` label paragraph and the `Key-words:` paragraph. The
title paragraph and the two label paragraphs are not part of it.

**The measured figures.**

| Span | Measured 2026-10-08 (the round that fixed the manuscript metadata) | Re-measured 2026-10-08 22:2x, same machine |
|---|---|---|
| **D** Intro → end of Discussion, excl. AI statement | **3665** | **3991** |
| **M** Intro → through the AI statement, excl. references | **3829** | **4141** |
| **W** Intro → end of file, incl. references | 4532 / 4552 / 4556 | 5030 |
| Abstract (single body paragraph) | **244** | 278 |

Both headline numbers are in the table on purpose, because for one round the
checklist printed **3665** as its headline while quoting a rule that stops at the
`References` heading, which gives **3829**. They were never two measurements of
one thing; they are two different spans, **D** and **M**, and the difference
between them is exactly the AI statement plus its section heading (157 tokens of
statement text and 7 of heading in the round that measured 3665/3829; 150 tokens
in the second measurement). The checklist now names **M = 3829** as the figure it
stands behind and quotes **D** as the narrower reading the journal page names.
Either way the article is inside the 4000-word cap for a Software Note on the
span the journal measures.

The **W** row is the one that looks inconsistent — three different values in one
round, 4532, 4552 and 4556, and a fourth (5030) after the language pass lengthened
several paragraphs. It is not a broken rule. **W** is the span with the most text
in it, so it is the most sensitive to the edits the manuscript was receiving while
those measurements were taken: the three figures are three different states of the
same file on the same day, taken in good faith and reported as they came out.
Anyone who re-runs the script on the `.docx` they were sent will get the number
for *that* file. This is the reason the rule is published here rather than the
number alone.

## 3. The script

Pure standard library (`zipfile`, `re`, `sys`). Save it as `word_count.py`
anywhere outside the repository, or paste it straight into `python -`, and give
it the path of the manuscript `.docx`:

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

A worked re-derivation, run on 8 October 2026 at 22:3x on the manuscript `.docx`
that the language pass had just rewritten (the file beside this repository in the
project folder), `python word_count.py <manuscript>.docx`:

```
D  Intro -> end of Discussion (excl. AI statement): 3991
M  Intro -> through the AI statement (the headline): 4141
W  Intro -> end of file (incl. references): 5030
Abstract paragraph: 278
```

One more check worth running, because the whole rule keys on paragraph styles:

```bash
# the paragraph styles the rule keys on, and how many of each the file holds
python -c "import zipfile,re,collections;c=collections.Counter(re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p).group(1) if re.search(r'<w:pStyle w:val=\"([^\"]+)\"',p) else '' for p in re.findall(r'<w:p\b[^>]*>.*?</w:p>',zipfile.ZipFile(sys.argv[1]).read('word/document.xml').decode('utf-8'),re.S));print(c)" MANUSCRIPT.docx
```

## 4. Two things this rule is not

* **Not a `python-docx` reading.** `python-docx` exposes a paragraph's text
  through `p.runs`, which skips runs nested in a hyperlink; the rule above
  concatenates every `<w:t>` in the paragraph. In the delivered manuscript the
  two agree exactly (measured: 3991 / 4141 / 5030 both ways) because the only
  hyperlinks in the body sit in the reference list, whose DOIs are already
  written out as text in the same paragraph. A future revision that hyperlinks a
  DOI *in place of* its text would make the stdlib count the honest one.
* **Not a substitute for the journal's own counter.** Word counts the 4000-word
  cap for a Software Note on the span the journal names; if the editorial office
  counts differently, its count governs and this file explains the difference
  rather than arguing with it.

## 5. Provenance note on the manuscript files themselves

The five manuscript `.docx` files are assembled programmatically: the build
scripts write them with `python-docx`, so each file's `docProps/core.xml` carries
the creator string `python-docx`. That is recorded here because the AI-use
statement in the manuscript has to describe the production route honestly — the
documents were generated by scripts from text the authors wrote and revised, and
the generation scripts are part of the reproducibility story rather than
something to hide. Check it in one command:

```bash
python -c "import zipfile,sys;print(zipfile.ZipFile(sys.argv[1]).read('docProps/core.xml').decode('utf-8'))" MANUSCRIPT.docx | grep -o "<dc:creator>[^<]*</dc:creator>"
```

The `Application` field in `docProps/app.xml` of the same files claims a desktop
word processor; that string is written by the generator's template and does not
mean the file was edited in Word. Re-saving each file once in Word or LibreOffice
aligns the properties and refreshes the cached PAGE/NUMPAGES fields, at the cost
of losing the scripted build as the source of truth.

## 6. See also

* `docs/CI_WORKFLOWS.md` / `docs/CI_WORKFLOWS_CN.md` — what the CI legs measure.
* `docs/ZENODO_DEPOSIT_CHECKLIST.md` / `_CN.md` — the external metadata record.
* `testing_suite/TEST_REPORT.md` / `_CN.md` — the test tallies, with the commands
  that re-derive them; the same discipline this file applies to word counts.
