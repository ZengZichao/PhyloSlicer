# Zenodo deposit checklist — the steps only the authors can run

**Applies to:** PhyloSlicer 0.1.0 · **Written:** 8 October 2026
**Audience:** the two authors. Nothing in this file can be done from inside the
repository, and nothing in it should be done from inside the repository: it is a
record of the *external* metadata that the package, the manuscript and the cover
letter assert, and of how to make the public record match.

**Why it exists.** Four documents in the submission state facts about "the deposit
metadata" — that it names both authors and that it records a 30 September release.
The live record does not say those things. One field was mis-read the other way in
the first pass of this file, and the mis-read then travelled into the three
submission documents: the licence *is* registered on the live record (DataCite
kernel-4 `rightsList`, SPDX `bsd-3-clause`), so "the deposit registers no licence" is
false and must not be printed — §1 explains why the measurement said `None`. The
local `CITATION.cff` and `.zenodo.json` are correct; a registry is not a file, so no
amount of editing in the tree closes the remaining gap. This file is the gap, field
by field, with the command that measures it.

---

## 1. Verify before and after (one command, no credentials)

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

Read `rightsList`, not `rights`. The kernel-4 payload Zenodo writes has no `rights`
key at all, so `record.get("rights")` prints `None` even when the licence is
registered correctly — which is how the first pass of this file, and the sentence in
the three submission documents that copied it, came to report a licence-less deposit.
The licence is registered.

The version DOI above is the citable object the manuscript names; the concept DOI
resolves to the newest version of the series and therefore does *not* pin 0.1.0.
Run the same command against the concept DOI too — its metadata currently mirrors
the 0.1.0 record, which is why a note that cites the concept DOI inherits every
defect below.

**Live values, re-verified by this pass on 2026-10-08** (DataCite kernel-4 JSON,
`state: findable`, `created` and `registered` both `2026-09-30T12:13:27Z`,
`publisher: Zenodo`, `language: en`, `types.resourceTypeGeneral: Software`,
`types.ris: COMP`, `schemaVersion: http://datacite.org/schema/kernel-4`):

| JSON path | Value as it stands |
|---|---|
| `data.attributes.titles` | `[{"title": "PhyloSlicer"}]` |
| `data.attributes.creators` | one entry — the corresponding author, with ORCID iD and full affiliation |
| `data.attributes.contributors` | one entry — the second author, `contributorType: "Supervisor"`, with their ORCID iD |
| `data.attributes.dates` | `[{"date": "2026-09-29", "dateType": "Issued"}]` |
| `data.attributes.publicationYear` | `2026` |
| `data.attributes.rightsList` | one entry — `BSD 3-Clause "New" or "Revised" License`, `rightsIdentifier: bsd-3-clause`, SPDX scheme. The licence is on the citable object (`data.attributes.rights` is absent from the payload, not empty) |
| `data.attributes.fundingReferences` | `[]` |
| `data.attributes.version` | `0.1.0` |
| `data.attributes.relatedIdentifiers` | `IsSupplementedBy` the repository URL, `IsDerivedFrom` doi:10.1111/ecog.07364, `IsVersionOf` the concept DOI |

## 2. Field-by-field: current live value → target value

| Field | Live value now | Target | Verify with |
|---|---|---|---|
| `titles[0].title` | `PhyloSlicer` | the full title string in `CITATION.cff`, i.e. `PhyloSlicer: temporal slicing of phylogenies and cumulative-diversity rate analysis for macroecology in Python` | `titles` |
| `creators` | 1 creator | 2 creators, in the paper's author order, each with ORCID iD **and** affiliation — copy them verbatim from the `authors:` block of `CITATION.cff` / the `creators` array of `.zenodo.json`, both of which already carry them | `creators` length and each `nameIdentifiers` |
| `contributors` | the second author sits here as `Supervisor` | delete this entry; a co-author is a creator, not a supervisor. Keep `contributors` for genuine non-authors only (there are none) | `contributors` |
| `dates[dateType=Issued]` | `2026-09-29` | the true release day the note, `CITATION.cff` (`date-released`) and `.zenodo.json` (`publication_date`) all state: `2026-09-30`. DataCite derives `Issued` from Zenodo's *Publication date* field, so change that field, not the DOI | `dates` |
| `publicationYear` | `2026` | unchanged | `publicationYear` |
| `rightsList` | the BSD-3-Clause entry above, already registered | keep as is — nothing to add. This is the field that answers the Software Notes rule "the licence must be clearly stated in the software *and* the note" on the citable object, and it is already met, so do not re-open it and do not print "the deposit registers no licence" | `rightsList` |
| `fundingReferences` | `[]` | three grants, if the authors choose to register them: two National Natural Science Foundation of China awards and one National Key Research and Development Program of China award. Use only the numbers on the title page — do not invent or transcribe from memory; the repository's own anonymisation table (`tools/make_review_archive.py`, deliberately excluded from the review archive) holds the same three strings and is the internal cross-check | `fundingReferences` |
| `version` | `0.1.0` | stays `0.1.0` unless the editor's version ruling forces a re-tag — see section 3 | `version` |
| `description` | the deposited abstract | regenerate from `.zenodo.json`'s `description` so the record and the file agree word for word | `descriptions` |
| `relatedIdentifiers` | as listed above | unchanged, plus — if a new version is issued — `IsPreviousVersionOf`/`IsNewVersionOf` links between the two version DOIs, which Zenodo writes automatically on a "new version" deposit | `relatedIdentifiers` |

## 3. What the 0.1.0 version ruling does to the DOIs

The cover letter asks the editor to rule on whether a `0.1.0` release satisfies the
Software Notes expectation of a stable, citable version. The metadata route depends
on that ruling, and the two outcomes are not symmetric:

* **If `0.1.0` stands as the version of record.** A published Zenodo version is
  immutable once cited: its version DOI and its DataCite metadata cannot be
  corrected in place. Then either (i) the four documents that assert "the deposit
  metadata records both authors and a 30 September release" are wrong in print, and
  the availability section must add the disclosure sentence — the archived record
  carries the short title, one creator and the 29 September issue date, while its
  licence (BSD-3-Clause) *is* registered, and the author list of record is the one in
  `CITATION.cff`; or (ii) the authors accept the disclosure *and* additionally
  publish a metadata-corrected version of the same code, which fixes the citable
  object for future readers while the note still cites `0.1.0`. Option (i) alone is
  the honest minimum; option (ii) is the better record.
* **If the editor rules that `0.1.0` is not enough.** Then the fix is a re-tag
  (`1.0.0`), which necessarily mints a **new version DOI under the same concept
  DOI**. In that case these local files change together with the deposit, in this
  order: `.zenodo.json` (`version`, `publication_date`), the tag and release,
  `CITATION.cff` (`version`, `date-released`, the `identifiers:` entry for the new
  version DOI — the concept DOI entry stays), `CHANGELOG.md` / `CHANGELOG_CN.md`,
  the version banners in `README.md` / `README_CN.md`, and only then the
  manuscript's version and date sentences and its self-citation entry. Do not
  delete or de-publish the `0.1.0` record: existing citations resolve through it,
  and the new version's `IsPreviousVersionOf` link is what keeps the lineage
  legible.

Either way, **never** present the concept DOI as the pinned citation: it always
resolves to the newest version of the series. A note that must pin the code the
paper shipped with names a version DOI, and `CITATION.cff` says so in its own
`identifiers:` descriptions.

## 4. Order of operations on Zenodo

1. Re-run section 1 and save the output — that is the "before" evidence.
2. Open the `0.1.0` record and use **New version** (this preserves the lineage and
   mints a new version DOI; the alternative, editing the draft of a record that has
   not yet been published, only exists if you are preparing a fresh deposit).
   Confirm the new draft inherits files, then correct each field in section 2.
   Zenodo's own licence picker is the only acceptable way to set the licence; typing a
   licence name into the description does not register a right. The `0.1.0` record
   already carries `bsd-3-clause` through the picker, so a "new version" deposit only
   has to keep it.
3. Fill **Funding** from the journal's Funder Registry entries, attaching the grant
   numbers exactly as they appear on the title page. A wrong number in a registry
   is worse than an empty one.
4. Publish, then re-run section 1 on the **new** version DOI and on the concept DOI,
   and confirm: two creators, no `Supervisor` contributor, the long title, the
   `rightsList` entry still reading BSD-3-Clause, `Issued` = the release day, three
   `fundingReferences`.
5. Only after step 4 agrees with the repository, the four documents' sentences about
   "the deposit metadata" are true. If the authors choose to leave the `0.1.0` record
   as published instead, add the disclosure sentence from section 3 to the
   availability section and to the cover letter in the same revision.

## 5. What this file deliberately does not do

It does not change the registry, and it does not restate author names, ORCID iDs or
grant numbers here: the review archive ships this file, and the archive's identity
guard rejects author-identifying literals in any reviewer-visible document. Those
literals live in `CITATION.cff`, `.zenodo.json`, `LICENSE` and `NOTICE.md`, which the
packager scrubs before the archive is built — which is precisely why the local files
must be the single source of truth this checklist points at.

## 6. See also

* `CITATION.cff` / `.zenodo.json` — the corrected local metadata.
* `CHANGELOG.md` — the version of record and its release date.
* `docs/CI_WORKFLOWS.md` — what is verified automatically, and what never can be.
* `docs/WORD_COUNT.md` — the companion "publish the rule, not just the number" note.
