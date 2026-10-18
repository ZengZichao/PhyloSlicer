"""Identity literals the review-archive guard must remove or reject.

This module lives in ``tools/``, which the packager deliberately excludes from the
reviewer archive.  It used to sit in ``tests/test_review_archive.py``, which IS
shipped: the fragment-spelling that kept the guard from tripping over its own
source also let a reviewer reassemble both authors’ names, e-mails, ORCID iDs,
institution, GitHub login, thesis folder and grant numbers from a file handed to
them.  Keeping the table out of the shipped tree is the fix; the test imports it
lazily and skips when ``tools/`` is absent, as it already does for every other
packager-driven case.
"""


def _frag(*parts: str) -> str:
    """Join fragments so this file neither hits the guard nor is scrubbed."""
    return "".join(parts)


IDENTITY_PATTERNS = {
    "author_surname": _frag("Ze", "ng"),
    "author_given_name": _frag("Zi", "chao"),
    "author_github_login": _frag("Ze", "ng", "Zi", "chao"),
    "author_local_account": _frag("ze", "ng", "zich", "ao"),
    "author_email": _frag("ze", "ng", "zich", "ao", "@", "sj", "tu", ".edu.cn"),
    "coauthor_email": _frag("wyz", "@", "sj", "tu", ".edu.cn"),
    "email_domain": _frag("sj", "tu", ".", "edu", ".cn"),
    "coauthor_surname": _frag("Wa", "ng"),
    "coauthor_given_name": _frag("Yin", "zha", "o"),
    "author_cjk_name": _frag("曾", "子", "超"),
    "institution": _frag("Sha", "nghai", " Jia", "o To", "ng"),
    "institution_tail": _frag("Jia", "o To", "ng"),
    "institution_cjk": _frag("上海", "交通", "大学"),
    "author_orcid": _frag("0000-0001-655", "3-", "970X"),
    "coauthor_orcid": _frag("0000-0001-871", "5-", "1386"),
    "zenodo_record": _frag("2305", "7620"),
    "zenodo_version_record": _frag("2305", "8460"),
    "grant_1": _frag("4242", "2209"),
    "grant_2": _frag("4227", "2354"),
    "grant_3": _frag("2023YFC", "3108600"),
    "home_users_dir": _frag("/", "Users", "/"),
    "checkout_path": _frag("博士", "课题"),
    "codeowners_ref": _frag("@Ze", "ng", "Zi", "chao"),
}
