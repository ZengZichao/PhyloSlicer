"""Newick & NEXUS readers/writers (lightweight, dependency-free)."""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np

from ..core.tree import TreeArray

__all__ = [
    "read_newick",
    "read_nexus",
    "write_newick",
    "write_nexus",
    "parse_newick",
    "split_newick_records",
]

_SAFE_NAME = re.compile(r"^[A-Za-z0-9_.\-|]+$")


def _strip_comments(text: str) -> str:
    """Remove square-bracket comments (incl. NHX annotations), quote-aware.

    A ``[...]`` run is only a comment when it starts outside a quoted label:
    the previous blanket ``re.sub(r"\\[[^\\]]*\\]", "")`` also ate brackets that
    are part of a species name (cultivar / strain designations such as
    ``'Beta vulgaris [N.E. Br.]'``), silently truncating the label.  Quotes are
    recognised only at token boundaries, for the reason given in
    :func:`split_newick_records`.
    """
    boundary = {"(", ")", ",", ":", ";", "'", "[", "]", ""}
    out: list[str] = []
    i = 0
    n = len(text)
    in_quote = False
    prev = ""
    while i < n:
        ch = text[i]
        if in_quote:
            out.append(ch)
            if ch == "'":
                in_quote = False
            prev = ch
            i += 1
            continue
        if ch == "'" and prev in boundary:
            in_quote = True
            out.append(ch)
            prev = ch
            i += 1
            continue
        if ch == "[":
            end = text.find("]", i + 1)
            if end == -1:  # unterminated comment: keep the text as-is
                out.append(ch)
                prev = ch
                i += 1
                continue
            i = end + 1
            prev = "]"
            continue
        out.append(ch)
        prev = ch
        i += 1
    return "".join(out)


def split_newick_records(text: str) -> list[str]:
    """Split a possibly multi-tree Newick string on top-level semicolons.

    Quoted labels (``'...'``, including the doubled-quote escape ``''``) and
    bracket comments (``[...]``) are respected, so a semicolon inside a quoted
    species name or a comment does not split a tree.  Records that carry no
    tree structure (an empty or ``;``-only chunk produced by a stray double
    separator) are dropped rather than handed to the parser.

    A ``'`` only opens a quoted label at a token boundary (after ``(``, ``,``,
    ``:``, ``)``, ``;`` or the start of the text).  That keeps bare apostrophes
    inside unquoted labels - ``a'b``, which this package's own parser accepts -
    from being mistaken for an unterminated quote.  ``[`` needs no such guard:
    a comment may legitimately start straight after a label (``) #HND[&&NHX:``).
    """
    boundary = {"(", ")", ",", ":", ";", "'", "[", "]", ""}
    records: list[str] = []
    buf: list[str] = []
    in_quote = False
    in_comment = False
    prev = ""
    for ch in text:
        buf.append(ch)
        if in_comment:
            if ch == "]":
                in_comment = False
        elif in_quote:
            if ch == "'":
                in_quote = False
        elif ch == "'" and prev in boundary:
            in_quote = True
        elif ch == "[":
            in_comment = True
        elif ch == ";":
            record = "".join(buf).strip()
            buf = []
            if record.strip(";\n\r\t "):
                records.append(record)
        prev = ch
    if in_quote:
        raise ValueError(
            'unterminated quoted label: the closing "\'" is missing, so the '
            "remaining text cannot be parsed as Newick"
        )
    tail = "".join(buf).strip()
    if tail.strip(";\n\r\t "):
        records.append(tail if tail.endswith(";") else tail + ";")
    return records


def parse_newick(text: str, label_prefix: str = "tip_") -> TreeArray:
    """Parse a single Newick string into a :class:`TreeArray`."""
    s = _strip_comments(text.strip())
    if not s.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    s = s[:-1]

    tip_labels: list[str] = []
    edges: list[tuple[int, int]] = []
    lengths: list[float] = []
    state = {"pos": 0, "internal_count": 0}
    n = len(s)

    def skip_ws() -> None:
        while state["pos"] < n and s[state["pos"]] in " \t\r\n":
            state["pos"] += 1

    def read_name() -> str:
        skip_ws()
        i = state["pos"]
        if i < n and s[i] == "'":
            state["pos"] = i + 1
            parts: list[str] = []
            while True:
                try:
                    j = s.index("'", state["pos"])
                except ValueError as exc:
                    raise ValueError(
                        f"unterminated quoted label starting at position {i}: "
                        f"{s[i : min(i + 20, n)]!r}"
                    ) from exc
                parts.append(s[state["pos"] : j])
                if j + 1 < n and s[j + 1] == "'":
                    # Newick escapes a literal quote by doubling it: 'a''b' -> a'b
                    parts.append("'")
                    state["pos"] = j + 2
                    continue
                state["pos"] = j + 1
                break
            return "".join(parts)
        start = i
        while state["pos"] < n and s[state["pos"]] not in "(),:;":
            state["pos"] += 1
        return s[start : state["pos"]].strip()

    def read_length() -> float:
        skip_ws()
        if state["pos"] < n and s[state["pos"]] == ":":
            state["pos"] += 1
            skip_ws()
            start = state["pos"]
            while state["pos"] < n and s[state["pos"]] not in "(),;":
                state["pos"] += 1
            try:
                return float(s[start : state["pos"]])
            except ValueError as exc:
                raise ValueError(f"invalid branch length near {s[start : state['pos']]!r}") from exc
        return 0.0

    def close_internal(frame: tuple[int, list]) -> tuple[int, float]:
        """Emit the edges of a finished internal node, return (id, own length)."""
        my_id, children, own_length = frame
        for child, clen in children:
            edges.append((my_id, child))
            lengths.append(clen)
        return my_id, own_length

    # Iterative shift-reduce parse of the same grammar the recursive descent
    # used (( '(' child{',' child} ')' label? length? ) | leaf ).  Recursion
    # depth equals tree depth here, so a 10^4-tip caterpillar - a perfectly
    # ordinary BEAST output - raised RecursionError.
    stack: list[tuple[int, list, float]] = []
    root_pair: tuple[int, float] | None = None
    while True:
        skip_ws()
        if state["pos"] >= n:
            raise ValueError("malformed Newick: unexpected end of string")
        ch = s[state["pos"]]
        if ch == "(":
            my_id = -(state["internal_count"] + 1)  # temp negative internal id
            state["internal_count"] += 1
            state["pos"] += 1
            stack.append((my_id, [], 0.0))
            continue
        if ch == ")":
            if not stack:
                raise ValueError(f"malformed Newick: unbalanced ')' at position {state['pos']}")
            state["pos"] += 1
            read_name()  # optional internal label / support value
            my_id, children, _ = stack[-1]
            stack.pop()
            completed = (my_id, children, read_length())
            pair = close_internal(completed)
        else:
            name = read_name()
            own_length = read_length()
            if ch == "," or ch == ";":
                raise ValueError(f"malformed Newick near position {state['pos']}")
            tip_idx = len(tip_labels)
            tip_labels.append(name if name else f"{label_prefix}{tip_idx + 1}")
            pair = (tip_idx, own_length)
        if stack:
            parent_id, children, plen = stack[-1]
            children.append(pair)
            stack[-1] = (parent_id, children, plen)
            skip_ws()
            nxt = s[state["pos"]] if state["pos"] < n else ""
            if nxt == ",":
                state["pos"] += 1
                continue
            if nxt in (")", ";"):
                continue
            raise ValueError(f"malformed Newick near position {state['pos']}")
        # completed the root node
        root_pair = pair
        break

    if root_pair is None:
        raise ValueError("malformed Newick: no tree found")
    skip_ws()
    if state["pos"] != n:
        raise ValueError(f"trailing characters after Newick tree: {s[state['pos'] :][:20]!r}")

    n_tips = len(tip_labels)
    if n_tips < 2:
        raise ValueError(
            f"a sliceable tree needs at least 2 tips, but the Newick string declares "
            f"{n_tips} ({s[:60]!r})"
        )
    if not edges:
        raise ValueError("the Newick string declares no edges")
    # temp internal ids are -1, -2, ... in creation order and the newick root
    # is always -1 (outermost parens parsed first); number internals starting
    # at n_tips with the root first (TreeArray convention).
    temp_internal = sorted({p for p, _ in edges} | {c for _, c in edges if c < 0}, reverse=True)
    remap = {temp: n_tips + i for i, temp in enumerate(temp_internal)}
    out_edges = np.asarray([[remap[p], remap.get(c, c)] for p, c in edges], dtype=np.int64)
    out_lengths = np.asarray(lengths, dtype=np.float64)
    root_edge = float(root_pair[1])
    if root_edge:
        warnings.warn(
            f"the Newick string carries a root stem of length {root_edge:g} above the "
            "root node; TreeArray follows the ape::keep.tip(root.edge = 0) convention "
            f"and excludes it from total_pd(), so it is recorded as root_edge rather "
            f"than as an edge (write_newick() restores it on output)",
            stacklevel=2,
        )
    return TreeArray(
        tip_labels=tip_labels, edges=out_edges, lengths=out_lengths, root_edge=root_edge
    )


def _looks_like_newick(text: str) -> bool:
    """Whether ``text`` carries Newick structure rather than being a path.

    The previous heuristic ("the string contains ``(``") misclassified both a
    path whose directory name contains a parenthesis and a legitimate
    bracket-free single-tip tree.  Structural characters are checked instead,
    after existing files have already been claimed by :func:`read_newick`.
    """
    return any(ch in text for ch in "();:")


def read_newick(source) -> TreeArray:
    """Read a tree from a Newick string or a file containing one tree.

    Resolution order for a ``str`` / :class:`~pathlib.Path` argument:

    1. an existing file (checked first, so paths containing ``(``, ``:`` or
       spaces work);
    2. otherwise a string carrying Newick structure (``(``, ``)``, ``:``,
       ``;``) is parsed directly;
    3. otherwise it is treated as a path and a :class:`FileNotFoundError` with
       the attempted path is raised.

    Files with multiple trees: only the first tree is returned; use
    :func:`read_nexus` (or :func:`phyloslicer.uncertainty.iter_trees`) for
    explicit multi-tree support.
    """
    if isinstance(source, Path):
        if not source.is_file():
            raise FileNotFoundError(f"tree file not found: {str(source)[:200]}")
        text = source.read_text()
    elif isinstance(source, str):
        candidate = source.strip()
        if candidate and _is_existing_file(candidate):
            text = Path(candidate).read_text()
        elif _looks_like_newick(candidate):
            text = candidate
        else:
            raise FileNotFoundError(
                f"{source!r:.200} is neither an existing file nor a Newick string "
                "(no '(', ')', ':' or ';' found); pass a path to a tree file or a "
                "Newick string"
            )
    else:
        raise TypeError("source must be a Newick string or a file path")
    _reject_nexus(text, source)
    records = split_newick_records(text)
    if not records:
        raise ValueError(f"no tree found in source: {str(source)[:80]!r}")
    return parse_newick(records[0])


def _reject_nexus(text: str, source) -> None:
    """Turn "I fed you a NEXUS file" into a named, actionable error.

    ``posterior`` accepts ``.nex`` through :func:`read_nexus`, so a NEXUS path
    reaching :func:`read_newick` is a plausible user mistake; truncating the
    header at the first ``;`` used to report "a sliceable tree needs at least 2
    tips", which names neither the format nor the remedy.
    """
    head = text.lstrip()[:400].upper()
    if head.startswith("#NEXUS") or "BEGIN TREES" in head:
        raise ValueError(
            f"{str(source)[:120]!r} is a NEXUS file, not a single Newick tree; "
            "read it with TreeArray.from_nexus(path) or "
            "phyloslicer.io.read_nexus(path), which return every TREE statement"
        )


def _is_existing_file(candidate: str) -> bool:
    """``Path.is_file`` that never raises on odd strings (NUL bytes, long names)."""
    if len(candidate) > 255 or "\x00" in candidate:
        return False
    try:
        return Path(candidate).is_file()
    except OSError:  # pragma: no cover - ENAMETOOLONG, EINVAL, ...
        return False


def _is_existing_dir(candidate: str) -> bool:
    """``Path.is_dir`` with the same never-raise guarantee as above."""
    if len(candidate) > 255 or "\x00" in candidate:
        return False
    try:
        return Path(candidate).is_dir()
    except OSError:  # pragma: no cover
        return False


def read_nexus(path) -> list[TreeArray]:
    """Read all trees from a NEXUS file's TREES block.

    Tree statements may span multiple lines: the newick string is
    accumulated until the terminating ``;``.
    """
    lines = Path(path).read_text().splitlines()
    trees: list[TreeArray] = []
    in_trees = False
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if not in_trees:
            if stripped.upper().startswith("BEGIN TREES"):
                in_trees = True
            i += 1
            continue
        if stripped.upper().startswith("END"):
            break
        if stripped.upper().startswith("TREE"):
            eq = stripped.find("=")
            if eq != -1:
                newick = stripped[eq + 1 :].strip()
                # multi-line statement: keep appending until ';'
                while not newick.endswith(";") and i + 1 < len(lines):
                    i += 1
                    newick += " " + lines[i].strip()
                newick = re.sub(r"^\[[^\]]*\]\s*", "", newick)  # e.g. "[&R]"
                if not newick.endswith(";"):
                    newick += ";"
                trees.append(parse_newick(newick))
        i += 1
    if not trees:
        raise ValueError(f"no TREE entries found in NEXUS file {path}")
    return trees


def _format_name(name: str) -> str:
    """Quote a label that needs it, escaping embedded quotes Newick-style.

    A quoted Newick label escapes ``'`` by doubling it; wrapping without
    escaping produced ``'a'b'`` for the very common botanical/strain label
    ``a'b``, which no Newick reader (including this one) can parse.
    """
    if _SAFE_NAME.match(name):
        return name
    return "'" + name.replace("'", "''") + "'"


def write_newick(tree: TreeArray, length_fmt: str = "%.12g") -> str:
    """Serialise a :class:`TreeArray` to a Newick string (root has no stem length).

    Iterative on purpose.  The previous recursive implementation blew the
    recursion limit on strongly unbalanced topologies (a 10^4-tip caterpillar
    nests 10^4 calls) and built each clade by string concatenation inside that
    recursion, which is quadratic in the node count for such trees.  Here every
    clade is emitted through an explicit token stack and the pieces are joined
    once, so the cost is linear.
    """
    children = tree.children_of()
    lengths_by_child = {int(c): float(seg) for c, seg in zip(tree.edges[:, 1], tree.lengths)}
    n_tips = tree.n_tips
    out: list[str] = []
    stack: list[object] = [(tree.root, True)]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            out.append(item)
            continue
        node, is_root = item  # type: ignore[misc]
        kids = children.get(node, [])
        if not kids:
            if node >= n_tips:
                raise ValueError(f"internal node {node} has no children; cannot serialise")
            out.append(_format_name(tree.tip_labels[node]))
            out.append(f":{length_fmt % lengths_by_child[node]}")
            continue
        out.append("(")
        # push in reverse emission order: the closing token first (so it pops
        # last), then the kids newest-first with commas between them
        stem = f":{length_fmt % tree.root_edge}" if tree.root_edge else ""
        if is_root:
            stack.append(f"){stem}")
        else:
            stack.append(f"):{length_fmt % lengths_by_child[node]}")
        for i, kid in enumerate(reversed(kids)):
            if i:
                stack.append(",")
            stack.append((int(kid), False))
    return "".join(out) + ";"


def write_nexus(trees: list[TreeArray], tree_names: list[str] | None = None) -> str:
    """Serialise multiple trees to a NEXUS string."""
    if tree_names is None:
        tree_names = [f"tree_{i + 1}" for i in range(len(trees))]
    body = "\n".join(f"    TREE {name} = {write_newick(t)}" for name, t in zip(tree_names, trees))
    return "#NEXUS\n\nBEGIN TREES;\n" + body + "\nEND;\n"
