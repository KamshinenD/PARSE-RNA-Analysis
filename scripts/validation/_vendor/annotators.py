"""Shared by the consensus scripts so they are self-contained; two external imports
are replaced by the local shim below.

Canonical residue/pair ids + per-annotator ingestion for the referee panel.

PI directive (2026-08-17): every annotator's output is converted into ONE
canonical representation at ingestion, every time it is read:

  - :class:`ResidueId` — ``(chain, number, icode)`` with the insertion code
    held SEPARATE from the integer residue number (PARSE ids like
    ``A-G-190K`` carry it appended to the number; BGSU Hub unit ids carry it
    in field 8; both are split apart here).
  - :data:`PairKey` — order-normalized: the lexicographically smaller
    ResidueId always comes first, and the LW label is flipped
    (``cWH`` <-> ``cHW``) whenever the listed direction is reversed.

One ingestion function per annotator (:func:`parse_pairs`,
:func:`fr3d_pairs`, :func:`tool_pairs` for rnaview / barnaba / rnapolis /
bpnet / fr3d_local) returns canonical records only; NO consumer parses raw
annotator ids itself. Ids that cannot be canonicalized raise
:class:`ValueError` — ingestion fails loudly, never silently drops a pair.

Vote tallying stays in the scripts (each adjudicator tallies differently);
notebooks read the scripts' JSON outputs under ``data/validation/``.

Annotator data files (see ``scripts/validation/fetch_fr3d_basepairs.py`` and
``run_third_annotators.py`` for provenance): published FR3D lives as one
BGSU-Hub CSV per structure under ``data/comparison/fr3d/``; every other tool
is a single ``data/comparison/{tool}_pairs.json`` whose number fields carry
any insertion code appended (``"34L"``), per DATA_GUIDE.md.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any

# --- vendored dependencies (were parse_analysis.loaders / parse_analysis.paths) ---
from pathlib import Path as _Path
from types import SimpleNamespace as _NS

_DATA = _Path(__file__).resolve().parents[3] / "data"


def _data_file(rel: str) -> _Path:
    path = _DATA / rel
    if not path.exists():
        raise FileNotFoundError(path)
    return path


paths = _NS(data_file=_data_file)
loaders = _NS(LW_COMPARISON_CLASSES=("cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
                                     "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"))
# ---------------------------------------------------------------------------------

#: Class marker for "tool reports the contact but assigns no clean LW class"
#: (RNAView's '.'/'?' edges, unmapped BPNET types, ...).
UNCLASSIFIED = "?"

#: Tools stored as a single ``data/comparison/{tool}_pairs.json``.
JSON_TOOLS = ("rnaview", "barnaba", "rnapolis", "bpnet", "fr3d_local")

#: Reasons meaning "PARSE rejected DSSR's candidate pair outright" (vs. a
#: class disagreement or a PARSE-only detection).
REJECTION_REASONS = frozenset(
    {"no_strong_hbond", "stacking", "not_candidate", "not_selected", "adjacent"}
)

_NUMBER = re.compile(r"^(-?\d+)([A-Za-z]?)$")


@dataclass(frozen=True, order=True)
class ResidueId:
    """Canonical residue identity: multi-char-safe chain, int number, icode.

    The insertion code is a separate field, never folded into ``number``
    (``A-G-190K`` -> ``ResidueId("A", 190, "K")``). Ordering is the tuple
    order ``(chain, number, icode)``, which defines the canonical pair-key
    direction.
    """

    chain: str
    number: int
    icode: str = ""


#: An order-normalized pair key: ``key[0] <= key[1]`` always; labels stored
#: with the key are expressed in the ``key[0] -> key[1]`` direction.
PairKey = tuple[ResidueId, ResidueId]

#: One structure's annotations: order-normalized pair key -> the LW classes
#: the tool assigned it (a set: some tools emit multiple records per pair).
PairSets = dict[PairKey, set[str]]


def flip(lw: str) -> str:
    """Return the LW class as seen from the other residue (``cWS`` -> ``cSW``).

    Near-pair labels keep their ``n`` prefix (``ncWS`` -> ``ncSW``); labels
    that are not 3 characters after any prefix are returned unchanged.
    """
    if lw.startswith("n"):
        return "n" + flip(lw[1:])
    if len(lw) != 3:
        return lw
    return lw[0] + lw[2] + lw[1]


def flip_label(lw: str) -> str:
    """Flip an LW label (or each candidate of a pipe-joined ambiguous label)."""
    if lw == UNCLASSIFIED:
        return lw
    return "|".join(flip(part) for part in lw.split("|"))


def residue(chain: str, number: str, icode: str = "") -> ResidueId:
    """Canonicalize one residue id, splitting any appended insertion code.

    Args:
        chain: Chain id (multi-char safe); must be non-empty.
        number: Residue number, optionally with the insertion code appended
            (``"190K"``, ``"-3"``).
        icode: Explicit insertion code, for annotators that keep it in its
            own field. Mutually exclusive with an appended one.

    Raises:
        ValueError: If the id cannot be canonicalized (empty chain,
            non-numeric number, two insertion codes, malformed icode).
    """
    match = _NUMBER.match(number)
    if not chain or match is None:
        raise ValueError(f"uncanonicalizable residue id: chain={chain!r} number={number!r}")
    appended = match.group(2)
    if appended and icode:
        raise ValueError(f"two insertion codes for chain={chain!r} number={number!r}: "
                         f"{appended!r} and {icode!r}")
    icode = icode or appended
    if icode and not icode.isalpha():
        raise ValueError(f"malformed insertion code {icode!r} for chain={chain!r} "
                         f"number={number!r}")
    return ResidueId(chain, int(match.group(1)), icode)


def parse_res_id(res_id: str) -> ResidueId:
    """Canonicalize a PARSE ``chain-resname-number`` id (``"A-G-190K"``).

    The first token is the chain, the second the residue name, and everything
    after rejoins as the number — so negative numbers (``"A-G--3"``) and
    insertion-code suffixes (``"A-G-190K"``) survive.

    Raises:
        ValueError: If the id does not have the 3-token shape or the number
            cannot be canonicalized.
    """
    parts = res_id.split("-")
    if len(parts) < 3:
        raise ValueError(f"uncanonicalizable PARSE res_id: {res_id!r}")
    return residue(parts[0], "-".join(parts[2:]))


def number_string(rid: ResidueId) -> str:
    """The canonical number-field string for the pair-list JSONs (``"190K"``)."""
    return f"{rid.number}{rid.icode}"


def pair_entry(res1: ResidueId, res2: ResidueId, label: str) -> tuple[PairKey, str]:
    """Order-normalize one directed annotation into (PairKey, label).

    The lexicographically smaller ResidueId comes first; the label is flipped
    when the listed direction is reversed, so it is always expressed in the
    ``key[0] -> key[1]`` direction.
    """
    if res2 < res1:
        return (res2, res1), flip_label(label)
    return (res1, res2), label


def lookup(pairs: PairSets, res1: ResidueId, res2: ResidueId) -> set[str]:
    """All classes a tool assigned one residue pair, in ``res1 -> res2`` direction."""
    if res2 < res1:
        return {flip_label(lw) for lw in pairs.get((res2, res1), set())}
    return set(pairs.get((res1, res2), set()))


def _ingest(rows: Iterator[tuple[ResidueId, ResidueId, str]]) -> PairSets:
    """Collect directed (res1, res2, label) rows into order-normalized PairSets."""
    pairs: PairSets = defaultdict(set)
    for res1, res2, label in rows:
        key, canonical_label = pair_entry(res1, res2, label)
        pairs[key].add(canonical_label)
    return dict(pairs)


def fr3d_pairs() -> dict[str, PairSets]:
    """Load every published-FR3D basepair CSV under ``data/comparison/fr3d/``.

    Source: one BGSU RNA 3D Hub CSV per structure (``"unit_id","lw","unit_id"``
    rows, both directions listed, true pairs only — the Hub feed excludes
    near pairs). Unit ids are ``pdb|model|chain|comp|number[|atom|alt|icode
    |sym]``; the insertion code is field 8. Model-1 units only. An empty CSV
    yields an empty PairSets (the Hub knows the structure but reports no
    basepairs), whereas a structure with no CSV at all is absent from the
    result — absence of the file is not evidence of absence of pairs.

    Raises:
        ValueError: On any non-blank row that cannot be canonicalized.
    """
    out: dict[str, PairSets] = {}
    for csv in sorted(paths.data_file("comparison/fr3d").glob("[0-9A-Z]*.csv")):
        out[csv.stem] = _ingest(_hub_rows(csv.name, csv.read_text(errors="replace")))
    return out


def _hub_rows(name: str, text: str) -> Iterator[tuple[ResidueId, ResidueId, str]]:
    """Canonical rows of one Hub CSV; malformed non-blank lines raise."""
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = [part.strip().strip('"') for part in line.split('","')]
        if len(fields) != 3:
            raise ValueError(f"{name}: malformed Hub row {line!r}")
        res1, res2 = _hub_unit(name, fields[0]), _hub_unit(name, fields[2])
        if res1 is not None and res2 is not None:
            yield res1, res2, fields[1]


def _hub_unit(name: str, unit: str) -> ResidueId | None:
    """Canonicalize one Hub unit id; None for non-model-1 units."""
    fields = unit.split("|")
    if len(fields) < 5:
        raise ValueError(f"{name}: malformed Hub unit id {unit!r}")
    if fields[1] != "1":
        return None
    icode = fields[7] if len(fields) >= 8 else ""
    return residue(fields[2], fields[4], icode)


def tool_pairs(tool: str) -> dict[str, PairSets]:
    """Load one locally-run annotator's pair sets from its single JSON.

    Source: ``data/comparison/{tool}_pairs.json`` written by
    ``scripts/validation/run_third_annotators.py`` (number fields carry any
    insertion code appended: ``"34L"``). Structures whose run did not succeed
    (any ``status != "ok"``: a specific ``fail_*`` diagnosis, a timeout, ...)
    are absent from the result, so per-record panel
    sizes reflect true coverage. RNAView and barnaba are icode-blind (their
    output formats drop insertion codes); their pairs always have
    ``icode == ""``.

    Raises:
        ValueError: If ``tool`` is not a known JSON-backed annotator, or any
            pair row cannot be canonicalized.
    """
    if tool not in JSON_TOOLS:
        raise ValueError(f"Unknown annotator {tool!r}; expected one of {JSON_TOOLS}")
    with paths.data_file(f"comparison/{tool}_pairs.json").open() as handle:
        raw = json.load(handle)
    out: dict[str, PairSets] = {}
    for pdb, record in raw.items():
        if record["status"] != "ok":
            continue
        out[pdb] = _ingest(
            (residue(c1, n1), residue(c2, n2), lw or UNCLASSIFIED)
            for c1, n1, c2, n2, lw in record["pairs"]
        )
    return out


def parse_pairs() -> dict[str, PairSets]:
    """Load PARSE's own pair sets (``data/comparison/parse_pairs.json``).

    Ambiguous calls keep their pipe-joined candidate label (``"tHH|tWW"``);
    consumers detect ambiguity as ``"|" in label``.

    Raises:
        ValueError: If a row's is_ambiguous flag contradicts its label shape,
            or any id cannot be canonicalized.
    """
    with paths.data_file("comparison/parse_pairs.json").open() as handle:
        raw = json.load(handle)
    out: dict[str, PairSets] = {}
    for pdb, record in raw.items():
        if record["status"] != "ok":
            continue
        out[pdb] = _ingest(_parse_rows(pdb, record["pairs"]))
    return out


def _parse_rows(pdb: str, rows: list[Any]) -> Iterator[tuple[ResidueId, ResidueId, str]]:
    """Canonical rows of one parse_pairs.json record; flag mismatches raise."""
    for c1, n1, c2, n2, lw, ambiguous in rows:
        if bool(ambiguous) != ("|" in (lw or "")):
            raise ValueError(f"{pdb}: is_ambiguous={ambiguous!r} contradicts label {lw!r}")
        yield residue(c1, n1), residue(c2, n2), lw or UNCLASSIFIED


def vote_sides(disputed: str, candidates: tuple[str, ...]) -> tuple[set[str], set[str]] | None:
    """Direction-agnostic label sets supporting DSSR vs PARSE for one record.

    A referee label supports DSSR if it equals the disputed class or its
    flip; it supports PARSE if it equals a PARSE candidate (or its flip)
    that is NOT direction-equivalent to the disputed class. Returns
    ``(dssr_side, parse_side)``, or None when every PARSE candidate is
    direction-equivalent to the disputed class — a DEGENERATE record that
    must be excluded from the class-disagreement tally.
    """
    dssr = {disputed, flip(disputed)}
    parse = {v for cand in candidates if cand not in dssr for v in (cand, flip(cand))}
    if not parse:
        return None
    return dssr, parse


@dataclass(frozen=True)
class DisagreementRecord:
    """One deduplicated PARSE-vs-DSSR disagreement, tagged for adjudication.

    Attributes:
        pdb_id: 4-character upper-case PDB id.
        res1: Canonical ResidueId of the first partner (file order).
        res2: Canonical ResidueId of the second partner (file order).
        source: ``"dssr_only"`` or ``"python_only"``.
        reason: The comparison pipeline's raw reason string.
        bucket: ``"rejected_by_parse"``, ``"class_disagreement"``,
            ``"dssr_partial_class"``, ``"novel"``, or ``"other"``.
        disputed_class: The LW class under dispute (the comparison file the
            record came from; per-record class fields are unpopulated).
        parse_candidates: ALL of PARSE's alternative classes for
            ``class_disagreement`` records (every candidate of
            ``python_assigns:X|Y``), else empty.
    """

    pdb_id: str
    res1: ResidueId
    res2: ResidueId
    source: str
    reason: str
    bucket: str
    disputed_class: str
    parse_candidates: tuple[str, ...]


def _bucket(reason: str) -> str:
    """Map a raw reason string onto its adjudication bucket."""
    if reason in REJECTION_REASONS:
        return "rejected_by_parse"
    if reason.startswith("python_assigns:"):
        return "class_disagreement"
    if reason.startswith("dssr_assigns:"):
        return "dssr_partial_class"
    if reason == "novel":
        return "novel"
    return "other"


def disagreement_records() -> Iterator[DisagreementRecord]:
    """Iterate every deduplicated disagreement record across all 12 classes.

    Source: the ``differences`` lists of ``data/comparison/{class}.json``
    (never the files' stale ``summary`` counters — ledger C30). Records are
    deduplicated on (pdb, unordered residue pair, source) WITHIN a class
    file; the same key appearing in two different class files would be two
    distinct disputes collapsing silently, so that raises instead. 123,022
    records survive dedup in the archived snapshot.

    Raises:
        ValueError: If a dedup key recurs across different class files, or
            any res_id cannot be canonicalized.
    """
    seen: dict[tuple[str, frozenset[str], str], str] = {}
    for lw_class in loaders.LW_COMPARISON_CLASSES:
        with paths.data_file(f"comparison/{lw_class}.json").open() as handle:
            differences = json.load(handle)["differences"]
        for row in differences:
            key = (row["pdb_id"], frozenset((row["res_id1"], row["res_id2"])), row["source"])
            if key in seen:
                if seen[key] != lw_class:
                    raise ValueError(f"dedup collision across class files "
                                     f"{seen[key]}/{lw_class}: {key}")
                continue
            seen[key] = lw_class
            yield _record(row, lw_class)


def _record(row: dict[str, Any], lw_class: str) -> DisagreementRecord:
    """Canonicalize one raw differences row."""
    reason = row.get("reason") or "?"
    candidates: tuple[str, ...] = ()
    if reason.startswith("python_assigns:"):
        candidates = tuple(reason.split(":", 1)[1].split("|"))
    return DisagreementRecord(
        pdb_id=row["pdb_id"].upper()[:4],
        res1=parse_res_id(row["res_id1"]),
        res2=parse_res_id(row["res_id2"]),
        source=row["source"],
        reason=reason,
        bucket=_bucket(reason),
        disputed_class=lw_class,
        parse_candidates=candidates,
    )
