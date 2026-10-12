"""Offline runtime adapters, adapted from the project's archived annotator runner.

Retains its chain aliases, grouped runs, output validation, and residue mapping.
All inputs come from LOCAL_CIFS; barnaba always uses those same coordinates via
synthetic PDB conversion. No download or annotation-result cache is consulted.
Installation paths are configured by run_one.py in each fresh worker.
"""
import argparse
import json
import math
import os
import re
import shutil
import string
import subprocess
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
COMPARISON_DIR = REPO / "data" / "comparison"
LOCAL_CIFS = Path(".")
RNAVIEW_HOME = Path.home() / "local" / "installs" / "RNAView"
RNAVIEW_BIN = RNAVIEW_HOME / "bin" / "rnaview"
PDB_CACHE = Path(tempfile.gettempdir()) / "parse_barnaba_pdbs"

RV_LINE = re.compile(
    r"^\s*\S+,\s+(\S+):\s+(-?\d+\w?)\s+(\S)-(\S)\s+(-?\d+\w?)\s+(\S+):\s+(\S)/(\S)\s+(cis|tran)"
)

#: Single-char chain aliases RNAView can carry through its .out format.
ALIAS_POOL = string.ascii_uppercase + string.digits + string.ascii_lowercase


def candidate_ids() -> list[str]:
    ids = set()
    for f in sorted(COMPARISON_DIR.glob("[ct][WHS][WHS].json")):
        for r in json.load(open(f)).get("differences", []):
            ids.add(r["pdb_id"].upper()[:4])
    return sorted(ids)


def get_cif(pdb: str, workdir: Path) -> Path | None:
    local = LOCAL_CIFS / f"{pdb}.cif"
    return local if local.is_file() else None


class TooManyChains(Exception):
    """More nucleic chains than single-char aliases; RNAView cannot cover."""


def _atom_site_tags(lines: list[str]) -> list[str]:
    """The _atom_site column names, in loop order."""
    tags: list[str] = []
    for line in lines:
        if line.startswith("_atom_site."):
            tags.append(line.split(".", 1)[1].strip())
        elif tags and line.startswith(("ATOM", "HETATM")):
            break
    return tags


def _atom_site_scan(lines: list[str]) -> tuple[int, int, list[str], set[str]]:
    """(auth_asym col, label_atom col, all chains in row order, nucleic chains)."""
    tags = _atom_site_tags(lines)
    chains: list[str] = []
    nucleic: set[str] = set()
    auth_i = atom_i = -1
    for line in lines:
        if tags and line.startswith(("ATOM", "HETATM")):
            if auth_i < 0:
                auth_i, atom_i = tags.index("auth_asym_id"), tags.index("label_atom_id")
            tok = line.split()
            if len(tok) != len(tags):
                raise ValueError(f"atom_site row has {len(tok)} tokens, expected {len(tags)}")
            if tok[auth_i] not in chains:
                chains.append(tok[auth_i])
            if tok[atom_i].strip('"') == "C1'":
                nucleic.add(tok[auth_i])
    return auth_i, atom_i, chains, nucleic


def _build_alias(chains: list[str], nucleic: set[str]) -> dict[str, str]:
    """Injective single-char aliases; nucleic chains have priority on the pool.

    Non-nucleic chains beyond the pool share one overflow char that is never
    a nucleic alias (they cannot appear in base-pair output).
    """
    used: set[str] = set()
    alias: dict[str, str] = {}
    for group in (sorted(nucleic), [c for c in chains if c not in nucleic]):
        for chain in group:
            if len(chain) == 1 and chain not in used:
                free = chain
            else:
                free = next((ch for ch in ALIAS_POOL if ch not in used), None)
            if free is None:
                if chain in nucleic:
                    raise TooManyChains(f"{len(nucleic)} nucleic chains > pool")
                nucleic_aliases = {alias[c] for c in nucleic}
                free = next(ch for ch in ALIAS_POOL if ch not in nucleic_aliases)
            alias[chain] = free
            used.add(free)
    return alias


def chain_groups(names: list[str], size: int = len(ALIAS_POOL)) -> list[list[str]]:
    """Chain groups of <= size such that every PAIR of names shares a group.

    One group when everything fits the alias pool; otherwise names are
    partitioned into blocks of size//2 and every pair of blocks becomes one
    group, so any two chains co-occur in the group joining their blocks.
    The covering property is verified programmatically and raises if broken.
    """
    if len(names) <= size:
        return [list(names)]
    half = size // 2
    blocks = [names[i:i + half] for i in range(0, len(names), half)]
    groups = [blocks[i] + blocks[j] for i, j in combinations(range(len(blocks)), 2)]
    uncovered = {frozenset(p) for p in combinations(names, 2)}
    for group in groups:
        uncovered -= {frozenset(p) for p in combinations(group, 2)}
    if uncovered:
        raise AssertionError(f"chain grouping left {len(uncovered)} pairs uncovered")
    return groups


def _flip_lw(lw: str | None) -> str | None:
    """LW label as seen from the other residue (cWS -> cSW)."""
    return lw if lw is None or len(lw) != 3 else lw[0] + lw[2] + lw[1]


def _dedupe_pairs(pairs: list[list]) -> list[list]:
    """Drop duplicate rows across group runs, on order-normalized keys."""
    seen = set()
    out = []
    for c1, n1, c2, n2, lw in pairs:
        if (c2, n2) < (c1, n1):
            key = (c2, n2, c1, n1, _flip_lw(lw))
        else:
            key = (c1, n1, c2, n2, lw)
        if key not in seen:
            seen.add(key)
            out.append([c1, n1, c2, n2, lw])
    return out


def alias_cif_chains(text: str) -> tuple[str, dict[str, str]] | None:
    """Rewrite auth_asym_id to single-char aliases; None if not needed.

    Returns (rewritten cif text, reverse map alias -> original chain).

    Raises:
        TooManyChains: If the nucleic chains alone exceed the alias pool.
        ValueError: If an atom_site row cannot be tokenized cleanly.
    """
    lines = text.splitlines()
    auth_i, _atom_i, chains, nucleic = _atom_site_scan(lines)
    if auth_i < 0 or all(len(c) == 1 for c in chains):
        return None
    alias = _build_alias(chains, nucleic)
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            tok = line.split()
            if len(tok) > auth_i and tok[auth_i] in alias:
                tok[auth_i] = alias[tok[auth_i]]
                line = " ".join(tok)
        out.append(line)
    reverse = {a: c for c, a in alias.items() if c in nucleic or len(c) == 1}
    return "\n".join(out) + "\n", reverse


#: Ring atoms RNAView uses to recognise a nucleotide, and its distance cutoffs.
NT_RING = ("N1", "C2", "C6")
NT_CUTOFFS = ((0, 1, 2.0), (0, 2, 2.0), (1, 2, 3.0))


def _rnaview_nucleotides(lines: list[str], cols: dict[str, int]) -> list[tuple]:
    """Residue keys RNAView calls a nucleotide, in file order.

    Mirrors residue_ident() in RNAView's fpair_sub.c: a residue is a base when
    it carries N1, C2 and C6 with d(N1,C2) <= 2.0, d(N1,C6) <= 2.0 and
    d(C2,C6) <= 3.0 — a purely geometric test that also accepts ligands with a
    pyridine ring such as LLP (pyridoxal-lysine).
    """
    order: list[tuple] = []
    ring: dict[tuple, dict[str, tuple[float, float, float]]] = {}
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        tok = line.split()
        key = (tok[cols["auth_asym_id"]], tok[cols["auth_seq_id"]],
               tok[cols["auth_comp_id"]], tok[cols["pdbx_PDB_ins_code"]])
        if key not in ring:
            ring[key] = {}
            order.append(key)
        name = tok[cols["label_atom_id"]].strip('"')
        if name in NT_RING and name not in ring[key]:
            ring[key][name] = (float(tok[cols["Cartn_x"]]), float(tok[cols["Cartn_y"]]),
                               float(tok[cols["Cartn_z"]]))
    out = []
    for key in order:
        atoms = ring[key]
        if len(atoms) == len(NT_RING) and all(
            math.dist(atoms[NT_RING[i]], atoms[NT_RING[j]]) <= cut
            for i, j, cut in NT_CUTOFFS
        ):
            out.append(key)
    return out


def drop_lone_nucleotides(text: str) -> str | None:
    """The cif without residues RNAView sees as a one-nucleotide chain.

    RNAView strips every chain left holding a single nucleotide ("rid of
    ligand", rnaview.c:800) and rebuilds its residue table, but does NOT
    rebuild chain_idx; rna_dna_bb_angle() then indexes the shortened table
    through the stale chain_idx and segfaults (4NIA, 5X6B, 6YL5, 7PD3, 7M2T).
    Removing those residues from the input cannot change the base pairs
    RNAView reports, because RNAView itself discards them before any pairing
    analysis runs. Chains are RNAView's own: maximal runs of one auth_asym_id
    over the nucleotide residues in file order; removal is iterated because
    dropping a run can merge its neighbours. Returns None when nothing to
    drop (no second attempt is warranted).
    """
    lines = text.splitlines()
    tags = _atom_site_tags(lines)
    needed = ("auth_asym_id", "auth_seq_id", "auth_comp_id", "pdbx_PDB_ins_code",
              "label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if any(t not in tags for t in needed):
        return None
    cols = {t: tags.index(t) for t in needed}
    nucleotides = _rnaview_nucleotides(lines, cols)
    drop: set[tuple] = set()
    while True:
        runs: list[list[tuple]] = []
        for key in nucleotides:
            if key in drop:
                continue
            if runs and runs[-1][-1][0] == key[0]:
                runs[-1].append(key)
            else:
                runs.append([key])
        lone = {run[0] for run in runs if len(run) == 1}
        if not lone:
            break
        drop |= lone
    if not drop:
        return None
    keep = [
        line for line in lines
        if not (line.startswith(("ATOM", "HETATM")) and _residue_key(line, cols) in drop)
    ]
    return "\n".join(keep) + "\n"


def _residue_key(line: str, cols: dict[str, int]) -> tuple:
    tok = line.split()
    return (tok[cols["auth_asym_id"]], tok[cols["auth_seq_id"]],
            tok[cols["auth_comp_id"]], tok[cols["pdbx_PDB_ins_code"]])


def _parse_rnaview_out(text: str, reverse: dict[str, str] | None) -> list[list]:
    """Extract [chain1, num1, chain2, num2, lw] rows from an RNAView .out."""
    pairs, in_bp = [], False
    for line in text.splitlines():
        if line.startswith("BEGIN_base-pair"):
            in_bp = True
            continue
        if line.startswith("END_base-pair"):
            break
        if not in_bp or "stack" in line:
            continue
        m = RV_LINE.match(line)
        if not m:
            continue
        c1, n1, _b1, _b2, n2, c2, e1, e2, ct = m.groups()
        if e1 in "+-" and e2 in "+-":
            lw = "cWW"
        elif e1 in "WHS" and e2 in "WHS":
            lw = ("c" if ct == "cis" else "t") + e1 + e2
        else:
            lw = None  # '.', '?' unidentified edge: detected, unclassified
        c1, c2 = c1.rstrip(":"), c2.rstrip(":")
        if reverse is not None:
            c1, c2 = reverse[c1], reverse[c2]  # KeyError = alias leak, fail loudly
        pairs.append([c1, n1, c2, n2, lw])
    return pairs


def _run_rnaview(work: Path, target: Path) -> tuple[str | None, str, bool]:
    """(RNAView .out text or None, status, died-by-signal); output is judged, not exit code.

    RNAView frees its residue table with stale bounds and so segfaults during
    exit cleanup on large structures LONG AFTER writing complete results;
    exit status alone therefore separates nothing. Completeness is read off
    the file instead: RNAView always closes the pair block with END_base-pair,
    including when it found zero pairs, so

      - no file, or an empty one  -> fail_rnaview_crash (died before output)
      - file without END_base-pair -> fail_rnaview_truncated (died mid-write)
      - END_base-pair present      -> ok, whatever the exit status

    An ok result that died by signal is flagged so the dirty exit is never
    silent.
    """
    try:
        proc = subprocess.run(
            [str(RNAVIEW_BIN), "--cif", target.name],
            cwd=work, capture_output=True, timeout=600,
            env={**os.environ, "RNAVIEW": str(RNAVIEW_HOME)},
        )
    except subprocess.TimeoutExpired:
        return None, "timeout", False
    signalled = proc.returncode != 0
    out = work / f"{target.name}.out"
    if not out.exists() or out.stat().st_size == 0:
        return None, "fail_rnaview_crash", signalled
    text = out.read_text(errors="replace")
    if "END_base-pair" not in text:
        return None, "fail_rnaview_truncated", signalled
    return text, "ok", signalled


def _group_cif(lines: list[str], auth_i: int, alias: dict[str, str]) -> str:
    """The cif restricted to the aliased chains (all other chains dropped)."""
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            tok = line.split()
            if tok[auth_i] not in alias:
                continue
            tok[auth_i] = alias[tok[auth_i]]
            line = " ".join(tok)
        out.append(line)
    return "\n".join(out) + "\n"


def _rnaview_group_pass(pdb: str, lines: list[str], auth_i: int, names: list[str],
                        size: int, work: Path) -> tuple[list[list], bool] | None:
    """One covering pass at a given group size; None if any group run failed."""
    pairs: list[list] = []
    signalled = False
    for g, group in enumerate(chain_groups(names, size)):
        alias = _build_alias(group, set(group))
        target = work / f"{pdb}_s{size}_g{g}.cif"
        target.write_text(_group_cif(lines, auth_i, alias))
        out_text, _status, sig = _run_rnaview(work, target)
        signalled |= sig
        if out_text is None:
            return None
        pairs += _parse_rnaview_out(out_text, {a: c for c, a in alias.items()})
    return pairs, signalled


def _ok_record(via: str | None, pairs: list[list], signalled: bool = False) -> dict:
    """An ok result, tagging the non-default route and any dirty exit."""
    record: dict = {"status": "ok"}
    if via is not None:
        record["via"] = via
    if signalled:
        record["exit"] = "signal_after_complete_output"
    record["pairs"] = pairs
    return record


def _rnaview_grouped(pdb: str, text: str, work: Path, via: str | None) -> tuple[str, dict]:
    """Overlapping chain-group runs for structures beyond the alias pool.

    RNAView segfaults on some chain subsets; a pass with any failed group is
    retried whole with half the group size, down to pairwise runs (size 3
    makes blocks of 1, i.e. every pair of chains its own run).
    """
    lines = text.splitlines()
    auth_i, _atom_i, _chains, nucleic = _atom_site_scan(lines)
    tag = "cif_groups" if via is None else f"{via}_cif_groups"
    size = len(ALIAS_POOL)
    while size >= 2:
        try:
            done = _rnaview_group_pass(pdb, lines, auth_i, sorted(nucleic), size, work)
        except KeyError:
            return pdb, {"status": "fail_alias_leak", "pairs": []}
        if done is not None:
            pairs, signalled = done
            return pdb, _ok_record(tag, _dedupe_pairs(pairs), signalled)
        size //= 2
    return pdb, {"status": "fail_rnaview_crash", "pairs": []}


def _rnaview_attempt(pdb: str, cif: Path, text: str, work: Path,
                     via: str | None) -> tuple[str, dict]:
    """One RNAView run of `text` (the verbatim cif when `via` is None)."""
    target = work / f"{pdb}.cif"
    try:
        aliased = alias_cif_chains(text)
    except TooManyChains:
        return _rnaview_grouped(pdb, text, work, via)
    except ValueError:
        return pdb, {"status": "fail_alias", "pairs": []}
    reverse = None
    if aliased is None:
        if via is None and cif != target:
            shutil.copy(cif, target)
        elif via is not None:
            target.write_text(text)
    else:
        aliased_text, reverse = aliased
        target.write_text(aliased_text)
    out_text, status, signalled = _run_rnaview(work, target)
    if out_text is None:
        return pdb, {"status": status, "pairs": []}
    try:
        pairs = _parse_rnaview_out(out_text, reverse)
    except KeyError:
        return pdb, {"status": "fail_alias_leak", "pairs": []}
    return pdb, _ok_record(via, pairs, signalled)


def rnaview_one(pdb: str) -> tuple[str, dict]:
    """RNAView on the cif, retried without RNAView's own lone-nucleotide bug.

    A crash on the verbatim file is retried on the file with lone-nucleotide
    chains removed (drop_lone_nucleotides); that route is recorded in "via".
    """
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        cif = get_cif(pdb, work)
        if cif is None:
            return pdb, {"status": "no_cif", "pairs": []}
        text = cif.read_text(errors="replace")
        result = _rnaview_attempt(pdb, cif, text, work, None)
        if result[1]["status"] == "ok":
            return result
        stripped = drop_lone_nucleotides(text)
        if stripped is None:
            return result
        return _rnaview_attempt(pdb, cif, stripped, work, "lone_nt_stripped")


BARNABA_LW = {"WC": "cWW", "GU": "cWW"}


def _barnaba_lw(ann: str) -> str | None:
    """barnaba annotation string -> canonical LW label (None if unclassified)."""
    if ann[:2] in BARNABA_LW:
        return BARNABA_LW[ann[:2]]
    if len(ann) == 3 and ann[0] in "WHS" and ann[1] in "WHS" and ann[2] in "ct":
        return ann[2] + ann[0] + ann[1]
    return None


def _barnaba_rows(chain_ids: list[str], pairings: list, res: list) -> list[list]:
    """[chain, num, chain, num, lw] rows from one bb.annotate() result."""
    rows = []
    for (i, j), ann in zip(*pairings[0], strict=True):
        _b1, n1, ci1 = res[i].rsplit("_", 2)
        _b2, n2, ci2 = res[j].rsplit("_", 2)
        rows.append([chain_ids[int(ci1)], n1, chain_ids[int(ci2)], n2, _barnaba_lw(ann)])
    return rows


def barnaba_one(pdb: str) -> tuple[str, dict]:
    warnings.filterwarnings("ignore")
    return _barnaba_rescue(pdb)


def _synth_atom_line(serial: int, at, chain: str, resfield: str, num: int, icode: str) -> str:
    """One fixed-column ATOM record; serials wrap modulo the 5-digit field."""
    aname = at.name if len(at.name) == 4 else f" {at.name:<3}"
    alt = at.altloc if at.altloc and at.altloc != "\x00" else " "
    element = at.element.name.rjust(2)[:2]
    return (
        f"ATOM  {(serial - 1) % 99999 + 1:>5} {aname[:4]}{alt}{resfield}"
        f"{chain}{num:>4}{icode}   "
        f"{at.pos.x:8.3f}{at.pos.y:8.3f}{at.pos.z:8.3f}"
        f"{at.occ:6.2f}{at.b_iso:6.2f}          {element}\n"
    )


def _write_synth_pdb(path: Path, chains: list, alias: dict[str, str]) -> dict[str, dict[int, str]]:
    """Write a nucleic-only synthetic PDB; returns per-alias renumber maps.

    Chains whose residue numbers overflow the 4-column field are renumbered
    1..n and the map (alias chain -> {new: original number string}) returned
    so parsed output can be mapped back.
    """
    num_maps: dict[str, dict[int, str]] = {}
    serial = 0
    with path.open("w") as fh:
        for name, residues in chains:
            aliased = alias[name]
            renumber = any(not -999 <= r.seqid.num <= 9999 for r in residues)
            if renumber:
                if len(residues) > 9999:
                    raise ValueError(f"chain {name}: {len(residues)} residues overflow PDB")
                num_maps[aliased] = {}
            for k, r in enumerate(residues, 1):
                num = k if renumber else r.seqid.num
                if renumber:
                    num_maps[aliased][num] = str(r.seqid.num)
                resfield = r.name[:4] if len(r.name) >= 4 else f"{r.name:>3} "
                icode = r.seqid.icode.strip() or " "
                for at in r:
                    serial += 1
                    fh.write(_synth_atom_line(serial, at, aliased, resfield, num, icode))
            fh.write("TER\n")
        fh.write("END\n")
    return num_maps


def _nucleic_chain_residues(cif: Path) -> list:
    """[(auth chain name, [gemmi residues with a C1' atom])] for model 1."""
    import gemmi

    structure = gemmi.read_structure(str(cif))
    out = []
    for chain in structure[0]:
        residues = [r for r in chain if any(a.name == "C1'" for a in r)]
        if residues:
            out.append((chain.name, residues))
    return out


def _map_back(row: list, reverse: dict[str, str], num_maps: dict) -> list:
    """Alias-chain row -> original ids; unmappable ids raise (KeyError)."""
    ca1, n1, ca2, n2, lw = row
    if ca1 in num_maps:
        n1 = num_maps[ca1][int(n1)]
    if ca2 in num_maps:
        n2 = num_maps[ca2][int(n2)]
    return [reverse[ca1], n1, reverse[ca2], n2, lw]


def _barnaba_rescue(pdb: str) -> tuple[str, dict]:
    """barnaba on synthetic PDB(s) built from the local mmCIF with gemmi."""
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        cif = get_cif(pdb, work)
        if cif is None:
            return pdb, {"status": "no_cif", "pairs": []}
        try:
            chains = _nucleic_chain_residues(cif)
        except Exception as e:  # noqa: BLE001
            return pdb, {"status": f"fail_synth_{type(e).__name__}", "error": str(e), "pairs": []}
        if not chains:
            return pdb, {"status": "no_nucleic", "pairs": []}
        import barnaba as bb
        import mdtraj as md

        groups = chain_groups([name for name, _ in chains])
        pairs: list[list] = []
        for g, group in enumerate(groups):
            alias = _build_alias(group, set(group))
            synth = work / f"{pdb}_g{g}.pdb"
            num_maps = _write_synth_pdb(synth, [c for c in chains if c[0] in group], alias)
            try:
                chain_ids = [c.chain_id for c in md.load_topology(str(synth)).chains]
                if chain_ids != [alias[name] for name in group]:
                    return pdb, {"status": "fail_mdtraj_chain_split", "pairs": []}
                _stack, pairings, res = bb.annotate(str(synth))
            except Exception as e:  # noqa: BLE001
                return pdb, {"status": f"fail_barnaba_{type(e).__name__}", "error": str(e), "pairs": []}
            reverse = {a: c for c, a in alias.items()}
            pairs += [_map_back(r, reverse, num_maps)
                      for r in _barnaba_rows(chain_ids, pairings, res)]
        via = "synthetic_pdb" if len(groups) == 1 else "synthetic_pdb_groups"
        return pdb, {"status": "ok", "via": via, "pairs": _dedupe_pairs(pairs)}


def _rnapolis_num(nt) -> str:
    """Residue number with any insertion code appended ('34L')."""
    icode = nt.icode or ""
    return f"{nt.number}{'' if icode in ('', '?') else icode.strip()}"


def rnapolis_one(pdb: str) -> tuple[str, dict]:
    warnings.filterwarnings("ignore")
    with tempfile.TemporaryDirectory() as td:
        cif = get_cif(pdb, Path(td))
        if cif is None:
            return pdb, {"status": "no_cif", "pairs": []}
        try:
            from rnapolis import annotator
            from rnapolis import parser as rp_parser

            with open(cif) as fh:
                structure = rp_parser.read_3d_structure(fh)
            base_pairs = annotator.extract_base_interactions(structure).base_pairs
        except Exception as e:  # noqa: BLE001
            return pdb, {"status": f"fail_{type(e).__name__}", "error": str(e), "pairs": []}
    pairs = []
    for bp in base_pairs:
        lw = bp.lw.value if bp.lw else None
        if lw is not None and not (
            len(lw) == 3 and lw[0] in "ct" and lw[1] in "WHS" and lw[2] in "WHS"
        ):
            lw = None
        pairs.append([str(bp.nt1.chain), _rnapolis_num(bp.nt1),
                      str(bp.nt2.chain), _rnapolis_num(bp.nt2), lw])
    return pdb, {"status": "ok", "pairs": pairs}


BPNET_BIN = Path.home() / "local" / "installs" / "bpnet" / "src" / "bpnet"
BPNET_SYSFILES = Path.home() / "local" / "installs" / "bpnet" / "sysfiles"


def _bpnet_num(num: str, ins: str) -> str:
    """Residue number with bpnet's PDB_ins_code column appended ('?' = none)."""
    return num if ins == "?" else f"{num}{ins}"


def bpnet_one(pdb: str) -> tuple[str, dict]:
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        cif = get_cif(pdb, work)
        if cif is None:
            return pdb, {"status": "no_cif", "pairs": []}
        target = work / cif.name
        if cif.parent != work:
            shutil.copy(cif, target)
        try:
            subprocess.run(
                [str(BPNET_BIN), target.name],
                cwd=work, capture_output=True, timeout=600,
                env={**os.environ, "NUCLEIC_ACID_DIR": f"{BPNET_SYSFILES}/"}, check=True,
            )
        except subprocess.TimeoutExpired:
            return pdb, {"status": "timeout", "pairs": []}
        out = work / f"{pdb}.out"
        if not out.exists():
            return pdb, {"status": "fail", "pairs": []}
        return pdb, {"status": "ok", "pairs": _parse_bpnet_out(out.read_text(errors="replace"))}


def _parse_bpnet_out(text: str) -> list[list]:
    """Extract [chain1, num1, chain2, num2, lw] rows from a bpnet .out file."""
    pairs = []
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        tok = line.split()
        # self = 5 tokens (serial, resnum, name, ins, chain); then one or
        # more partner chunks of 8 (serial, resnum, name, ins, chain,
        # edge:edge+CT, BP|TP, evalue)
        if len(tok) < 13:
            continue
        n1, c1 = _bpnet_num(tok[1], tok[3]), tok[4]
        for k in range(5, len(tok) - 7, 8):
            n2, c2 = _bpnet_num(tok[k + 1], tok[k + 3]), tok[k + 4]
            ptype = tok[k + 5]
            if ":" not in ptype:
                continue
            e1, e2ct = ptype.split(":", 1)
            e1, e2, ct = e1[:1].upper(), e2ct[:1].upper(), e2ct[-1:].upper()
            if e1 in "WHS" and e2 in "WHS" and ct in "CT":
                lw = ("c" if ct == "C" else "t") + e1 + e2
            else:
                lw = None
            pairs.append([c1, n1, c2, n2, lw])
    return pairs


def fr3d_local_one(pdb: str) -> tuple[str, dict]:
    """Local, version-pinned FR3D-Python annotation on the supplied mmCIF."""
    import sys

    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        cif = get_cif(pdb, work)
        if cif is None:
            return pdb, {"status": "no_cif", "pairs": []}
        target = work / f"{pdb}.cif"
        if cif != target:
            shutil.copy(cif, target)
        try:
            subprocess.run(
                [sys.executable, "-m", "fr3d.classifiers.NA_pairwise_interactions",
                 "-c", "basepair", "-i", str(work), "-o", str(work), f"{pdb}.cif"],
                capture_output=True, timeout=1800, check=True,
            )
        except subprocess.TimeoutExpired:
            return pdb, {"status": "timeout", "pairs": []}
        out = next((p for p in [work / f"{pdb}_basepair.txt", work / f"{pdb}_basepairs.txt"] if p.exists()), work / f"{pdb}_basepair.txt")
        if not out.exists():
            return pdb, {"status": "fail", "pairs": []}
        pairs = []
        for line in out.read_text(errors="replace").splitlines():
            tok = line.split("\t")
            if len(tok) < 3:
                continue
            f1, f2 = tok[0].split("|"), tok[2].split("|")
            if len(f1) < 5 or len(f2) < 5 or f1[1] != "1" or f2[1] != "1":
                continue
            lw = tok[1]
            if lw.startswith("n"):
                lw = "n" + lw[1] + lw[2:].upper()
            elif len(lw) == 3:
                lw = lw[0] + lw[1:].upper()
            n1 = f1[4] + (f1[7] if len(f1) >= 8 else "")  # unit-id field 8 = icode
            n2 = f2[4] + (f2[7] if len(f2) >= 8 else "")
            pairs.append([f1[2], n1, f2[2], n2, lw])
        return pdb, {"status": "ok", "pairs": pairs}
