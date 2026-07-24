"""Small helpers vendored from pair_finder so the build scripts are self-contained."""
import re

# --- from pair_finder.scoring.issues -----------------------------------------
REPORT_GROUPS = {
    "Pair In-Plane Displacement": ["shear", "stretch"],
    "Pair Non-Coplanarity":       ["stagger", "buckle"],
    "Propeller Twist":            ["propeller"],
    "Pair Opening":               ["opening"],
    "H-Bond Distance":            ["distance"],
    "H-Bond Angles":              ["hbond_angles"],
    "Incorrect H-Bond Count":     ["incorrect_hbond_count"],
}

# --- from pair_finder.core.identifiers ---------------------------------------
_SEQ_RE = re.compile(r"^(-?\d+)([A-Za-z]*)$")

def parse_res_seq(res_id):
    """Parse 'A-G-20A' -> ('A', 20, 'A'); handles negative seq nums + insertion codes."""
    if not res_id:
        return None
    parts = res_id.split("-")
    if len(parts) < 3:
        return None
    m = _SEQ_RE.match("-".join(parts[2:]))
    if not m:
        return None
    return parts[0], int(m.group(1)), m.group(2)
