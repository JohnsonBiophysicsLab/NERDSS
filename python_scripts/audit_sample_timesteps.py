#!/usr/bin/env python3
"""Audit pairwise-reaction timesteps in sample_inputs.

The calculation intentionally applies the 3D equation supplied for this audit
to every pairwise declaration, including 2D and 3D-to-2D examples.  Results are
therefore a literal screen against that equation, not a replacement for the
separate 2D criterion used by NERDSS.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
MOL_EXPR_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)")


@dataclass(frozen=True)
class Molecule:
    name: str
    diffusion: tuple[float, float, float]
    rotational_diffusion: tuple[float, float, float]
    coordinates: dict[str, tuple[float, float, float]]

    @property
    def d_scalar(self) -> float:
        """Scalar D supplied on the molecule's mobile axes.

        All current sample inputs use the same nonzero value on every mobile
        axis.  Averaging only those axes preserves that supplied scalar for 3D,
        membrane, and sphere-surface templates.
        """

        active = [abs(value) for value in self.diffusion if abs(value) > 1e-15]
        return sum(active) / len(active) if active else 0.0

    @property
    def dr_scalar(self) -> float:
        active = [
            abs(value)
            for value in self.rotational_diffusion
            if abs(value) > 1e-15
        ]
        return sum(active) / len(active) if active else 0.0

    def lever_arm(self, interface: str) -> float | None:
        com = self.coordinates.get("com")
        point = self.coordinates.get(interface.lower())
        if com is None or point is None:
            return None
        return math.dist(com, point)


@dataclass
class Reaction:
    input_path: Path
    line: int
    equation: str
    mol_a: str
    iface_a: str
    mol_b: str
    iface_b: str
    sigma: float | None
    sigma_source: str = "explicit"


@dataclass
class AuditRow:
    reaction: Reaction
    configured_dt: float | None
    volume: float | None
    count_a: int | None
    geometry: str
    rho_a: float | None
    d_trans: float | None
    a_a: float | None
    a_b: float | None
    d_rot: float | None
    dt_trans: float | None
    dt_rot: float | None
    trans_status: str
    rot_status: str
    issue: str = ""


def remove_comment(line: str) -> str:
    return line.split("#", 1)[0]


def split_top_level(text: str, delimiter: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    start = 0
    for index, character in enumerate(text):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif character == delimiter and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    parts.append(text[start:].strip())
    return parts


def site_name(token: str) -> str:
    return re.split(r"[~!]", token.strip(), maxsplit=1)[0].strip()


def site_state(token: str) -> str | None:
    match = re.search(r"~([^!~]+)", token)
    return match.group(1) if match else None


def parse_mol_expression(expression: str) -> tuple[str, list[str]]:
    match = MOL_EXPR_RE.fullmatch(expression.strip())
    if not match:
        # The compartment pseudo-reactant has the same surface grammar, but
        # retaining this fallback gives a useful unresolved audit row.
        name = expression.split("(", 1)[0].strip()
        return name, []
    sites = [part.strip() for part in split_top_level(match.group(2), ",") if part.strip()]
    return match.group(1), sites


def product_molecules(rhs: str) -> list[tuple[str, list[str]]]:
    return [
        (match.group(1), [
            part.strip()
            for part in split_top_level(match.group(2), ",")
            if part.strip()
        ])
        for match in MOL_EXPR_RE.finditer(rhs)
    ]


def reacting_interface(
    lhs_sites: list[str], product_sites: list[str] | None
) -> str:
    if product_sites is not None:
        product_by_name = {site_name(token): token for token in product_sites}
        for lhs_token in lhs_sites:
            name = site_name(lhs_token)
            rhs_token = product_by_name.get(name)
            if rhs_token is None:
                continue
            lhs_bound = "!" in lhs_token
            rhs_bound = "!" in rhs_token
            if (not lhs_bound and rhs_bound) or site_state(lhs_token) != site_state(rhs_token):
                return name
    # For catalytic pairwise state changes, one partner is unchanged.  Its
    # listed interface is still the encounter site.
    return site_name(lhs_sites[0]) if lhs_sites else "?"


def parse_reactions(path: Path) -> list[Reaction]:
    lines = path.read_text(errors="replace").splitlines()
    candidates: list[tuple[int, str, float | None]] = []
    in_reactions = False
    index = 0
    while index < len(lines):
        cleaned = remove_comment(lines[index]).strip()
        compact = re.sub(r"\s+", "", cleaned).lower()
        if compact == "startreactions":
            in_reactions = True
            index += 1
            continue
        if compact == "endreactions":
            in_reactions = False
            index += 1
            continue
        if not in_reactions or not cleaned or ("->" not in cleaned and "<->" not in cleaned):
            index += 1
            continue

        arrow = "<->" if "<->" in cleaned else "->"
        lhs = cleaned.split(arrow, 1)[0]
        if len(split_top_level(lhs, "+")) != 2:
            index += 1
            continue

        sigma: float | None = None
        lookahead = index + 1
        while lookahead < len(lines):
            parameter_line = remove_comment(lines[lookahead]).strip()
            parameter_compact = re.sub(r"\s+", "", parameter_line).lower()
            if parameter_compact == "endreactions" or (
                parameter_line and ("->" in parameter_line or "<->" in parameter_line)
            ):
                break
            sigma_match = re.match(rf"(?i)^sigma\s*=\s*({NUMBER})", parameter_line)
            if sigma_match:
                sigma = float(sigma_match.group(1))
            lookahead += 1
        candidates.append((index + 1, cleaned, sigma))
        index += 1

    reactions: list[Reaction] = []
    for line_number, equation, sigma in candidates:
        arrow = "<->" if "<->" in equation else "->"
        lhs, rhs = [part.strip() for part in equation.split(arrow, 1)]
        lhs_species = split_top_level(lhs, "+")
        if len(lhs_species) != 2:
            continue
        mol_a, sites_a = parse_mol_expression(lhs_species[0])
        mol_b, sites_b = parse_mol_expression(lhs_species[1])
        products = product_molecules(rhs)
        product_a = products[0][1] if len(products) >= 1 else None
        product_b = products[1][1] if len(products) >= 2 else None
        reactions.append(
            Reaction(
                input_path=path,
                line=line_number,
                equation=equation,
                mol_a=mol_a,
                iface_a=reacting_interface(sites_a, product_a),
                mol_b=mol_b,
                iface_b=reacting_interface(sites_b, product_b),
                sigma=sigma,
            )
        )

    # Conditional-rate declarations omit geometry.  NERDSS attaches them to
    # the pre-existing physical interface pair; recover that pair's sigma.
    explicit_by_pair: dict[
        tuple[tuple[str, str], tuple[str, str]], list[tuple[int, float]]
    ] = defaultdict(list)
    for reaction in reactions:
        if reaction.sigma is not None:
            explicit_by_pair[physical_pair(reaction)].append((reaction.line, reaction.sigma))
    for reaction in reactions:
        if reaction.sigma is not None:
            continue
        matches = explicit_by_pair.get(physical_pair(reaction), [])
        preceding = [item for item in matches if item[0] < reaction.line]
        if preceding:
            reaction.sigma = preceding[-1][1]
            reaction.sigma_source = f"inherited from line {preceding[-1][0]}"
        elif matches:
            reaction.sigma = min(matches, key=lambda item: abs(item[0] - reaction.line))[1]
            reaction.sigma_source = "inherited from matching interface pair"
        else:
            reaction.sigma = 1.0
            reaction.sigma_source = "NERDSS default"
    return reactions


def physical_pair(
    reaction: Reaction,
) -> tuple[tuple[str, str], tuple[str, str]]:
    pair = sorted(
        (
            (reaction.mol_a.lower(), reaction.iface_a.lower()),
            (reaction.mol_b.lower(), reaction.iface_b.lower()),
        )
    )
    return pair[0], pair[1]


def parse_vector(text: str, key: str) -> tuple[float, float, float] | None:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*\[\s*({NUMBER})\s*,\s*({NUMBER})\s*,\s*({NUMBER})\s*\]",
        text,
    )
    if not match:
        return None
    return tuple(float(match.group(index)) for index in range(1, 4))  # type: ignore[return-value]


def parse_molecule(path: Path) -> Molecule | None:
    text = path.read_text(errors="replace")
    name_match = re.search(r"(?im)^\s*name\s*=\s*([A-Za-z_][A-Za-z0-9_]*)", text)
    diffusion = parse_vector(text, "D")
    rotational_diffusion = parse_vector(text, "Dr")
    if not name_match or diffusion is None or rotational_diffusion is None:
        return None
    coordinates: dict[str, tuple[float, float, float]] = {}
    coord_re = re.compile(
        rf"^\s*([A-Za-z_][A-Za-z0-9_]*)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})(?:\s|$)"
    )
    for line in text.splitlines():
        match = coord_re.match(remove_comment(line))
        if match:
            coordinates[match.group(1).lower()] = tuple(
                float(match.group(index)) for index in range(2, 5)
            )  # type: ignore[assignment]
    return Molecule(
        name=name_match.group(1),
        diffusion=diffusion,
        rotational_diffusion=rotational_diffusion,
        coordinates=coordinates,
    )


def parse_molecules(directory: Path) -> dict[str, Molecule]:
    molecules: dict[str, Molecule] = {}
    for path in sorted(directory.glob("*.mol")):
        molecule = parse_molecule(path)
        if molecule:
            molecules[molecule.name.lower()] = molecule
    return molecules


def section_lines(text: str, section: str) -> Iterable[str]:
    active = False
    for line in text.splitlines():
        compact = re.sub(r"\s+", "", remove_comment(line)).lower()
        if compact == f"start{section.lower()}":
            active = True
            continue
        if compact == f"end{section.lower()}":
            active = False
            continue
        if active:
            yield remove_comment(line).strip()


def parse_counts(texts: Iterable[str]) -> dict[str, int | None]:
    counts: dict[str, int | None] = {}
    for text in texts:
        for line in section_lines(text, "molecules"):
            if not line:
                continue
            if ":" not in line:
                counts[line.strip().lower()] = None
                continue
            name, value = line.split(":", 1)
            values = re.findall(r"(?:^|,)\s*(\d+)", value)
            counts[name.strip().lower()] = sum(int(number) for number in values) if values else None
    return counts


def scalar_parameter(text: str, key: str) -> float | None:
    match = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*({NUMBER})", text)
    return float(match.group(1)) if match else None


def parse_timestep(text: str) -> float | None:
    return scalar_parameter(text, "timeStep")


def parse_volume(text: str) -> float | None:
    sphere = re.search(r"(?im)^\s*isSphere\s*=\s*true\b", text) is not None
    radius = scalar_parameter(text, "sphereR")
    if sphere and radius is not None:
        return 4.0 * math.pi * radius**3 / 3.0
    box = re.search(
        rf"(?im)^\s*waterBox\s*=\s*\[\s*({NUMBER})\s*,\s*({NUMBER})\s*,\s*({NUMBER})\s*\]",
        text,
    )
    if box:
        return math.prod(float(box.group(index)) for index in range(1, 4))
    return None


def allowed_timestep(rho: float, sigma: float, diffusion: float) -> float:
    if diffusion <= 0.0:
        return math.inf
    spacing = (3.0 / (4.0 * math.pi * rho) + sigma**3) ** (1.0 / 3.0) - sigma
    return spacing**2 / (56.0 * diffusion)


def geometry(mol_a: Molecule, mol_b: Molecule) -> str:
    a_mobile_z = abs(mol_a.diffusion[2]) > 1e-15
    b_mobile_z = abs(mol_b.diffusion[2]) > 1e-15
    if a_mobile_z and b_mobile_z:
        return "3D"
    if a_mobile_z or b_mobile_z:
        return "3D-to-2D"
    return "2D"


def audit_reaction(
    reaction: Reaction,
    timestep: float | None,
    volume: float | None,
    counts: dict[str, int | None],
    molecules: dict[str, Molecule],
) -> AuditRow:
    issues: list[str] = []
    mol_a = molecules.get(reaction.mol_a.lower())
    mol_b = molecules.get(reaction.mol_b.lower())
    count_a = counts.get(reaction.mol_a.lower())
    if timestep is None:
        issues.append("timeStep missing")
    if volume is None:
        issues.append("volume missing")
    if count_a is None:
        issues.append(f"N_A missing for {reaction.mol_a}")
    elif count_a == 0:
        issues.append(f"initial N_A=0 for {reaction.mol_a}")
    if mol_a is None:
        issues.append(f"{reaction.mol_a}.mol unavailable")
    if mol_b is None:
        issues.append(f"{reaction.mol_b}.mol unavailable")

    rho = count_a / volume if count_a and volume else None
    row_geometry = geometry(mol_a, mol_b) if mol_a and mol_b else "?"
    d_trans = mol_a.d_scalar + mol_b.d_scalar if mol_a and mol_b else None
    a_a = mol_a.lever_arm(reaction.iface_a) if mol_a else None
    a_b = mol_b.lever_arm(reaction.iface_b) if mol_b else None
    if mol_a and a_a is None:
        issues.append(f"interface {reaction.mol_a}.{reaction.iface_a} unavailable")
    if mol_b and a_b is None:
        issues.append(f"interface {reaction.mol_b}.{reaction.iface_b} unavailable")
    d_rot = (
        (2.0 / 3.0)
        * (a_a**2 * mol_a.dr_scalar + a_b**2 * mol_b.dr_scalar)
        if mol_a and mol_b and a_a is not None and a_b is not None
        else None
    )

    dt_trans = (
        allowed_timestep(rho, reaction.sigma, d_trans)
        if rho is not None and reaction.sigma is not None and d_trans is not None
        else None
    )
    dt_rot = (
        allowed_timestep(rho, reaction.sigma, d_rot)
        if rho is not None and reaction.sigma is not None and d_rot is not None
        else None
    )

    def status(limit: float | None) -> str:
        if timestep is None or limit is None:
            return "UNRESOLVED"
        return "TOO BIG" if timestep > limit * (1.0 + 1e-12) else "OK"

    return AuditRow(
        reaction=reaction,
        configured_dt=timestep,
        volume=volume,
        count_a=count_a,
        geometry=row_geometry,
        rho_a=rho,
        d_trans=d_trans,
        a_a=a_a,
        a_b=a_b,
        d_rot=d_rot,
        dt_trans=dt_trans,
        dt_rot=dt_rot,
        trans_status=status(dt_trans),
        rot_status=status(dt_rot),
        issue="; ".join(dict.fromkeys(issues)),
    )


def fmt(value: float | None, significant: int = 4) -> str:
    if value is None:
        return "—"
    if math.isinf(value):
        return "∞"
    if value == 0:
        return "0"
    magnitude = abs(value)
    if 1e-3 <= magnitude < 1e4:
        decimals = max(0, significant - 1 - math.floor(math.log10(magnitude)))
        return f"{value:.{decimals}f}"
    return f"{value:.{significant - 1}e}"


def compress_lines(lines: list[int]) -> str:
    values = sorted(set(lines))
    if not values:
        return "—"
    if len(values) <= 5:
        return ", ".join(str(value) for value in values)
    return f"{values[0]}–{values[-1]} ({len(values)} declarations)"


def result_summary(rows: list[AuditRow], field: str) -> str:
    statuses = [getattr(row, field) for row in rows]
    failures = statuses.count("TOO BIG")
    unresolved = statuses.count("UNRESOLVED")
    if failures:
        result = f"TOO BIG: {failures}/{len(rows)}"
    elif unresolved == len(rows):
        result = f"UNRESOLVED: {unresolved}/{len(rows)}"
    else:
        result = "OK"
    if unresolved and unresolved != len(rows):
        result += f"; {unresolved} unresolved"
    return result


def finite_min(rows: list[AuditRow], field: str) -> float | None:
    values = [getattr(row, field) for row in rows if getattr(row, field) is not None]
    if not values:
        return None
    return min(values)


def group_rows(rows: list[AuditRow]) -> list[list[AuditRow]]:
    grouped: dict[tuple[object, ...], list[AuditRow]] = defaultdict(list)
    for row in rows:
        key = (
            row.reaction.input_path,
            row.reaction.mol_a,
            row.reaction.iface_a,
            row.reaction.mol_b,
            row.reaction.iface_b,
            row.reaction.sigma,
            row.configured_dt,
            row.count_a,
            row.rho_a,
            row.geometry,
            row.d_trans,
            row.a_a,
            row.a_b,
            row.d_rot,
            row.dt_trans,
            row.dt_rot,
            row.trans_status,
            row.rot_status,
            row.issue,
        )
        grouped[key].append(row)
    return sorted(
        grouped.values(),
        key=lambda group: (
            str(group[0].reaction.input_path),
            group[0].reaction.line,
        ),
    )


def report_link(report_path: Path, target: Path, line: int | None = None) -> str:
    relative = target.relative_to(report_path.parent)
    label = target.relative_to(report_path.parent.parent).as_posix()
    suffix = f"#L{line}" if line is not None else ""
    return f"[{label}{f':{line}' if line is not None else ''}]({relative.as_posix()}{suffix})"


def build_report(root: Path, report_path: Path, rows: list[AuditRow], input_count: int) -> str:
    by_input: dict[Path, list[AuditRow]] = defaultdict(list)
    for row in rows:
        by_input[row.reaction.input_path].append(row)

    trans_failures = sum(row.trans_status == "TOO BIG" for row in rows)
    rot_failures = sum(row.rot_status == "TOO BIG" for row in rows)
    trans_unresolved = sum(row.trans_status == "UNRESOLVED" for row in rows)
    rot_unresolved = sum(row.rot_status == "UNRESOLVED" for row in rows)
    inputs_trans_fail = sum(
        any(row.trans_status == "TOO BIG" for row in input_rows)
        for input_rows in by_input.values()
    )
    inputs_rot_fail = sum(
        any(row.rot_status == "TOO BIG" for row in input_rows)
        for input_rows in by_input.values()
    )

    output: list[str] = [
        "# Sample-input pairwise-reaction timestep audit",
        "",
        "## Result",
        "",
        (
            f"Audited **{len(rows)} pairwise declarations** in **{len(by_input)} input files** "
            f"({input_count} active `.inp` files total). The configured timestep exceeds the "
            f"translational limit for **{trans_failures} declarations in {inputs_trans_fail} "
            f"input files**, and exceeds the rotational-only estimate for **{rot_failures} "
            f"declarations in {inputs_rot_fail} input files**. There are {trans_unresolved} "
            f"translational and {rot_unresolved} rotational unresolved declarations."
        ),
        "",
        "`TOO BIG` means `configured timeStep > calculated Δt_AB`. `OK` includes an infinite "
        "limit when the relevant diffusion estimate is zero.",
        "",
        "## Unresolved declarations",
        "",
        "| Declarations | Input-data issue |",
        "|---:|---|",
    ]
    for issue, count in Counter(row.issue for row in rows if row.issue).most_common():
        output.append(f"| {count} | {issue} |")

    output.extend(
        [
            "",
            "## Method and assumptions",
            "",
            "For every source declaration with two reactant species, this audit uses",
            "",
            r"$$\Delta t_{A,B}=\frac{1}{56D_{AB}}\left[\left(\frac{3}{4\pi\rho_A}+\sigma_{AB}^3\right)^{1/3}-\sigma_{AB}\right]^2,$$",
            "",
            r"first with $D_{AB}=D_A+D_B$, then with $D_{AB}=\tfrac{2}{3}(a_A^2D_{r,A}+a_B^2D_{r,B})$.",
            "",
            "- `rho_A = N_A / V` uses the initial total copy number of the first-listed reactant "
            "and the water-box or sphere volume. Total template copies are a conservative upper "
            "bound for state-conditioned reactions.",
            "- `a` is the Euclidean COM-to-reacting-interface distance read from the local `.mol` file.",
            "- The scalar `D` or `Dr` is the mean of the nonzero vector components. In the current "
            "samples those components are equal, so this is the value supplied on each mobile axis.",
            "- Missing `sigma` on a conditional-rate declaration inherits the matching physical "
            "interface pair's explicit value; otherwise the NERDSS default `sigma = 1 nm` is used.",
            "- The provided 3D equation is applied literally to 3D, 3D-to-2D, and 2D declarations. "
            "NERDSS documents a separate 2D criterion, so surface rows are a requested-equation "
            "screen rather than a complete 2D stability assessment.",
            "- `add.inp` reactions use the timestep, volume, and existing populations from the sibling "
            "`parms.inp`, plus the populations added by `add.inp`.",
            "- Notebook checkpoints and the `#pucadyil#` editor-backup directory are excluded.",
            "",
            "Reference: [NERDSS and Time-Step Selection](https://ionerdss.readthedocs.io/en/latest/nerdss_user_guide_time_step.html).",
            "",
            "## Per-input summary",
            "",
            "| Input | Pairwise declarations | timeStep (µs) | Min translational Δt (µs) | Translational result | Min rotational Δt (µs) | Rotational result | Controlling line(s) |",
            "|---|---:|---:|---:|---|---:|---|---|",
        ]
    )

    for path in sorted(by_input, key=lambda item: item.as_posix()):
        input_rows = by_input[path]
        min_trans = finite_min(input_rows, "dt_trans")
        min_rot = finite_min(input_rows, "dt_rot")
        trans_control = [
            row.reaction.line
            for row in input_rows
            if min_trans is not None and row.dt_trans == min_trans
        ]
        rot_control = [
            row.reaction.line
            for row in input_rows
            if min_rot is not None and row.dt_rot == min_rot
        ]
        controlling = []
        if trans_control:
            controlling.append(f"T: {compress_lines(trans_control)}")
        if rot_control:
            controlling.append(f"R: {compress_lines(rot_control)}")
        output.append(
            "| "
            + " | ".join(
                [
                    report_link(report_path, path),
                    str(len(input_rows)),
                    fmt(input_rows[0].configured_dt),
                    fmt(min_trans),
                    result_summary(input_rows, "trans_status"),
                    fmt(min_rot),
                    result_summary(input_rows, "rot_status"),
                    "; ".join(controlling) if controlling else "—",
                ]
            )
            + " |"
        )

    output.extend(
        [
            "",
            "## Grouped per-reaction audit",
            "",
            "Numerically identical declarations are grouped only within the same input and "
            "reacting interface pair. The declaration count and source lines preserve coverage "
            "of every pairwise reaction.",
            "",
            "| Input lines | Pair | Geometry | n | N_A | rho_A (nm^-3) | sigma (nm) | D_A+D_B (nm^2/µs) | Translational Δt / result (µs) | a_A, a_B (nm) | Estimated rotational D (nm^2/µs) | Rotational Δt / result (µs) | Notes |",
            "|---|---|---|---:|---:|---:|---:|---:|---|---|---:|---|---|",
        ]
    )

    for group in group_rows(rows):
        row = group[0]
        lines = [item.reaction.line for item in group]
        path_cell = report_link(report_path, row.reaction.input_path)
        line_cell = f"{path_cell}<br>lines {compress_lines(lines)}"
        pair = (
            f"`{row.reaction.mol_a}.{row.reaction.iface_a}` + "
            f"`{row.reaction.mol_b}.{row.reaction.iface_b}`"
        )
        trans = f"{fmt(row.dt_trans)} / **{row.trans_status}**"
        rot = f"{fmt(row.dt_rot)} / **{row.rot_status}**"
        sigma_note = "" if row.reaction.sigma_source == "explicit" else row.reaction.sigma_source
        notes = "; ".join(part for part in (sigma_note, row.issue) if part) or "—"
        output.append(
            "| "
            + " | ".join(
                [
                    line_cell,
                    pair,
                    row.geometry,
                    str(len(group)),
                    str(row.count_a) if row.count_a is not None else "—",
                    fmt(row.rho_a),
                    fmt(row.reaction.sigma),
                    fmt(row.d_trans),
                    trans,
                    f"{fmt(row.a_a)}, {fmt(row.a_b)}",
                    fmt(row.d_rot),
                    rot,
                    notes,
                ]
            )
            + " |"
        )

    output.extend(
        [
            "",
            "## Interpretation cautions",
            "",
            "The equation depends on the instantaneous density of the first reactant. This report "
            "uses initial populations because those are the populations available in the inputs. "
            "Rows with `N_A = 0` are unresolved for later simulation times because creation "
            "reactions can raise their density. Likewise, assembly changes complex COM positions "
            "and diffusion constants; the rotational estimate here intentionally uses the "
            "unbound molecule-template geometry and constants requested for this screen.",
            "",
        ]
    )
    return "\n".join(output)


def discover_inputs(sample_root: Path) -> list[Path]:
    return sorted(
        path
        for path in sample_root.rglob("*.inp")
        if ".ipynb_checkpoints" not in path.parts
        and not any(part.startswith("#") and part.endswith("#") for part in path.parts)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Defaults to sample_inputs/TIMESTEP_AUDIT.md under --root",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    sample_root = root / "sample_inputs"
    output_path = (args.output or sample_root / "TIMESTEP_AUDIT.md").resolve()

    inputs = discover_inputs(sample_root)
    rows: list[AuditRow] = []
    for input_path in inputs:
        reactions = parse_reactions(input_path)
        if not reactions:
            continue
        own_text = input_path.read_text(errors="replace")
        context_path = input_path
        texts_for_counts = [own_text]
        if parse_timestep(own_text) is None and input_path.name.lower() == "add.inp":
            sibling = input_path.with_name("parms.inp")
            if sibling.exists():
                context_path = sibling
                texts_for_counts.insert(0, sibling.read_text(errors="replace"))
        context_text = context_path.read_text(errors="replace")
        timestep = parse_timestep(context_text)
        volume = parse_volume(context_text)
        counts = parse_counts(texts_for_counts)
        molecules = parse_molecules(input_path.parent)
        rows.extend(
            audit_reaction(reaction, timestep, volume, counts, molecules)
            for reaction in reactions
        )

    report = build_report(root, output_path, rows, len(inputs))
    output_path.write_text(report)
    print(f"Wrote {output_path}")
    print(f"Active input files: {len(inputs)}")
    print(f"Input files with pairwise reactions: {len({row.reaction.input_path for row in rows})}")
    print(f"Pairwise declarations: {len(rows)}")
    print(f"Translational TOO BIG: {sum(row.trans_status == 'TOO BIG' for row in rows)}")
    print(f"Rotational TOO BIG: {sum(row.rot_status == 'TOO BIG' for row in rows)}")
    print(f"Translational unresolved: {sum(row.trans_status == 'UNRESOLVED' for row in rows)}")
    print(f"Rotational unresolved: {sum(row.rot_status == 'UNRESOLVED' for row in rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
