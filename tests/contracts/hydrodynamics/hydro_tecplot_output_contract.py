#!/usr/bin/env python3
"""Two-step RT serialization contract; no solver tuning or physical accuracy gate."""
import csv
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid


FILES = ("run.dict", "geometry.dict", "mesh.dict",
         "interface.dict", "materials.dict", "hydro.dict")
QUOTED = r'"(?:\\.|[^"\\])*"'


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def number(token):
    token = token.strip()
    if token.lower() in ("true", "false"):
        return float(token.lower() == "true")
    try:
        value = float(token)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def equal_float(actual, expected, label):
    require(math.isfinite(actual) and math.isfinite(expected), label + ": nonfinite")
    require(abs(actual - expected) <= 8 * max(math.ulp(actual), math.ulp(expected)),
            label + f": {actual!r} differs from {expected!r}")


def variables(text):
    match = re.search(r"\bVARIABLES\s*=\s*((?:" + QUOTED + r"\s*,?\s*)+)", text)
    require(match is not None, "missing Tecplot VARIABLES")
    return [json.loads(value) for value in re.findall(QUOTED, match.group(1))]


def attribute(header, name):
    match = re.search(r"\b" + re.escape(name) + r"\s*=\s*(" + QUOTED + r"|[^\s,]+)", header)
    require(match is not None, "missing Tecplot attribute " + name)
    value = match.group(1)
    return json.loads(value) if value.startswith('"') else value


def auxiliary(text, dataset=False):
    prefix = "DATASETAUXDATA" if dataset else "AUXDATA"
    return {name: json.loads(value) for name, value in
            re.findall(r"^\s*" + prefix + r"\s+(\w+)\s*=\s*(" + QUOTED + r")\s*$",
                       text, flags=re.MULTILINE)}


def csv_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        require(reader.fieldnames is not None, "missing CSV header: " + str(path))
        return reader.fieldnames, rows


VOLUME_FIELDS = ["cell_id", "cell_capacity", "phase_volume_fraction", "liquid_volume", "gas_volume",
                 "plic_normal_X", "plic_normal_Y", "plic_normal_Z", "plic_offset", "plic_resolved",
                 "velocity_X", "velocity_Y", "velocity_Z", "pressure", "kappa_HF_force", "kappa_HF_forceAvailable",
                 "kappa_effective_force", "kappa_effective_forceAvailable", "capillary_curvature_method"]


CURVATURE_METHODS = "0=unavailable;1=heightFunction;2=segmentFit;3=underresolvedIsolatedFragmentZero"
CURVATURE_MEANING = "HF field is physical HF/fit estimate; effective field is actual force closure; method3 has physical unavailable, effective0 available"
INTERFACE_FIELDS = ["NX", "NY", "KappaHFForce", "KappaHFForceAvailable",
                    "KappaEffectiveForce", "KappaEffectiveForceAvailable", "CapillaryCurvatureMethod"]


def curvature_metadata(metadata):
    require(metadata.get("CapillaryCurvatureMethods") == CURVATURE_METHODS, "explicit frozen curvature method codes")
    require(metadata.get("CapillaryCurvatureMeaning") == CURVATURE_MEANING, "physical estimate/effective closure meaning is explicit")


def interface_dat(path, unavailable=False):
    text = path.read_text(encoding="utf-8")
    tags = list(re.finditer(r"^\s*AUXDATA\s+\w+\s*=.*$", text, flags=re.MULTILINE))
    require(tags, "interface has native provenance")
    header, tokens = text[:tags[-1].end()], text[tags[-1].end():].split()
    names = variables(header); nodes, elements = int(attribute(header, "NODES")), int(attribute(header, "ELEMENTS"))
    require(names == ["X", "Y"] + INTERFACE_FIELDS, "interface appends effective force and explicit method to original fields")
    require(attribute(header, "ZONETYPE") == "FELINESEG" and attribute(header, "DATAPACKING") == "BLOCK", "native interface line BLOCK")
    require(len(tokens) == len(names) * nodes + 2 * elements, "actual appended interface blocks/connectivity dimensions")
    require(all(number(x) is not None for x in tokens[:len(names)*nodes]), "finite interface field payload")
    arrays = {name: [float(x) for x in tokens[i*nodes:(i+1)*nodes]] for i, name in enumerate(names)}
    connectivity = tokens[len(names)*nodes:]
    require(all(re.fullmatch(r"\d+", x) and 1 <= int(x) <= nodes for x in connectivity), "actual interface native-node connectivity")
    curvature_metadata(auxiliary(header))
    require(all(x in (0, 1, 2, 3) for x in arrays["CapillaryCurvatureMethod"]), "interface method tags are frozen enum values")
    if unavailable:
        require(all(x == 0 for name in INTERFACE_FIELDS[2:] for x in arrays[name]), "sigma-zero/remap interface has no old physical/effective force or method")
    return arrays


def volume_dat(path, frame, conservation, previous_time):
    """Independent DAT geometry/inventory checks; no VTU is read or recreated."""
    text = path.read_text(encoding="utf-8")
    aux_matches = list(re.finditer(r"^\s*AUXDATA\s+\w+\s*=.*$", text, flags=re.MULTILINE))
    require(aux_matches, "volume DAT must declare accepted-state provenance")
    end = aux_matches[-1].end()
    header, tokens = text[:end], text[end:].split()
    require(len(re.findall(r"^\s*ZONE\b", header, flags=re.MULTILINE)) == 1, "one native volume zone")
    require(attribute(header, "ZONETYPE") == "FEQUADRILATERAL" and attribute(header, "DATAPACKING") == "BLOCK", "native FEQUAD BLOCK")
    nodes, cells = int(attribute(header, "NODES")), int(attribute(header, "ELEMENTS"))
    require(nodes > 0 and cells > 8 * 32, "test retains native interface-band AMR")
    names = variables(header)
    require(names == ["X", "Y"] + VOLUME_FIELDS, "exact native appended cell-centered fields")
    cc = re.search(r"\bVARLOCATION\s*=\s*\(\[([\s\d,]+)\]\s*=\s*CELLCENTERED\s*\)", header)
    require(cc is not None and [int(value) for value in re.findall(r"\d+", cc.group(1))] == list(range(3, len(names) + 1)), "only X/Y are nodal")
    types = re.search(r"\bDT\s*=\s*\((.*?)\)", header, flags=re.DOTALL)
    require(types is not None and types.group(1).split() == ["DOUBLE"] * len(names), "all explicit DOUBLE variables")
    count = 2 * nodes + len(VOLUME_FIELDS) * cells
    require(len(tokens) == count + 4 * cells, "actual blocks and quad connectivity dimensions")
    require(all(number(token) is not None for token in tokens[:count]), "every numeric volume payload is finite")
    arrays, offset = {}, 0
    for name in names:
        size = nodes if name in ("X", "Y") else cells
        arrays[name] = [float(token) for token in tokens[offset:offset + size]]
        offset += size
    require(all(re.fullmatch(r"\d+", token) for token in tokens[count:]), "integer one-based connectivity")
    connection = [int(token) - 1 for token in tokens[count:]]
    require(all(0 <= value < nodes for value in connection), "actual quad node bounds")
    identifiers = arrays["cell_id"]
    require(len(set(identifiers)) == cells and all(value.is_integer() and value >= 0 for value in identifiers), "unique native cell identities")
    require(all(value in (0, 1) for value in arrays["plic_resolved"]), "actual PLIC resolution flags")
    rectangles, widths = [], set()
    for index in range(cells):
        vertices = connection[4 * index:4 * index + 4]
        require(len(set(vertices)) == 4, "quad has four distinct native nodes")
        polygon = [(arrays["X"][v], arrays["Y"][v]) for v in vertices]
        xs, ys = sorted(set(x for x, _ in polygon)), sorted(set(y for _, y in polygon))
        require(len(xs) == len(ys) == 2 and set(polygon) == {(x, y) for x in xs for y in ys}, "actual Cartesian leaf rectangle")
        require(0 <= xs[0] < xs[1] <= 1 and 0 <= ys[0] < ys[1] <= 4, "native cells remain inside original RT domain")
        dx, dy = xs[1] - xs[0], ys[1] - ys[0]
        widths.update((dx, dy))
        signed = .5 * math.fsum(polygon[j][0] * polygon[(j + 1) % 4][1] - polygon[(j + 1) % 4][0] * polygon[j][1] for j in range(4))
        require(signed != 0, "nondegenerate ordered quad")
        capacity, liquid, gas, alpha = (arrays[name][index] for name in ("cell_capacity", "liquid_volume", "gas_volume", "phase_volume_fraction"))
        require(capacity > 0 and liquid >= 0 and gas >= 0 and 0 <= alpha <= 1, "retained independent phase stocks")
        equal_float(abs(signed), dx * dy, "quad area from real connectivity")
        equal_float(capacity, dx * dy, "real leaf area equals capacity")
        equal_float(liquid + gas, capacity, "independent phase closure")
        equal_float(alpha * capacity, liquid, "alpha and liquid stock")
        rectangles.append((xs[0], xs[1], ys[0], ys[1]))
    equal_float(math.fsum(arrays["cell_capacity"]), 4., "original RT domain measure")
    require(len(widths) >= 2, "actual coarse and fine cell widths")
    equal_float(max(widths) / min(widths), 2., "true static L1 2:1 grid")
    for index, a in enumerate(rectangles):
        for b in rectangles[index + 1:]:
            require(min(a[1], b[1]) <= max(a[0], b[0]) or min(a[3], b[3]) <= max(a[2], b[2]), "leaf rectangle interiors cannot overlap")
    total_liquid, total_gas = math.fsum(arrays["liquid_volume"]), math.fsum(arrays["gas_volume"])
    equal_float(total_liquid, float(conservation["negative_phase_volume"]), "DAT stocks and current conservation total")
    equal_float(3 * total_liquid + total_gas, float(conservation["total_geometric_mass"]), "original RT phase densities and actual mass")
    metadata = auxiliary(header)
    curvature_metadata(metadata)
    require(int(metadata["AcceptedStep"]) == int(frame["step"]) and int(metadata["GeometryIdentity"]) == int(frame["geometry_identity"]), "accepted DAT provenance")
    clocks = {"VTK_TimeValue", "VTK_accepted_step", "VTK_geometry_identity", "VTK_pressure_sample_time", "VTK_hf_force_time", "VTK_hf_force_geometry_identity"}
    require(clocks <= metadata.keys(), "all prior native FieldData clocks retained as VTK AUX")
    require(int(metadata["VTK_accepted_step"]) == int(frame["step"]) and int(metadata["VTK_geometry_identity"]) == int(frame["geometry_identity"]), "native accepted clock and identity")
    actual_time = float(conservation["time"])
    equal_float(float(frame["time"]), actual_time, "frame equals actual conservation clock")
    equal_float(float(attribute(header, "SOLUTIONTIME")), actual_time, "DAT solution time")
    equal_float(float(metadata["VTK_TimeValue"]), actual_time, "native accepted time")
    expected_pressure_time = 0. if int(frame["step"]) == 0 else .5 * (previous_time + actual_time)
    equal_float(float(metadata["VTK_pressure_sample_time"]), expected_pressure_time, "real midpoint pressure time")
    require(metadata["VTK_hf_force_time"].lower() == "nan" and int(metadata["VTK_hf_force_geometry_identity"]) == 0, "zero sigma has no invented HF source clock")
    require(all(value == 0 for value in arrays["kappa_HF_force"]) and all(value == 0 for value in arrays["kappa_HF_forceAvailable"]), "zero sigma curvature explicitly unavailable")
    require(all(value == 0 for name in ("kappa_effective_force", "kappa_effective_forceAvailable", "capillary_curvature_method") for value in arrays[name]), "zero sigma has no invented effective curvature or closure method")
    explanation = " ".join(metadata.values()).lower()
    require("placeholder" in explanation and "available" in explanation, "unavailable curvature placeholder disclosed")
    return arrays["X"], arrays["Y"], connection, arrays["cell_id"], arrays["cell_capacity"]


def load_conservation(output):
    return [decoded_row(row) for row in table_history(output / "series/conservation.dat", "conservation")]


def check_volume_output(output, expected_steps, completed=True):
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest["status"] == ("COMPLETED" if completed else "FAILED"), "actual completion/failure status")
    # Candidate A must declare its public DAT-only volume format explicitly.
    require(manifest.get("volume_format") == "tecplotDat" and manifest.get("frame_index") == "meta/frames.csv", "manifest declares native DAT volume format")
    require(manifest.get("capillary_method") == "alphaHeightFunction", "manifest retains the existing physical capillary method")
    require("capillary_curvature=" in (output / "meta/build_info.txt").read_text(encoding="utf-8"), "build metadata declares physical/effective curvature semantics")
    require(not list(output.rglob("*.vtu")) and not list(output.rglob("*.pvd")) and not (output / "meta/tecplot_frames.csv").exists(), "VOF emits no VTU/PVD/duplicate Tecplot index")
    fields, frames = csv_rows(output / "meta/frames.csv")
    require(fields == ["frame", "step", "time", "geometry_identity", "volume", "interface"], "one authoritative frames schema")
    require([int(row["step"]) for row in frames] == list(expected_steps), "initial and every accepted frame")
    require(manifest["complete_frames"] == len(frames), "manifest real frame count")
    require(all(path.is_file() for path in (output / "fields").iterdir()), "accepted fields are flat without category directories")
    for frame in frames:
        for key in ("volume", "interface"):
            require(not frame[key] or Path(frame[key]).parent == Path("fields"), "frame index uses flat fields paths")
    conservation = {int(row["step"]): row for row in load_conservation(output)}
    mesh = None
    for frame in frames:
        step = int(frame["step"])
        require(frame["volume"].endswith(".dat"), "authoritative volume path points directly to DAT")
        require(step in conservation, "every frame has actual accepted conservation row")
        previous_time = 0. if step == 0 else float(conservation[step - 1]["time"])
        current = volume_dat(output / frame["volume"], frame, conservation[step], previous_time)
        if frame["interface"]:
            interface_dat(output / frame["interface"], unavailable=True)
        require(mesh is None or mesh == current, "static AMR actual nodes, connectivity and capacities remain fixed")
        mesh = current
    require(int(frames[-1]["step"]) == manifest["accepted_step"] and int(frames[-1]["geometry_identity"]) == manifest["geometry_identity"], "last physical frame is the accepted prefix")
    equal_float(float(frames[-1]["time"]), float(manifest["accepted_time"]), "manifest accepted prefix time")
    return manifest, frames

SERIES_SCHEMAS = {'conservation': ['step', 'time', 'dt', 'geometry_identity', 'negative_phase_volume', 'initial_volume', 'relative_volume_drift', 'total_geometric_mass', 'divergence_max', 'dual_R_max', 'native_max_speed', 'kappa_max', 'phase_volume_closure', 'mass_closure', 'liquid_p0_moment_x', 'liquid_p0_moment_y', 'liquid_p0_moment_xx', 'liquid_p0_moment_yy', 'liquid_p0_moment_xy', 'vof_work_kinetic_change', 'vof_work_capillary', 'vof_work_body', 'vof_work_viscous', 'vof_work_pressure', 'vof_work_half_R', 'vof_work_half_R_mapped_primal', 'vof_work_half_R_geometry_transfer', 'vof_work_half_R_spatial_exchange', 'vof_work_time_mass', 'vof_work_rounded_phase_mass', 'vof_work_operator_rate_roundoff', 'vof_work_raw_minus_canonical_outflow_half', 'vof_work_half_R_decomposition_residual', 'vof_work_momentum_equation', 'vof_work_kinetic_budget_residual', 'vof_work_closure_residual', 'vof_work_closure_relative_residual', 'vof_stable_gas_duals', 'vof_stable_gas_half_R_work', 'vof_stable_gas_spatial_exchange_half_R_work', 'vof_stable_gas_max_abs_integrated_R_over_mid_mass', 'vof_peak_gas_native_face', 'vof_peak_gas_abs_end_normal_velocity', 'vof_peak_gas_integrated_R_over_mid_mass', 'vof_peak_gas_mass_change_over_mid_mass', 'vof_peak_gas_outward_integral_over_mid_mass', 'vof_peak_gas_geometry_transfer_over_mid_mass', 'vof_peak_gas_spatial_exchange_over_mid_mass'], 'energy': ['step', 'time', 'kinetic', 'surface', 'total', 'body_power', 'cumulative_body_work', 'cumulative_viscous_loss', 'work_corrected_total', 'capillary_energy_defect', 'interface_power', 'interface_power_available', 'capillary_power', 'capillary_power_available'], 'solver': ['step', 'time', 'stage', 'status', 'AC_stages', 'linear_iterations', 'linear_refinements', 'algebraic_residual', 'momentum_residual', 'reason'], 'performance': ['step', 'time', 'stage', 'wall_seconds', 'cells', 'native_velocity_dofs', 'pressure_dofs'], 'vof_native_velocity': ['step', 'accepted_time', 'geometry_identity', 'force_time', 'force_geometry_identity', 'global_abs_normal_velocity', 'global_normal_velocity', 'global_velocity_dof', 'global_native_face', 'global_normal_axis', 'global_face_x', 'global_face_y', 'global_negative_cell', 'global_positive_cell', 'global_negative_alpha', 'global_positive_alpha', 'global_midpoint_capillary_load', 'global_midpoint_force_kappa', 'pure_liquid_abs_normal_velocity', 'pure_liquid_normal_velocity', 'pure_liquid_velocity_dof', 'pure_liquid_native_face', 'pure_liquid_normal_axis', 'pure_liquid_face_x', 'pure_liquid_face_y', 'pure_liquid_negative_cell', 'pure_liquid_positive_cell', 'pure_liquid_negative_alpha', 'pure_liquid_positive_alpha', 'pure_liquid_midpoint_capillary_load', 'pure_liquid_midpoint_force_kappa', 'pure_gas_abs_normal_velocity', 'pure_gas_normal_velocity', 'pure_gas_velocity_dof', 'pure_gas_native_face', 'pure_gas_normal_axis', 'pure_gas_face_x', 'pure_gas_face_y', 'pure_gas_negative_cell', 'pure_gas_positive_cell', 'pure_gas_negative_alpha', 'pure_gas_positive_alpha', 'pure_gas_midpoint_capillary_load', 'pure_gas_midpoint_force_kappa', 'interface_abs_normal_velocity', 'interface_normal_velocity', 'interface_velocity_dof', 'interface_native_face', 'interface_normal_axis', 'interface_face_x', 'interface_face_y', 'interface_negative_cell', 'interface_positive_cell', 'interface_negative_alpha', 'interface_positive_alpha', 'interface_midpoint_capillary_load', 'interface_midpoint_force_kappa', 'geometry_pure_gas_abs_normal_velocity', 'geometry_pure_gas_normal_velocity', 'geometry_pure_gas_velocity_dof', 'geometry_pure_gas_native_face', 'geometry_pure_gas_normal_axis', 'geometry_pure_gas_face_x', 'geometry_pure_gas_face_y', 'geometry_pure_gas_negative_cell', 'geometry_pure_gas_positive_cell', 'geometry_pure_gas_negative_alpha', 'geometry_pure_gas_positive_alpha', 'geometry_pure_gas_midpoint_capillary_load', 'geometry_pure_gas_midpoint_force_kappa'], 'vof_csf': ['step', 'accepted_time', 'force_time', 'force_geometry_identity', 'force_kappa_min', 'force_kappa_max', 'loaded_dofs', 'height_cells', 'height_columns', 'height_column_cell_queries', 'height_alternate_directions', 'height_extended_cells', 'height_extended_slices', 'height_refined_cells', 'height_refined_slices', 'height_fit_cells', 'height_fit_segments'], 'events': ['event', 'attempted_step', 'accepted_step', 'accepted_time', 'old_dt', 'new_dt', 'stage', 'reason']}

SERIES_SCHEMAS["vof_csf"].append("unresolved_isolated_fragment_cells")

SERIES_SCHEMAS["time_step"] = ["step", "time", "begin", "actual_dt", "proposed_dt", "flow_rate",
    "capillary_rate", "interface_h", "flow_dt_bound", "capillary_dt_bound", "growth_dt_bound",
    "max_dt", "reason", "retries", "endpoint_aligned", "final_below_min", "minimum_clock_rounded"]

SERIES_SCHEMAS['regrid_decision'] = ['step', 'time', 'begin', 'end', 'same_end', 'source_geometry_identity', 'target_geometry_identity', 'source_leaves', 'target_leaves', 'mixed_cells', 'geometric_mixed_cells', 'unresolved_mixed_cells', 'pure_jump_faces', 'interface_seed_cells', 'interface_band_cells', 'flow_refine_cells', 'flow_coarsen_cells', 'maximum_h_abs_vorticity', 'additional_refine_cells', 'additional_hold_cells', 'requested', 'changed', 'clones', 'sweeps', 'indicator_seconds', 'planning_seconds', 'adaptation_seconds', 'decision_seconds', 'inventory_transfer_seconds', 'regrid_seconds']
SERIES_SCHEMAS['regrid'] = ['step', 'time', 'begin', 'end', 'same_end', 'source_geometry_identity', 'target_geometry_identity', 'source_leaves', 'target_leaves', 'refined_donors', 'refined_targets', 'coarsened_targets', 'unchanged_targets', 'resolved_refined_donors', 'unresolved_refined_donors', 'p0_refined_pieces', 'corrected_refined_donors', 'roundoff_closed_phases', 'source_liquid', 'target_liquid', 'source_gas', 'target_gas', 'donor_relative_liquid', 'donor_relative_gas', 'global_relative_liquid', 'global_relative_gas', 'raw_inventory_relative', 'raw_geometry_relative', 'target_phase_closure', 'history_reset', 'remap_solves', 'source_native_P_x', 'target_native_P_x', 'source_native_P_y', 'target_native_P_y', 'delta_native_P_x', 'delta_native_P_y', 'source_native_kinetic', 'target_native_kinetic', 'delta_native_kinetic', 'projection_delta_P_x', 'projection_delta_P_y', 'projection_delta_kinetic', 'projection_correction_energy', 'projection_energy_identity_residual', 'raw_relative_divergence', 'final_relative_divergence', 'raw_divergence_L1', 'final_divergence_L1', 'pressure_gauge_shift', 'pressure_gauge_relative_residual', 'source_leaf_candidates', 'face_spans', 'polynomial_integrals', 'wall_normals', 'axis_normals', 'decision_seconds', 'inventory_transfer_seconds', 'raw_flux_seconds', 'projection_seconds', 'fluid_remap_seconds', 'regrid_seconds']

CATEGORY_COLUMNS = {"stage", "status", "reason", "event"}
PERFORMANCE_STAGES = {"initialization", "io", "step_total", "VOF_AC", "vof_transport", "vof_material",
                      "vof_hf", "vof_frozen_ns", "vof_post_diagnostics", "vof_ac_stages", "vof_other_and_retries",
                      "vof_maxwell", "vof_phase_geometry", "vof_material_transport", "vof_dual_mass",
                      "vof_convection", "vof_viscosity"}


def table_history(path, logical):
    """Read a native wide POINT table, retaining raw tokens and reversible text."""
    fields = SERIES_SCHEMAS[logical]
    numeric_columns = [name for name in fields if name not in CATEGORY_COLUMNS]
    categories = [name for name in fields if name in CATEGORY_COLUMNS]
    text = path.read_text(encoding="utf-8")
    pieces = re.split(r"(?m)(?=^ZONE\b)", text)
    require(len(pieces) > 1, "native table has at least one real zone: " + str(path))
    require(variables(pieces[0]) == numeric_columns + [name + "Available" for name in numeric_columns], "native table preserves numeric schema/order and explicit availability: " + path.name)
    dataset = auxiliary(pieces[0], dataset=True)
    require(dataset.get("Generator") == "LAMP_TableHistory_v1", "native wide TableHistory generator")
    for index, name in enumerate(fields):
        require(dataset.get(f"Column{index}Name") == name and
                dataset.get(f"Column{index}Kind") == ("Category" if name in CATEGORY_COLUMNS else "Numeric"),
                "dataset preserves original full schema order/kind")
    clock_column = next((name for name in ("time", "accepted_time", "step", "accepted_step", "timestep", "attempted_step") if name in numeric_columns), "none")
    require(dataset.get("ClockColumn") == clock_column, "dataset declares the actual original clock column")
    require("nonincreasing" in dataset.get("ZonePolicy", "").lower(), "dataset declares duplicate/rollback clock boundaries")
    require("hex" in dataset.get("CategoryEncoding", "").lower(), "dataset declares reversible category bytes")
    policy = " ".join(dataset.values()).lower()
    require("available" in policy and "0" in policy and any(word in policy for word in ("blank", "undefined", "unavailable", "missing")), "missing-value mask/blanking rule is explicit")
    all_rows = []
    for zone_index, zone in enumerate(pieces[1:]):
        match = re.search(r"(?m)^\s*[+\-]?(?:\d|\.\d)", zone)
        require(match is not None, "native table zone has real data")
        header = zone[:match.start()]
        declaration = header.splitlines()[0]
        require(re.search(r"\bI\s*=", declaration) is None, "appendable native stream deliberately omits I")
        require(attribute(header, "DATAPACKING") == "POINT", "native wide table uses POINT")
        types = re.search(r"\bDT\s*=\s*\((.*?)\)", header, flags=re.DOTALL)
        require(types is not None and types.group(1).split() == ["DOUBLE"] * (2 * len(numeric_columns)), "native table explicit DOUBLE columns")
        metadata = auxiliary(header)
        classified = {}
        for index, name in enumerate(categories):
            prefix = f"Category{index}"
            require(metadata.get(prefix + "Name") == name, "classification original name/order")
            raw_hex = metadata.get(prefix + "ValueHex")
            require(raw_hex is not None and re.fullmatch(r"(?:[0-9a-f]{2})*", raw_hex) is not None, "classification bytes are reversible lowercase hex")
            raw = bytes.fromhex(raw_hex).decode("utf-8")
            require(prefix + "Value" in metadata, "classification has readable display")
            if not any(ord(c) < 32 or ord(c) == 127 for c in raw):
                require(metadata[prefix + "Value"] == raw, "readable category text retains original characters")
            classified[name] = raw
        require(len([key for key in metadata if re.fullmatch(r"Category\d+Name", key)]) == len(categories), "no synthetic category columns")
        previous = None
        data_rows = [line.split() for line in zone[match.start():].splitlines() if line.strip()]
        require(data_rows, "zone cannot be a synthetic empty history")
        for tokens in data_rows:
            require(len(tokens) == 2 * len(numeric_columns) and all(number(token) is not None for token in tokens), "every complete native wide row is finite")
            size = len(numeric_columns)
            payload, mask = tokens[:size], tokens[size:]
            require(all(token in ("0", "1") for token in mask), "availability is explicit Boolean numeric data")
            values = dict(zip(numeric_columns, payload))
            valid = dict(zip(numeric_columns, map(int, mask)))
            require(all(valid[name] or float(values[name]) == 0 for name in numeric_columns), "missing payload only zero with Available=0")
            clock_name = next((name for name in ("time", "accepted_time", "step", "attempted_step") if name in values), None)
            if clock_name and valid[clock_name]:
                clock = float(values[clock_name])
                require(previous is None or clock > previous, "one zone is a real strictly increasing clock segment")
                previous = clock
            all_rows.append({"tokens": values, "available": valid, "category": classified.copy(), "zone": zone_index})
    return all_rows


def decoded_row(row):
    return {name: token if row["available"][name] else "NA" for name, token in row["tokens"].items()}


def physical_series_only(output):
    """The default keeps physical histories and essential failure/solver records."""
    names = {path.stem for path in (output / "series").glob("*.dat")}
    allowed = {"conservation", "energy", "solver", "vof_csf", "electric", "maxwell_work",
               "charge_balance", "phase_boundary", "vof_boundary_work", "regrid_projection_work",
               "time_step", "regrid", "regrid_decision", "events"}
    require(names <= allowed, "physical profile emits no developer histories: " + str(names - allowed))
    require({"conservation", "energy", "solver"} <= names, "physical process and essential solver histories remain")
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("series_profile") == "physical", "default series profile is disclosed as physical")


def native_series(output, expected_steps, diagnostic=False):
    directory = output / "series"
    require(not list(directory.rglob("*.csv")), "current VOF series writes no CSV")
    require(not (directory / "performance.dat").exists(), "no duplicate unsplit performance history")
    require(not (directory / "time_step.dat").exists(), "fixed mode creates no adaptive history")
    require(not (directory / "events.dat").exists() and not (directory / "vof_csf.dat").exists(), "no-data events and sigma-zero HF histories are lazy")
    tables = {}
    required = ("conservation", "energy", "solver") + (("vof_native_velocity",) if diagnostic else ())
    for name in required:
        path = directory / (name + ".dat")
        require(path.is_file(), "native primary table exists: " + name)
        rows = table_history(path, name)
        require([int(float(row["tokens"]["step"])) for row in rows] == list(expected_steps), "actual initial and accepted table rows: " + name)
        require(all(row["available"]["step"] == 1 for row in rows), "real step availability")
        tables[name] = rows
    for row in tables["energy"]:
        require(all(row["available"][name] == 0 and float(row["tokens"][name]) == 0 for name in ("surface", "total", "work_corrected_total", "capillary_energy_defect")), "unavailable VOF area/energy is explicitly masked")
        require(all(row["available"][name] == 1 for name in ("kinetic", "body_power", "cumulative_body_work", "cumulative_viscous_loss")), "actual available energy diagnostics")
        require(float(row["tokens"]["interface_power_available"]) == 0 and row["available"]["interface_power_available"] == 1, "real false flag is distinct from undefined quantity")
    if diagnostic:
        for conservation, velocity in zip(tables["conservation"], tables["vof_native_velocity"]):
            equal_float(float(velocity["tokens"]["global_abs_normal_velocity"]), float(conservation["tokens"]["native_max_speed"]), "same accepted native velocity diagnostic")
        require((directory / "reconstruction_templates.dat").is_file(), "diagnostic profile retains reconstruction bank counters")
    else:
        physical_series_only(output)
    paths = {path.stem.removeprefix("performance_"): path for path in directory.glob("performance_*.dat")}
    require(set(paths) == (PERFORMANCE_STAGES if diagnostic else set()), "profile selects exact fixed provider performance channels")
    for stage, path in paths.items():
        rows = table_history(path, "performance")
        target_steps = [0] if stage == "initialization" else list(expected_steps) if stage in ("io", "step_total") else list(expected_steps)[1:]
        require([int(float(row["tokens"]["step"])) for row in rows] == target_steps, "stage keeps its real accepted samples: " + stage)
        require(all(row["category"]["stage"] == stage for row in rows), "performance filename and exact stage bytes")
        if 1 in target_steps and 2 in target_steps:
            require(len({row["zone"] for row in rows if int(float(row["tokens"]["step"])) in (1, 2)}) == 1, "stage steps 1/2 form a real multi-point curve")
    require(len({row["zone"] for row in tables["conservation"]}) == 1 and len(tables["conservation"]) == 3, "accepted history is not one singleton zone per step")
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("series_format") == "tecplotDat", "manifest declares native DAT-only series")
    require(manifest.get("series_profile") == ("diagnostic" if diagnostic else "physical"), "manifest declares the selected series profile")
    return sum(len(rows) for rows in tables.values())


def select_diagnostic_series(text):
    require(len(re.findall(r"\bOutput\s*\{", text)) == 1, "one copied output block selects diagnostic series")
    series = list(re.finditer(r"\bSeries\s*\{([^{}]*)\}", text))
    require(len(series) <= 1 and len(series) == len(re.findall(r"\bSeries\s*\{", text)),
            "at most one complete copied series block")
    if series:
        block = series[0]
        body, count = re.subn(r"\bprofile\s+(?:physical|diagnostic)\s*;", "profile diagnostic;", block[1])
        require(count == 1 and len(re.findall(r"\bprofile\b", block[1])) == 1,
                "one copied series profile selects diagnostic histories")
        return text[:block.start(1)] + body + text[block.end(1):]
    text, count = re.subn(r"(\bOutput\s*\{)", r"\1\n    Series { profile diagnostic; }", text)
    require(count == 1, "one copied output block adds diagnostic series")
    return text


def diagnostic_profile(executable, positive_task, environment):
    """Changing output verbosity preserves every accepted field and physical row."""
    task = positive_task.parent / (positive_task.name + "_diagnostic")
    shutil.copytree(positive_task / "system", task / "system")
    path = task / "system/run.dict"
    text = select_diagnostic_series(path.read_text(encoding="utf-8"))
    path.write_text(text, encoding="utf-8")
    process = subprocess.run([str(executable), "--threads", "2"], cwd=task, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", timeout=120)
    (task / "contract.log").write_text(process.stdout, encoding="utf-8")
    require(process.returncode == 0 and "[COMPLETED]" in process.stdout,
            "diagnostic copied RT status; see " + str(task / "contract.log"))
    output, physical = task / "output", positive_task / "output"
    check_volume_output(output, (0, 1, 2))
    native_series(output, (0, 1, 2), diagnostic=True)
    for name in ("conservation", "energy", "solver"):
        require(table_history(output / ("series/" + name + ".dat"), name) ==
                table_history(physical / ("series/" + name + ".dat"), name),
                "profile leaves physical history tokens, availability and categories identical: " + name)
    for relative in ("meta/frames.csv",):
        require((output / relative).read_bytes() == (physical / relative).read_bytes(),
                "profile preserves accepted field index")
    fields = {path.name: path.read_bytes() for path in (physical / "fields").iterdir()}
    require(fields == {path.name: path.read_bytes() for path in (output / "fields").iterdir()},
            "profile preserves every accepted field byte")
    return task


def failure_prefix(executable, positive_task, environment):
    """Test-only copied inputs force existing typed-CFL retry exhaustion."""
    task = positive_task.parent / (positive_task.name + "_courant_rejection")
    shutil.copytree(positive_task / "system", task / "system")
    path = task / "system/hydro.dict"
    text = path.read_text(encoding="utf-8")
    if re.search(r"\bmaximumLocalCourant\s+[^;]+;", text):
        text, count = re.subn(r"\bmaximumLocalCourant\s+[^;]+;", "maximumLocalCourant 1e-12;", text)
        require(count == 1, "one copied CFL control")
    else:
        text, count = re.subn(r"(\bHydrodynamics\s*\{)", r"\1\n    Numerics { maximumLocalCourant 1e-12; }", text)
        require(count == 1 and "Numerics" not in path.read_text(encoding="utf-8"), "insert copied CFL-only Numerics block")
    path.write_text(text, encoding="utf-8")
    process = subprocess.run([str(executable), "--threads", "2"], cwd=task, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", timeout=120)
    (task / "contract.log").write_text(process.stdout, encoding="utf-8")
    require(process.returncode != 0, "negative copied CFL control must fail")
    output = task / "output"
    failure = json.loads((output / "meta/failure.json").read_text(encoding="utf-8"))
    require(failure["accepted_step"] < 2 and failure["accepted_time"] < .004 and failure["accepted_state_preserved"], "failure retains a real accepted prefix")
    reason = failure["reason"]
    require("Courant" in reason, "negative failure must be the requested local CFL condition")
    check_volume_output(output, range(failure["accepted_step"] + 1), completed=False)
    require(not list((output / "series").glob("*.csv")), "failure exit remains DAT-only")
    events = table_history(output / "series/events.dat", "events")
    require(any(row["category"]["event"] == "retry" for row in events), "typed CFL enters existing StepRetry")
    final = [row for row in events if row["category"]["event"] == "failure"]
    require(final and any(row["category"]["reason"] == reason for row in final), "failure reason original bytes survive DAT CategoryValueHex")
    require(all(row["available"]["accepted_time"] == 1 and float(row["tokens"]["accepted_time"]) <= failure["accepted_time"] for row in events), "rejected attempts do not advance accepted time")
    physical_series_only(output)
    return task

def adaptive_history(executable, positive_task, environment, reject=False):
    """Actual accepted-clock serialization, including a later policy failure."""
    suffix = "_adaptive_policy_rejection" if reject else "_adaptive"
    task = positive_task.parent / (positive_task.name + suffix)
    shutil.copytree(positive_task / "system", task / "system")
    path = task / "system/run.dict"
    text = path.read_text(encoding="utf-8")
    text, count = re.subn(r"\bendTime\s+[^;]+;", "endTime .0067;", text)
    require(count == 1, "one adaptive test-only endpoint")
    minimum, target = (".002", "1e-12") if reject else (".000125", ".25")
    text, count = re.subn(r"\bminDt\s+[^;]+;", "minDt " + minimum + ";", text)
    require(count == 1, "copied retry minDt shares adaptive minimum")
    text, count = re.subn(r"(\bTime\s*\{)", lambda match: match[1] +
        "\n    Adaptive { minDt " + minimum + "; maxDt .004; targetCourant " + target +
        "; safety .9; maxGrowth 1.25; }", text)
    require(count == 1, "explicit test-only Adaptive block")
    path.write_text(text, encoding="utf-8")
    process = subprocess.run([str(executable), "--threads", "2"], cwd=task, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", timeout=120)
    (task / "contract.log").write_text(process.stdout, encoding="utf-8")
    require((process.returncode != 0) if reject else
            (process.returncode == 0 and "[COMPLETED]" in process.stdout),
            "adaptive copied RT status; see " + str(task / "contract.log"))
    output = task / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    accepted = manifest["accepted_step"]
    if reject:
        failure = json.loads((output / "meta/failure.json").read_text(encoding="utf-8"))
        require(accepted == 1 and failure["accepted_step"] == accepted and
                failure["accepted_state_preserved"], "policy failure preserves one true accepted sample")
        require("adaptive" in failure["reason"].lower() and "minDt" in failure["reason"],
                "failure is the requested adaptive lower bound")
    else:
        require(accepted >= 3, "short adaptive run records growth and final samples")
        equal_float(float(manifest["accepted_time"]), .0067, "adaptive true final time")
    check_volume_output(output, range(accepted + 1), completed=not reject)
    rows = table_history(output / "series/time_step.dat", "time_step")
    conservation = load_conservation(output)
    require(len(rows) == accepted and len(conservation) == accepted + 1,
            "one adaptive row per committed step, no initial or rejected row")
    events_path = output / "series/events.dat"
    events = table_history(events_path, "events") if events_path.exists() else []
    retry_counts = Counter(int(float(row["tokens"]["attempted_step"])) for row in events
                           if row["category"]["event"] == "retry")
    reasons = {"maxDt", "initialDt", "growth", "flow", "capillary", "endTime"}
    for index, (row, state) in enumerate(zip(rows, conservation[1:]), 1):
        values, valid = row["tokens"], row["available"]
        actual_time, previous_time = float(state["time"]), float(conservation[index - 1]["time"])
        require(int(float(values["step"])) == index, "adaptive actual accepted step")
        equal_float(float(values["time"]), actual_time, "adaptive accepted time from conservation")
        equal_float(float(values["begin"]), previous_time, "adaptive actual begin")
        actual_dt = actual_time - previous_time
        equal_float(float(values["actual_dt"]), actual_dt, "adaptive actual interval")
        equal_float(float(state["dt"]), actual_dt, "conservation and adaptive dt agree")
        require(valid["proposed_dt"] == 1 and float(values["proposed_dt"]) > 0,
                "true positive original proposal")
        require(float(values["actual_dt"]) <= float(values["proposed_dt"]) +
                8 * max(math.ulp(actual_dt), math.ulp(float(values["proposed_dt"]))),
                "actual interval cannot hide a reduced retry")
        require(row["category"]["reason"] in reasons, "exact proposal limiter category bytes")
        require(valid["flow_rate"] == 1 and float(values["flow_rate"]) >= 0,
                "measured zero flow is available")
        for name in ("capillary_rate", "interface_h", "capillary_dt_bound"):
            require(valid[name] == 0 and float(values[name]) == 0,
                    "sigma zero is an absent constraint: " + name)
        if float(values["flow_rate"]) == 0:
            require(valid["flow_dt_bound"] == 0 and float(values["flow_dt_bound"]) == 0,
                    "zero flow has no fabricated finite bound")
        else:
            require(valid["flow_dt_bound"] == 1, "finite actual flow bound is available")
            equal_float(float(values["flow_dt_bound"]), .9 * float(target) /
                        float(values["flow_rate"]), "flow bound uses report values")
        equal_float(float(values["max_dt"]), .004, "actual configured max bound")
        expected_growth = .002 if index == 1 else float(conservation[index - 1]["dt"]) * 1.25
        equal_float(float(values["growth_dt_bound"]), expected_growth, "original proposal growth bound")
        require(int(float(values["retries"])) == retry_counts[index], "actual retry count")
        require(all(valid[name] == 1 for name in ("step", "time", "begin", "actual_dt",
                    "proposed_dt", "growth_dt_bound", "max_dt", "retries", "endpoint_aligned", "final_below_min", "minimum_clock_rounded")),
                "clock, actual zero flags and finite controls remain available")
        require(float(values["endpoint_aligned"]) in (0., 1.) and
                float(values["minimum_clock_rounded"]) in (0., 1.) and
                float(values["final_below_min"]) == float(actual_time == .0067 and actual_dt < float(minimum)),
                "endpoint/minimum flags describe the actual accepted clock")
    require(not list((output / "series").glob("*.csv")), "adaptive native series remains DAT only")
    physical_series_only(output)
    return task


def dynamic_payload(path):
    text = path.read_text(encoding="utf-8")
    tags = list(re.finditer(r"^\s*AUXDATA\s+\w+\s*=.*$", text, flags=re.MULTILINE))
    require(tags, "dynamic field provenance")
    header, tokens = text[:tags[-1].end()], text[tags[-1].end():].split()
    names = variables(header)
    nodes, cells = int(attribute(header, "NODES")), int(attribute(header, "ELEMENTS"))
    require(names == ["X", "Y"] + VOLUME_FIELDS and nodes > 0 and cells > 0, "dynamic native volume schema")
    count = 2 * nodes + len(VOLUME_FIELDS) * cells
    require(len(tokens) == count + 4 * cells and all(number(x) is not None for x in tokens[:count]), "complete finite dynamic native field")
    arrays, offset = {}, 0
    for name in names:
        size = nodes if name in ("X", "Y") else cells
        arrays[name] = [float(x) for x in tokens[offset:offset + size]]
        offset += size
    connection = [int(x) - 1 for x in tokens[count:]]
    require(all(0 <= x < nodes for x in connection), "new mesh connectivity bounds")
    require(len(set(arrays["cell_id"])) == cells, "current native IDs unique")
    metadata = auxiliary(header)
    curvature_metadata(metadata)
    return cells, arrays, metadata


def dynamic_history(executable, project, positive_task, environment, reject=False):
    """One two-step surface-force fixture, plus its accepted-prefix failure."""
    task = positive_task.parent / (positive_task.name + ("_remap_failure" if reject else "_remap"))
    original = project / "cases/hydro_static_drop/system"
    hashes = {name: hashlib.sha256((original / name).read_bytes()).hexdigest() for name in FILES}
    shutil.copytree(original, task / "system")
    (task / "system/mesh.dict").write_text(
        "Mesh { resolution (8 8); initialMesh uniform; Regrid { minLevel 1; maxLevel 2; "
        "interfaceBandLayers 1; vorticityTolerance .5; coarsenRatio .25; } }\n", encoding="utf-8")
    path = task / "system/run.dict"
    text, count = re.subn(r"\bAdaptive\s*\{[^{}]*\}", "", path.read_text(encoding="utf-8"))
    require(count == 1, "copied capillary fixture has one adaptive block")
    text, count = re.subn(r"\bendTime\s+[^;]+;", "endTime .0004;", text)
    require(count == 1, "dynamic fixture two original dt steps")
    text, count = re.subn(r"\bwriteEverySteps\s+\d+;", "writeEverySteps 1;", text)
    require(count == 1, "dynamic fixture writes accepted mesh every step")
    text = select_diagnostic_series(text)
    if reject:
        text, count = re.subn(r"\bminDt\s+[^;]+;", "minDt .0002;", text)
        require(count == 1, "prefix failure retry lower bound")
        text, count = re.subn(r"(\bTime\s*\{)", lambda m: m[1] +
            "\n Adaptive { minDt .0002; maxDt .0002; targetCourant 1e-12; safety .9; maxGrowth 1; }", text)
        require(count == 1, "prefix failure explicit adaptive lower bound")
    path.write_text(text, encoding="utf-8")
    process = subprocess.run([str(executable), "--threads", "2"], cwd=task, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", timeout=120)
    (task / "contract.log").write_text(process.stdout, encoding="utf-8")
    require((process.returncode != 0) if reject else
            (process.returncode == 0 and "[COMPLETED]" in process.stdout), "dynamic capillary status; see " + str(task / "contract.log"))
    output = task / "output"
    manifest = json.loads((output / "meta/manifest.json").read_text(encoding="utf-8"))
    accepted = manifest["accepted_step"]
    require(accepted == (1 if reject else 2), "two-step dynamic success or one remapped accepted prefix")
    if reject:
        failure = json.loads((output / "meta/failure.json").read_text(encoding="utf-8"))
        require(failure["accepted_step"] == 1 and failure["accepted_state_preserved"] and
                "adaptive" in failure["reason"].lower() and "minDt" in failure["reason"], "failure preserves first real remapped step")
    _, frames = csv_rows(output / "meta/frames.csv")
    require([int(x["step"]) for x in frames] == list(range(accepted + 1)) and
            manifest["complete_frames"] == len(frames), "unique complete dynamic frames/accepted prefix")
    conservation = {int(x["step"]): x for x in load_conservation(output)}
    decisions = [decoded_row(x) for x in table_history(output / "series/regrid_decision.dat", "regrid_decision")]
    remaps = [decoded_row(x) for x in table_history(output / "series/regrid.dat", "regrid")]
    require([int(float(x["step"])) for x in decisions] == list(range(1, accepted + 1)), "one decision per enabled accepted End")
    changed = {int(float(x["step"])): x for x in decisions if float(x["changed"]) == 1}
    receipts = {int(float(x["step"])): x for x in remaps}
    require(1 in changed and set(receipts) == set(changed), "only actual topology changes have remap receipts")
    native = {int(float(x["tokens"]["step"])): x for x in table_history(output / "series/vof_native_velocity.dat", "vof_native_velocity")}
    energy = {int(float(x["tokens"]["step"])): x for x in table_history(output / "series/energy.dat", "energy")}
    hf = {int(float(x["tokens"]["step"])): x for x in table_history(output / "series/vof_csf.dat", "vof_csf")}
    for row in hf.values():
        require(row["available"]["unresolved_isolated_fragment_cells"] == 1 and float(row["tokens"]["unresolved_isolated_fragment_cells"]) == 0, "existing resolved dynamic fixture declares zero isolated-fragment closures")
    layout_sizes = {}
    for frame in frames:
        step = int(frame["step"]); cells, arrays, meta = dynamic_payload(output / frame["volume"]);layout_sizes[step] = cells
        require(int(meta["GeometryIdentity"]) == int(frame["geometry_identity"]), "dynamic frame owner")
        equal_float(float(meta["VTK_TimeValue"]), float(conservation[step]["time"]), "dynamic exact accepted field time")
        require(math.isclose(math.fsum(arrays["liquid_volume"]), float(conservation[step]["negative_phase_volume"]), rel_tol=2e-12), "new field independent liquid total")
        require(math.isclose(math.fsum(100*l+g for l,g in zip(arrays["liquid_volume"],arrays["gas_volume"])), float(conservation[step]["total_geometric_mass"]), rel_tol=2e-12), "new field independent mass total")
        if step not in receipts:
            continue
        receipt, decision = receipts[step], changed[step]
        require(int(float(receipt["target_geometry_identity"])) == int(frame["geometry_identity"]) and
                receipt["target_geometry_identity"] != receipt["source_geometry_identity"], "D1 authenticated distinct source and target IDs")
        require(float(receipt["same_end"]) == 1 and float(receipt["history_reset"]) == 1 and
                int(float(receipt["remap_solves"])) == 1, "real same End remap/history reset/projection count")
        equal_float(float(receipt["time"]), float(conservation[step]["time"]), "remap accepted time")
        equal_float(float(receipt["begin"]), float(conservation[step-1]["time"]), "remap full-clock begin")
        require(int(float(receipt["target_leaves"])) == cells and
                int(float(receipt["source_leaves"])) == int(float(decision["source_leaves"])), "receipt actual old/new leaves")
        for phase in ("liquid", "gas"):
            require(math.isclose(float(receipt["source_"+phase]), float(receipt["target_"+phase]), rel_tol=1e-12), "independent remap phase ledger")
            require(float(receipt["global_relative_"+phase]) <= 1e-12 and float(receipt["donor_relative_"+phase]) <= 1e-12, "own-phase report gates")
        require(math.isclose(float(receipt["delta_native_kinetic"]), float(receipt["target_native_kinetic"])-float(receipt["source_native_kinetic"]), rel_tol=2e-12, abs_tol=1e-24), "remap delta KE is separate native change")
        require(math.isclose(float(energy[step]["tokens"]["kinetic"]), float(receipt["target_native_kinetic"]), rel_tol=2e-12, abs_tol=1e-24), "energy uses new retained M and velocity")
        require(all(x == 0 for x in arrays["kappa_HF_force"]+arrays["kappa_HF_forceAvailable"]), "no old force array pasted on new cells")
        require(all(x == 0 for name in ("kappa_effective_force", "kappa_effective_forceAvailable", "capillary_curvature_method") for x in arrays[name]), "remap keeps old effective force and method unavailable on new cells")
        require(int(meta["VTK_hf_force_geometry_identity"]) == int(float(hf[step]["tokens"]["force_geometry_identity"])) != 0 and
                int(meta["VTK_hf_force_geometry_identity"]) != int(frame["geometry_identity"]), "original physical force owner retained")
        for group in ("global", "pure_liquid", "pure_gas", "interface", "geometry_pure_gas"):
            for field in ("midpoint_capillary_load", "midpoint_force_kappa"):
                require(native[step]["available"][group+"_"+field] == 0, "old per-DOF force is unavailable on new native velocity peak")
        require(native[step]["available"]["force_time"] == 1 and native[step]["available"]["force_geometry_identity"] == 1, "scalar force metadata remains available")
        for stage, expected_cells in (("VOF_AC",int(float(receipt["source_leaves"]))), ("vof_hf",int(float(receipt["source_leaves"]))), ("io",cells), ("step_total",cells)):
            rows = table_history(output / ("series/performance_"+stage+".dat"), "performance")
            value = next(x for x in rows if int(float(x["tokens"]["step"])) == step)
            require(int(float(value["tokens"]["cells"])) == expected_cells, "performance uses real source or target layout")
        if frame["interface"]:
            interface_dat(output / frame["interface"], unavailable=True)
            interface_text = (output / frame["interface"]).read_text(encoding="utf-8")
            require(auxiliary(interface_text).get("RemappedForceUnavailable", "").startswith("old physical force"), "new interface explicitly carries unavailable old curvature")
    require(layout_sizes[0] != layout_sizes[1], "actual first accepted topology change")
    require(not list((output / "series").glob("*.csv")), "dynamic output remains native DAT only")
    info = (output / "meta/build_info.txt").read_text(encoding="utf-8")
    require("never external work" in info and "original physical-step owner" in info, "physical/remap budget provenance disclosed")
    require("capillary_curvature=" in info and "method3 is unresolved isolated zero closure" in info, "output metadata does not label closure as physical curvature")
    require(all(hashlib.sha256((original/name).read_bytes()).hexdigest() == value for name,value in hashes.items()), "original capillary six inputs unchanged")
    return task


def main():
    require(len(sys.argv) == 4, "usage: hydro_tecplot_output_contract.py LAMP.exe projectRoot workRoot")
    executable, project, work = (Path(value).resolve() for value in sys.argv[1:])
    require(executable.is_file() and work.drive.upper() == "D:", "formal executable and D: working directory required")
    task = work / ("hydro_tecplot_" + uuid.uuid4().hex[:12])
    task.mkdir(parents=True, exist_ok=False)
    original = project / "cases/hydro_rayleigh_taylor/system"
    hashes = {name: hashlib.sha256((original / name).read_bytes()).hexdigest() for name in FILES}
    shutil.copytree(original, task / "system")
    # This fixture explicitly exercises fixed dt; adaptive histories below
    # add their own block to independent copies of this fixed input.
    run_path = task / "system/run.dict"
    fixed_text, adaptive_blocks = re.subn(r"\bAdaptive\s*\{[^{}]*\}", "",
        run_path.read_text(encoding="utf-8"))
    require(adaptive_blocks <= 1, "at most one production Adaptive block in RT copy")
    if adaptive_blocks:
        run_path.write_text(fixed_text, encoding="utf-8")
    require({path.name for path in (task / "system").iterdir()} == set(FILES), "actual RT six-file input")
    for name, pattern, replacement in (("mesh.dict", r"\bresolution\s*\(\s*\d+\s+\d+\s*\)", "resolution (8 32)"),
                                       ("run.dict", r"\bendTime\s+[^;]+;", "endTime .004;"),
                                       ("run.dict", r"\bwriteEverySteps\s+\d+;", "writeEverySteps 1;")):
        path = task / "system" / name
        value, count = re.subn(pattern, replacement, path.read_text(encoding="utf-8"))
        require(count == 1, "unique test-only input edit: " + pattern)
        path.write_text(value, encoding="utf-8")
    temporary = task / "_temp"
    temporary.mkdir()
    environment = os.environ.copy()
    environment.update(TEMP=str(temporary), TMP=str(temporary), OMP_NUM_THREADS="2", MKL_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.run([str(executable), "--threads", "2"], cwd=task, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, encoding="utf-8", errors="replace", timeout=120)
    (task / "contract.log").write_text(process.stdout, encoding="utf-8")
    require(process.returncode == 0 and "[COMPLETED]" in process.stdout, "two-step real RT run failed; see " + str(task / "contract.log"))
    progress = [line for line in process.stdout.splitlines() if line.startswith("[Step ")]
    require(progress and all(not any(label in line for label in ("surface_area=", "energy_area=", "energy_closure=", "plic_sweeps=", "force_geometry=")) for line in progress), "summary progress contains supported concise diagnostics")
    require("[ADVANCE]" not in process.stdout, "advance diagnostics remain detail-only")
    require("threads=2 requested_threads=2" in process.stdout and "OMP: Warning #96" not in process.stdout, "Core loops obey the configured two-thread runtime without fixed-12 requests")
    require(all(hashlib.sha256((original / name).read_bytes()).hexdigest() == value for name, value in hashes.items()), "original RT dictionaries remain unchanged")
    output = task / "output"
    manifest, frames = check_volume_output(output, (0, 1, 2))
    for index, frame in enumerate(frames):
        equal_float(float(frame["time"]), .002 * index, "two original-dt RT steps")
    accepted_rows = native_series(output, (0, 1, 2))
    diagnostic_task = diagnostic_profile(executable, task, environment)
    failure_task = failure_prefix(executable, task, environment)
    adaptive_task = adaptive_history(executable, task, environment)
    adaptive_failure_task = adaptive_history(executable, task, environment, reject=True)
    remap_task = dynamic_history(executable, project, task, environment)
    remap_failure_task = dynamic_history(executable, project, task, environment, reject=True)
    require(all(hashlib.sha256((original / name).read_bytes()).hexdigest() == value for name, value in hashes.items()), "positive/negative copies do not alter original RT inputs")
    print(f"HYDRO_TECPLOT_OUTPUT_PASS frames=3 accepted_table_rows={accepted_rows} work={task} diagnostic_work={diagnostic_task} failure_work={failure_task} adaptive_work={adaptive_task} adaptive_failure_work={adaptive_failure_task} remap_work={remap_task} remap_failure_work={remap_failure_task}")


if __name__ == "__main__":
    main()
