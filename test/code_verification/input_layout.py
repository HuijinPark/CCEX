"""The input layout of A64/A65 (2026-10-15) for the scripts of test/code_verification.

Many checks build an input from flat keys (outfile, backend, cluster_mode, ...) or edit one with them, as the inputs were
written before 2026-10-15. flat(d) lifts the blocks of a current input to those flat names, so an edit written for the
flat form still applies; nest(d) puts every flat name back into its block, as CCEX reads it now:

  configuration  gyrofile, bathfile, bath_shift, bath_state_file, nstate, seed, and Defect.axis_file / Defect.state_file
                 as defect_axis_file / defect_state_file; the per-bath-file names as lists (A64)
  execution      backend, gpu_options, propagator, evolution, evolution_update, gcce_refinement, column_evolution, profile
  cluster        cluster_mode -> mode, add_subclusters / max_clusters_per_order -> selection.{add_subclusters, max_per_order},
                 pcce, bundle, adaptive, probe_memory_gb -> adaptive.probe_memory_gb
  output         outfile -> file, savemode, cluster_contributions, cluster_only
  frame          coordinate_frame_rotation -> rotation
  Defect.species n_nuclear_spins -> n_onsite_spins, n_configurations -> n_axis_states

A94 (2026-10-16): nest also writes the forms A94 asks for: configuration.nstate N as s_range [1, N] (0: left out, the
ensemble), the Qubit.file form ({"file", "spin", "name"} in the Qubit block) as a Qubit.qubit entry, a Defect species'
axis as zfs.axis of the D form (left out without one), pulse.bath.alphams / betams as levels (pulse.bath.spin left out).

A refusal check of an old place keeps it: the names in d["__keep_old__"] (a list) are left where they are.
"""
CONFIGURATION = ("gyrofile", "bathfile", "bath_shift", "bath_state_file", "nstate", "s_range", "seed")   # nstate: s_range (A94)
EXECUTION = ("backend", "gpu_options", "propagator", "evolution", "evolution_update", "gcce_refinement",
             "column_evolution", "profile")
CLUSTER = {"cluster_mode": "mode", "pcce": "pcce", "bundle": "bundle", "adaptive": "adaptive"}
SELECTION = {"add_subclusters": "add_subclusters", "max_clusters_per_order": "max_per_order"}
OUTPUT = {"outfile": "file", "savemode": "savemode", "cluster_contributions": "cluster_contributions",
          "cluster_only": "cluster_only"}
SPECIES = {"n_nuclear_spins": "n_onsite_spins", "n_configurations": "n_axis_states"}


def nest(d):
    keep = set(d.pop("__keep_old__", ()))

    def block(name):
        b = d.get(name)
        if not isinstance(b, dict):
            b = d[name] = {}
        return b

    for k in CONFIGURATION:
        if k in d and k not in keep:
            block("configuration")[k] = d.pop(k)
    defect = d.get("Defect")
    if isinstance(defect, dict):
        for k in ("axis_file", "state_file"):
            if k in defect and k not in keep:
                block("configuration")["defect_" + k] = defect.pop(k)
        for sp in defect.get("species", []) if isinstance(defect.get("species"), list) else []:
            for old, new in SPECIES.items():
                if isinstance(sp, dict) and old in sp and old not in keep:
                    sp[new] = sp.pop(old)
            # A94 (G5): the symmetry axis is zfs.axis of the D form (left out where nothing used it)
            if isinstance(sp, dict) and "axis" in sp and "axis" not in keep:
                axis, zfs = sp.pop("axis"), sp.get("zfs")
                if isinstance(zfs, dict) and "D" in zfs:
                    zfs.setdefault("axis", axis)
    conf = d.get("configuration")
    if isinstance(conf, dict) and "nstate" in conf and "nstate" not in keep:      # A94
        n = conf.pop("nstate")
        if n:
            conf["s_range"] = [1, n]
    q = d.get("Qubit")
    if isinstance(q, dict) and "file" in q and "Qubit.file" not in keep:          # A94 (G12)
        old = q.get("qubit")
        entry = dict(old[0]) if isinstance(old, list) and len(old) == 1 and isinstance(old[0], dict) else {}
        entry.pop("xyz", None)
        entry["file"] = q.pop("file")
        for k in ("name", "spin"):
            if k in q:
                entry[k] = q.pop(k)
        q["qubit"] = [entry]
    pb = d.get("pulse", {}).get("bath") if isinstance(d.get("pulse"), dict) else None
    if isinstance(pb, dict) and "pulse.bath" not in keep:                         # A94 (G6)
        pb.pop("spin", None)
        if "alphams" in pb or "betams" in pb:
            pb["levels"] = [pb.pop("alphams", 0.5), pb.pop("betams", -0.5)]
    if isinstance(conf, dict):
        for k in ("bath_state_file", "defect_axis_file", "defect_state_file"):
            if isinstance(conf.get(k), str) and k not in keep:
                conf[k] = [conf[k]]
    for k in EXECUTION:
        if k in d and k not in keep:
            block("execution")[k] = d.pop(k)
    for old, new in CLUSTER.items():
        if old in d and old not in keep:
            v = d.pop(old)
            c = block("cluster")
            if new == "adaptive" and isinstance(c.get("adaptive"), dict) and isinstance(v, dict) \
                    and "probe_memory_gb" in c["adaptive"] and "probe_memory_gb" not in v:
                v = dict(v, probe_memory_gb=c["adaptive"]["probe_memory_gb"])   # the block replaced, its probe memory kept
            c[new] = v
    for old, new in SELECTION.items():
        if old in d and old not in keep:
            c = block("cluster")
            c.setdefault("selection", {})[new] = d.pop(old)
    if "probe_memory_gb" in d and "probe_memory_gb" not in keep:
        block("cluster").setdefault("adaptive", {})["probe_memory_gb"] = d.pop("probe_memory_gb")
    for old, new in OUTPUT.items():
        if old in d and old not in keep:
            block("output")[new] = d.pop(old)
    if "coordinate_frame_rotation" in d and "coordinate_frame_rotation" not in keep:
        block("frame")["rotation"] = d.pop("coordinate_frame_rotation")
    return d


def flat(d):
    """The inverse of nest for the A65 blocks (the configuration block stays: the scripts edit it as it is)."""
    e = d.pop("execution", None)
    if isinstance(e, dict):
        d.update(e)
    c = d.pop("cluster", None)
    if isinstance(c, dict):
        for old, new in CLUSTER.items():
            if new in c:
                d[old] = c.pop(new)
        sel = c.pop("selection", None)
        if isinstance(sel, dict):
            for old, new in SELECTION.items():
                if new in sel:
                    d[old] = sel.pop(new)
        if isinstance(d.get("adaptive"), dict) and "probe_memory_gb" in d["adaptive"]:
            d["probe_memory_gb"] = d["adaptive"].pop("probe_memory_gb")
            if not d["adaptive"]:
                d.pop("adaptive")
    o = d.get("output")
    if isinstance(o, dict):
        for old, new in OUTPUT.items():
            if new in o:
                d[old] = o.pop(new)
        if not o:
            d.pop("output")
    f = d.pop("frame", None)
    if isinstance(f, dict) and "rotation" in f:
        d["coordinate_frame_rotation"] = f["rotation"]
    defect = d.get("Defect")
    if isinstance(defect, dict) and isinstance(defect.get("species"), list):
        for sp in defect["species"]:
            for old, new in SPECIES.items():
                if isinstance(sp, dict) and new in sp:
                    sp[old] = sp.pop(new)
    return d
