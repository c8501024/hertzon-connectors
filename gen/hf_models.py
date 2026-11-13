"""Step 2: simulation models for Hertzon, and Hertzon projects of the detailed STEP models.

python gen/hf_models.py --hertzon /path/to/Hertzon   ->  hzn/<model>_HF.hzn, hzn/<model>_detailed.hzn

<model>_HF.hzn: simplified, fully parametric full-wave model of connector + launch + board line, built from
Hertzon primitives with design variables (b, a, bead, launch pin, board h, w, g ...), wave port 1 on the coax
reference plane and wave port 2 on the board line, a solution setup at f_max and a sweep.

What is simplified, for an FEM solve up to the connector's f_max (and why it does not change S):
  * "negative" modelling: only the field volumes (coax air/PTFE, bead, launch dielectric, channel, substrate,
    air above the line) are drawn; their outer faces are the metal (Perfect E, the Hertzon default for outer
    faces). Connector body, block, screws and threads are outside the field and are not drawn.
  * coupling thread, socket slots, chamfers and screw holes: outside the field, or below lambda/50 at f_max.
  * the socket is replaced by a solid pin (the mating plane is the port; the socket is behind it in the
    mated pair and its slots are below lambda/100).
  * kept because they set S11: bead with the necked pin (compensation), the step from the coax to the launch
    pin and launch dielectric, launch pin length on the trace, channel over the launch, GCPW gap and ground.
  * metals: centre conductor 'gold' (BeCu with >= 1.27 um Au, > 1.6 skin depths at 10 GHz); outer walls PEC.
    For insertion loss, assign Finite Conductivity (stainless steel 1.1e6 S/m for CRES housings, gold for
    plated brass) to the coax outer faces.

<model>_detailed.hzn: every body of step/<model>.step as an Import object with its material (relative path, so
keep the hzn/ and step/ folders side by side).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectors_db import INTERFACES, DIELECTRICS, cpwg, cpwg_width, line_for_board  # noqa: E402
from geometry import MODELS, derived, Params  # noqa: E402

LOCAL_MATERIALS = [
    # name, er, tand, color, description
    ("PTFE_connector", "2.05", "0.0003", "#f4f4f4", "solid PTFE in RF connectors (typ. 2.02-2.08)"),
    ("Ultem_1000", "3.15", "0.0025", "#c8732a", "PEI support bead (SABIC Ultem 1000, typical at 1-10 GHz)"),
    ("Kel-F_PCTFE", "2.35", "0.005", "#d9e6f2", "PCTFE support bead (typical)"),
]
HZ_DIEL = {"PTFE": "PTFE_connector", "Ultem 1000 (PEI)": "Ultem_1000", "Kel-F (PCTFE)": "Kel-F_PCTFE",
           "PEEK": "PEEK", "air": "vacuum"}
BOARD = {"Rogers RO4003C (tm)": 3.55}


# ------------------------------------------------------------------------------------------------------------
# line synthesis (no scipy): complete elliptic integral by the AGM, conductor-backed CPW (Ghione / Wadell)
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
def _hz(path: str):
    sys.path.insert(0, path)
    from hertzon.model.design import Project, Command, Assignment, FaceRef, Setup, Sweep, Report  # noqa: F401
    from hertzon.model import schema
    from hertzon.model.materials import Material
    return Project, Command, Assignment, FaceRef, Setup, Sweep, Report, schema, Material


def mm(v: float) -> str:
    return f"{v:.4g}mm"


def build_hf(p: Params, H) -> tuple[object, dict]:
    Project, Command, Assignment, FaceRef, Setup, Sweep, Report, schema, Material = H
    d = derived(p)
    f = INTERFACES[p.interface]
    fmax = {"SMA": 27.0 if p.template == "flange_block" else 18.0, "2.92mm": 40.0, "2.4mm": 50.0}[p.interface]
    er_b = BOARD["Rogers RO4003C (tm)"]
    h = p.board_h
    w, g = d["line_w"], d["line_g"]
    z_line, eeff = cpwg(w, g, h, er_b)
    proj = Project(f"{p.name}_HF")
    proj.materials = [Material(name=n, permittivity=e, loss_tangent=t, color=c, library="Project", description=s)
                      for n, e, t, c, s in LOCAL_MATERIALS]
    D = proj.designs[0]
    D.name = "EdgeLaunch"
    D.settings["MaterialOverride"] = "True"
    v = D.variables
    lam_min = 299.792458 / fmax / math.sqrt(er_b)          # mm, in the substrate
    v.update({"b": mm(f["outer_conductor_id_mm"]), "a": mm(f["center_conductor_od_mm"]), "h": mm(h),
              "w": mm(w), "g": mm(g),
              # the trace is a sheet at z = 0: the launch pin / tab rests on it (copper thickness neglected)
              "zax": mm(d["launch_pin"] / 2 if p.template == "flange_block" else d["tab_t_est"] / 2), "Lline": mm(max(3.0, 0.6 * lam_min * 1.5)),
              "Lport": mm(max(1.5, 0.35 * 299.792458 / fmax))})
    notes = {"line": f"GCPW on RO4003C h={h} mm: w={w:.3f} mm, g={g} mm (analytic open-board Z0 {z_line:.1f} ohm"
                     + ("; width calibrated to 50 ohm port impedance in the shielded model)" if p.template == "flange_block" else ")")}
    pin_mat = "gold"
    objs, pe = [], []

    def add(kind, name, material="vacuum", **props):
        o = D.add_object(kind, name=name, material=material, **{k: str(vv) for k, vv in props.items()})
        objs.append(o)
        return o

    bead_mat = HZ_DIEL[p.bead]
    if p.template == "flange_block":
        Lw, Lf, Lt = d["wall_est"], d["flange_t_est"], d["thread_len_est"]
        v.update({"Lw": mm(Lw), "Lcx": mm(Lf + Lt), "bl": mm(d["bead_len_est"]), "bz": mm(d["bead_z_est"]),
                  "neck": mm(d["neck"]), "dp": mm(d["launch_pin"]), "dd": mm(d["launch_diel"]),
                  "Lover": mm(d["launch_over_est"]), "wc": "w+2*g+0.3mm", "hc": mm(d["channel_h_est"]),
                  "Lch": mm(d["channel_len_est"]), "Wbox": "w+2*g+2*0.6mm", "Hair": "max(1.5mm, 6*h)"})
        cx_mat = "PTFE_connector" if p.interface == "SMA" else "vacuum"
        # coax from the reference plane (y = -(Lw+Lcx)) to the flange rear face (y = -Lw)
        add("CreateCylinder", "CoaxFront", cx_mat, CenterPosition="0,-(Lw+Lcx),zax", Axis="Y", Radius="b/2",
            Height="bz-bl/2")
        add("CreateCylinder", "CoaxRear", cx_mat, CenterPosition="0,-(Lw+Lcx)+bz+bl/2,zax", Axis="Y", Radius="b/2",
            Height="Lcx-bz-bl/2")
        add("CreateCylinder", "Bead", bead_mat, CenterPosition="0,-(Lw+Lcx)+bz-bl/2,zax", Axis="Y", Radius="b/2",
            Height="bl")
        add("CreateCylinder", "LaunchDiel", "PTFE_connector", CenterPosition="0,-Lw,zax", Axis="Y", Radius="dd/2",
            Height="Lw")
        add("CreateBox", "Channel", "vacuum", Position="-wc/2,0,0", XSize="wc", YSize="Lch", ZSize="hc")
        add("CreateBox", "Air", "vacuum", Position="-Wbox/2,Lch,0", XSize="Wbox", YSize="Lline", ZSize="Hair")
        add("CreateBox", "Substrate", "Rogers RO4003C (tm)", Position="-Wbox/2,0,-h", XSize="Wbox",
            YSize="Lch+Lline", ZSize="h")
        # centre conductor: front pin, neck in the bead, rear pin, taper, launch pin on the trace
        add("CreateCylinder", "Pin", pin_mat, CenterPosition="0,-(Lw+Lcx),zax", Axis="Y", Radius="a/2",
            Height="bz-bl/2")
        add("CreateCylinder", "PinNeck", pin_mat, CenterPosition="0,-(Lw+Lcx)+bz-bl/2,zax", Axis="Y",
            Radius="neck/2", Height="bl")
        add("CreateCylinder", "PinRear", pin_mat, CenterPosition="0,-(Lw+Lcx)+bz+bl/2,zax", Axis="Y",
            Radius="a/2", Height="Lcx-bz-bl/2-0.3mm")
        add("CreateCone", "PinTaper", pin_mat, CenterPosition="0,-Lw-0.3mm,zax", Axis="Y", LowerRadius="a/2",
            UpperRadius="dp/2", Height="0.3mm")
        add("CreateCylinder", "LaunchPin", pin_mat, CenterPosition="0,-Lw,zax", Axis="Y", Radius="dp/2",
            Height="Lw+Lover")
        trace = add("CreateRectangle", "Trace", Position="-w/2,0,0", Axis="Z", XSize="w", YSize="Lch+Lline")
        gl = add("CreateRectangle", "GroundL", Position="-Wbox/2,Lch,0", Axis="Z", XSize="Wbox/2-w/2-g",
                 YSize="Lline")
        gr = add("CreateRectangle", "GroundR", Position="w/2+g,Lch,0", Axis="Z", XSize="Wbox/2-w/2-g", YSize="Lline")
        glc = add("CreateRectangle", "GroundChL", Position="-wc/2,0,0", Axis="Z", XSize="wc/2-w/2-g", YSize="Lch")
        grc = add("CreateRectangle", "GroundChR", Position="w/2+g,0,0", Axis="Z", XSize="wc/2-w/2-g", YSize="Lch")
        pe = [trace.name, gl.name, gr.name, glc.name, grc.name]
        p1 = add("CreateCircle", "Port1", CenterPosition="0,-(Lw+Lcx),zax", Axis="Y", Radius="b/2")
        p2 = add("CreateRectangle", "Port2", Position="-Wbox/2,Lch+Lline,-h", Axis="Y", XSize="Wbox", ZSize="h+Hair")
        notes["box"] = ("Air, channel and substrate are closed by Perfect E walls: the block, the via fence of the "
                        "GCPW and a lid at Hair above the line. Wbox = w+2g+1.2 mm keeps box modes above f_max.")
        rad = []
    else:
        Lb, Lt = d["body_len_est"], d["thread_len_est"]
        Wa = max(d["body_w_est"], 2 * d["leg_x_est"] + d["leg_w_est"] + 1.0)
        v.update({"Lcx": mm(Lb + Lt), "bl": mm(d["bead_len_est"]), "bz": mm(d["bead_z_est"]), "neck": mm(d["neck"]),
                  "tw": mm(d["tab_w_est"]), "tt": mm(d["tab_t_est"]), "tl": mm(d["tab_len_est"]),
                  "lw": mm(d["leg_w_est"]), "lt": mm(d["leg_t_est"]), "ll": mm(d["leg_len_est"]),
                  "lx": mm(d["leg_x_est"]), "Wa": mm(Wa), "Hair": "max(3mm, 0.25*299.79mm/" + f"{fmax:g})",
                  "Lline": mm(max(6.0, d["leg_len_est"] + 3.0))})
        if p.interface == "SMA":
            add("CreateCylinder", "Coax", "PTFE_connector", CenterPosition="0,-Lcx,zax", Axis="Y", Radius="b/2",
                Height="Lcx")
        else:
            add("CreateCylinder", "CoaxFront", "vacuum", CenterPosition="0,-Lcx,zax", Axis="Y", Radius="b/2",
                Height="bz-bl/2")
            add("CreateCylinder", "CoaxRear", "vacuum", CenterPosition="0,-Lcx+bz+bl/2,zax", Axis="Y", Radius="b/2",
                Height="Lcx-bz-bl/2-0.9mm")
            add("CreateCylinder", "Bead", bead_mat, CenterPosition="0,-Lcx+bz-bl/2,zax", Axis="Y", Radius="b/2",
                Height="bl")
            add("CreateCylinder", "RearInsulator", "PTFE_connector", CenterPosition="0,-0.9mm,zax", Axis="Y",
                Radius="b/2", Height="0.9mm")
        add("CreateBox", "Air", "vacuum", Position="-Wa/2,0,0", XSize="Wa", YSize="Lline", ZSize="Hair")
        add("CreateBox", "Substrate", "Rogers RO4003C (tm)", Position="-Wa/2,0,-h", XSize="Wa", YSize="Lline",
            ZSize="h")
        if p.interface == "SMA":
            add("CreateCylinder", "Pin", pin_mat, CenterPosition="0,-Lcx,zax", Axis="Y", Radius="a/2", Height="Lcx")
        else:
            add("CreateCylinder", "Pin", pin_mat, CenterPosition="0,-Lcx,zax", Axis="Y", Radius="a/2",
                Height="bz-bl/2")
            add("CreateCylinder", "PinNeck", pin_mat, CenterPosition="0,-Lcx+bz-bl/2,zax", Axis="Y",
                Radius="neck/2", Height="bl")
            add("CreateCylinder", "PinRear", pin_mat, CenterPosition="0,-Lcx+bz+bl/2,zax", Axis="Y", Radius="a/2",
                Height="Lcx-bz-bl/2")
        add("CreateBox", "Tab", pin_mat, Position="-tw/2,0,zax-tt/2", XSize="tw", YSize="tl", ZSize="tt")
        for s, x in (("L", "-lx-lw/2"), ("R", "lx-lw/2")):
            add("CreateBox", f"Leg{s}", "gold", Position=f"{x},0,0", XSize="lw", YSize="ll", ZSize="lt")
        trace = add("CreateRectangle", "Trace", Position="-w/2,0,0", Axis="Z", XSize="w", YSize="Lline")
        gl = add("CreateRectangle", "GroundL", Position="-Wa/2,0,0", Axis="Z", XSize="Wa/2-w/2-g", YSize="Lline")
        gr = add("CreateRectangle", "GroundR", Position="w/2+g,0,0", Axis="Z", XSize="Wa/2-w/2-g", YSize="Lline")
        pe = [trace.name, gl.name, gr.name]
        # via fence: SMA 0.3 mm vias, pitch 0.8 mm, 0.35 mm from the gap; 2.92/2.4 mm 0.25 mm, pitch 0.5 mm, 0.2 mm
        dense = p.interface != "SMA"
        vp_ = 0.5 if dense else 0.8
        v.update({"vd": "0.25mm" if dense else "0.3mm", "vp": mm(vp_), "vo": "w/2+g+" + ("0.2mm" if dense else "0.35mm")})
        n_v = int(max(6.0, d["leg_len_est"] + 3.0) / vp_)
        for i in range(n_v):
            for s, sx in (("L", "-vo"), ("R", "vo")):
                add("CreateCylinder", f"Via{s}{i + 1}", "pec", CenterPosition=f"{sx},vp/2+{i}*vp,-h", Axis="Z",
                    Radius="vd/2", Height="h")
        p1 = add("CreateCircle", "Port1", CenterPosition="0,-Lcx,zax", Axis="Y", Radius="b/2")
        p2 = add("CreateRectangle", "Port2", Position="-Wa/2,Lline,-h", Axis="Y", XSize="Wa", ZSize="h+Hair")
        top = add("CreateRectangle", "RadTop", Position="-Wa/2,0,Hair", Axis="Z", XSize="Wa", YSize="Lline")
        sl = add("CreateRectangle", "RadSideL", Position="-Wa/2,0,-h", Axis="X", YSize="Lline", ZSize="h+Hair")
        sr = add("CreateRectangle", "RadSideR", Position="Wa/2,0,-h", Axis="X", YSize="Lline", ZSize="h+Hair")
        rad = [top.name, sl.name, sr.name]
        notes["box"] = ("Open board: Radiation on the top and side faces of the air box; the plane y = 0 (connector "
                        "body face), the substrate bottom and the coax walls are Perfect E. Bottom legs lie on "
                        "the ground and are not drawn.")

    def assign(lst, kind, name, objects, **props):
        pr = schema.defaults((schema.BOUNDARIES if lst is D.boundaries else schema.EXCITATIONS)[kind])
        pr.update({k: str(x) for k, x in props.items()})
        lst.append(Assignment(name, kind, pr, faces=[], objects=list(objects)))

    assign(D.boundaries, "Perfect E", "PerfE_Traces", pe)
    if rad:
        assign(D.boundaries, "Radiation", "Rad", rad)
    assign(D.excitations, "Wave Port", "Coax", [p1.name], PortNumber=1)
    assign(D.excitations, "Wave Port", "GCPW", [p2.name], PortNumber=2)
    s = Setup("Setup1")
    s.props.update({"Frequency": f"{fmax:g}GHz", "MaximumPasses": "12", "MaximumDeltaS": "0.02"})
    sw = Sweep("Sweep")
    sw.props.update({"SweepType": "Interpolating", "Distribution": "Linear Count", "Start": "0.5GHz",
                     "Stop": f"{fmax:g}GHz", "Count": str(int(fmax * 4) + 1)})
    s.sweeps.append(sw)
    D.setups.append(s)
    D.reports.append(Report("Return and insertion loss", "Rectangular", "Setup1 : Sweep",
                            ["dB(S(1,1))", "dB(S(2,1))", "dB(S(2,2))"]))
    notes.update({"f_max_GHz": fmax, "interface": p.interface, "template": p.template,
                  "coax_z0_ohm": f["z0_ohm"], "bead": p.bead, "bead_neck_mm": round(d["neck"], 3)})
    if p.template == "flange_block":
        from connectors_db import z_coax
        notes["launch_section_z0_ohm"] = round(z_coax(d["launch_diel"], d["launch_pin"], 2.05), 1)
    return proj, notes


def build_detailed(p: Params, meta: dict, H) -> object:
    Project, Command, Assignment, FaceRef, Setup, Sweep, Report, schema, Material = H
    proj = Project(f"{p.name}_detailed")
    proj.materials = [Material(name=n, permittivity=e, loss_tangent=t, color=c, library="Project", description=s)
                      for n, e, t, c, s in LOCAL_MATERIALS]
    D = proj.designs[0]
    D.name = "Detailed"
    rel = f"../step/{meta['step']}"
    for bdy in meta["bodies"]:
        mat = {"Teflon (tm)": "PTFE_connector"}.get(bdy["hertzon_material"], bdy["hertzon_material"])
        if mat == "stainless steel":
            mat = "stainless_steel"
        D.add_object("Import", name=bdy["name"], material=mat, File=rel, Body=str(bdy["index"]))
    return proj


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--hertzon", required=True, help="path of a Hertzon checkout (for its model classes)")
    ap.add_argument("--validate", action="store_true", help="run Hertzon's model check on every project")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    H = _hz(a.hertzon)
    out = ROOT / "hzn"
    out.mkdir(exist_ok=True)
    summary = {}
    for p in MODELS:
        if a.only and a.only not in p.name:
            continue
        proj, notes = build_hf(p, H)
        f_hf = out / f"{p.name}_HF.hzn"
        proj.save(f_hf)
        meta = json.loads((ROOT / "step" / f"{p.name}.json").read_text())
        f_det = out / f"{p.name}_detailed.hzn"
        build_detailed(p, meta, H).save(f_det)
        summary[p.name] = notes
        msg = ""
        if a.validate:
            from hertzon.model.design import Project
            from hertzon.model.validate import validate
            from hertzon.mcp.tools import tessellate
            res = []
            for f in (f_hf, f_det):
                pr = Project.load(f)
                de = pr.designs[0]
                tess, errors = tessellate(pr, de)
                bad = [c.message for c in validate(pr, de, tess, errors) if c.level == "error"]
                warn = [c.message for c in validate(pr, de, tess, errors) if c.level == "warning"]
                res.append(f"{f.name}: {'OK' if not bad else 'ERRORS ' + str(bad)}" + (f" warnings {warn}" if warn else ""))
            msg = " | ".join(res)
        print(p.name, notes["line"], msg)
    (out / "hf_notes.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
