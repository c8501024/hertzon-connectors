"""Step 1: detailed (mechanical) models of edge-launch connectors, with real helical coupling threads.

python gen/geometry.py                 -> step/<model>.step (+ <model>.json: body order, names, materials)
python gen/geometry.py --only sma_legs_tab_jack

Built with gmsh's OpenCASCADE kernel (the kernel Hertzon uses), so the STEP files load with
Modeler > Import STEP in Hertzon (one Import object per body; the JSON lists body index -> part -> material).

Coordinates (mm, as in Hertzon): the board lies in the XY plane with its top copper at z = 0 and the board in
y > 0; the board edge is the plane y = 0; the connector axis runs along -Y at x = 0, z = z_axis.

What is exact and what is estimated:
  * mating interface (outer conductor ID, centre conductor OD, male pin / socket bore, coupling thread form and
    pitch): from the interface standards (IEEE 287, MIL-STD-348, IEC 61169), see data/interfaces.json
  * launch pin and launch dielectric diameters (Southwest Microwave): from the catalogue
  * body, flange, block, legs, bead position: typical values; the vendors publish them only as drawings.
    For a vendor-exact body use the vendor STEP from data/connectors.json (cad.step_url) and keep this model's
    inside (coax, bead, launch), which is what the fields see.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

import gmsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectors_db import INTERFACES, DIELECTRICS, SWM_LAUNCH, IN, d_inner_for_z, line_for_board  # noqa: E402

o = None


# ------------------------------------------------------------------------------------------------------------
# parameters
# ------------------------------------------------------------------------------------------------------------
@dataclass
class Params:
    name: str
    interface: str                 # "SMA", "2.92mm", "2.4mm"
    template: str                  # "legs_tab" or "flange_block"
    board_h: float = 0.508         # substrate thickness (mm)
    t_cu: float = 0.035            # copper thickness
    # launch (flange_block): Southwest Microwave pin option ".010", ".007", ".005"
    launch_option: str = ".010"
    bead: str = "Ultem 1000 (PEI)"
    narrow: bool = False           # flange_block: narrow block (.350 in)
    n_legs: int = 4                # legs_tab: 2 (one above, one below) or 4
    thread_segments: int = 24      # points per turn of the helix
    origin: dict = field(default_factory=dict)

    @property
    def iface(self) -> dict:
        return INTERFACES[self.interface]


def derived(p: Params) -> dict:
    """All dimensions of the model (mm). Keys ending in _est are estimates."""
    f = p.iface
    b, a = f["outer_conductor_id_mm"], f["center_conductor_od_mm"]
    d = {"b": b, "a": a, "pin_male": f["male_pin_od_mm"], "thr_D": f["thread_major_mm"], "thr_P": f["thread_pitch_mm"]}
    d["socket_bore"] = d["pin_male"] + 0.0   # closed socket fingers grip the pin
    d["socket_depth_est"] = {"SMA": 2.3, "2.92mm": 2.3, "2.4mm": 1.9}[p.interface]
    d["thread_len_est"] = {"SMA": 5.6, "2.92mm": 5.6, "2.4mm": 5.0}[p.interface]
    d["bead_len_est"] = 0.9
    # bead centre behind the reference plane: behind the socket fingers
    d["bead_z_est"] = d["socket_depth_est"] + 0.5 + d["bead_len_est"] / 2
    d["neck"] = d_inner_for_z(b, DIELECTRICS[p.bead]["er"])
    d["line_w"], d["line_g"] = line_for_board(p.board_h, shielded=p.template == "flange_block",
                                              narrow=p.template == "legs_tab" and p.interface != "SMA")  # 50 ohm GCPW, RO4003C
    if p.template == "legs_tab":
        d["tab_w_est"] = 0.76
        d["tab_t_est"] = 0.30
        d["tab_len_est"] = 3.0
        d["z_axis"] = p.t_cu + d["tab_t_est"] / 2
        d["body_len_est"] = 3.0
        d["body_w_est"] = 6.35 if p.interface != "2.4mm" else 7.9
        d["leg_w_est"], d["leg_t_est"], d["leg_len_est"] = 1.0, 0.5, 3.8
        d["leg_x_est"] = 2.4 if p.interface != "2.4mm" else 2.9
    else:
        dp, dd = SWM_LAUNCH[p.launch_option]
        d["launch_pin"], d["launch_diel"] = dp * IN, dd * IN
        d["z_axis"] = p.t_cu + d["launch_pin"] / 2
        d["block_w_est"] = 0.350 * IN if p.narrow else 0.500 * IN
        d["block_len_est"] = 6.35            # along the board, behind the edge
        d["block_top_est"] = 3.2             # above the board
        d["block_bot_est"] = 1.6             # bottom clamping plate under the board
        # launch section (pin in PTFE, 72-77 ohm) kept short: 1.5 mm gave S11 -12 dB at 11 GHz in Hertzon,
        # 0.3 mm with a narrow, low channel gives about -17 dB or better to 40 GHz (2.92 mm, 0.203 mm board)
        d["wall_est"] = 0.3                  # launch dielectric length (block front wall ahead of the board edge)
        d["launch_over_est"] = 0.4           # launch pin length on the trace
        d["channel_w_est"] = max(0.6, d.get("line_w", 0.407) + 2 * d.get("line_g", 0.15) + 0.3)
        d["channel_h_est"] = max(0.4, 2.5 * p.board_h)
        d["channel_len_est"] = 1.8
        d["flange_t_est"] = 2.0
        d["flange_h_est"] = 7.6
        d["screw_dx_est"] = d["block_w_est"] / 2 - 1.7
    return d


# ------------------------------------------------------------------------------------------------------------
# geometry helpers (local frame: axis = +z, board edge plane z = 0, connector in z > 0, board in z < 0,
# board top at local y = -z_axis)
# ------------------------------------------------------------------------------------------------------------
def cyl(z0, z1, r, x=0.0, y=0.0):
    return (3, o.addCylinder(x, y, z0, 0, 0, z1 - z0, r))


def cone(z0, z1, r0, r1):
    return (3, o.addCone(0, 0, z0, 0, 0, z1 - z0, r0, r1))


def box(x0, y0, z0, x1, y1, z1):
    return (3, o.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0))


def cut(a, tools, keep_tools=False):
    out, _ = o.cut([a] if isinstance(a, tuple) else a, tools if isinstance(tools, list) else [tools],
                   removeObject=True, removeTool=not keep_tools)
    return out


def fuse(objs):
    if len(objs) == 1:
        return objs
    out, _ = o.fuse(objs[:1], objs[1:])
    return out


def external_thread(r_major: float, pitch: float, z0: float, z1: float, segs: int = 24) -> list:
    """ISO/UN 60 deg external thread between z0 and z1 (flat crest P/8, root P/4), as one solid.

    A triangle profile is swept along a helix (Frenet frame) around the minor-diameter core and the result is
    clipped to [z0, z1] with 45 deg run-outs (chamfers) at both ends."""
    H = math.sqrt(3) / 2 * pitch
    r_min = r_major - 5 / 8 * H
    L = z1 - z0
    nturn = L / pitch + 2
    n = int(nturn * segs) + 1
    pts = [o.addPoint(r_min * math.cos(2 * math.pi * i / segs), r_min * math.sin(2 * math.pi * i / segs),
                      z0 - pitch + pitch * i / segs) for i in range(n)]
    wire = o.addWire([o.addBSpline(pts)])
    ov = 0.05 * pitch   # profile overlaps the core so the fuse is clean
    zc = z0 - pitch
    prof = [(r_min - ov, zc - pitch * 0.5 + 1e-3), (r_major, zc - pitch / 16), (r_major, zc + pitch / 16),
            (r_min - ov, zc + pitch * 0.5 - 1e-3)]
    pp = [o.addPoint(r, 0, z) for r, z in prof]
    ll = [o.addLine(pp[i], pp[(i + 1) % 4]) for i in range(4)]
    s = o.addPlaneSurface([o.addCurveLoop(ll)])
    sweep = o.addPipe([(2, s)], wire, "Frenet")
    core = cyl(z0 - 1.5 * pitch, z1 + 1.5 * pitch, r_min)
    th = fuse([core] + sweep)
    ch = 0.6 * (r_major - r_min) + 0.05
    clip = fuse([cone(z0, z0 + ch, r_major - ch, r_major), cyl(z0 + ch, z1 - ch, r_major + 0.01),
                 cone(z1 - ch, z1, r_major, r_major - ch)])
    clip = fuse([clip[0], cyl(z0 + ch * 0.5, z1 - ch * 0.5, r_min)])  # keep the core under the run-outs
    out, _ = o.intersect(th, clip)
    return out


def socket(z_front: float, a: float, bore: float, depth: float, z_back: float, nslot: int = 4) -> list:
    """Female centre contact: tube a/bore with nslot slots of 0.15 mm, entry chamfer."""
    body = cyl(z_back, z_front, a / 2)
    hole = cyl(z_front - depth, z_front + 0.01, bore / 2)
    ent = cone(z_front - 0.12, z_front + 0.001, bore / 2, bore / 2 + 0.1)
    body = cut(body, [hole, ent])
    slots = []
    for k in range(nslot // 2):
        s = box(-a, -0.075, z_front - depth * 0.8, a, 0.075, z_front + 0.01)
        o.rotate([s], 0, 0, 0, 0, 0, 1, math.pi / nslot * 2 * k)
        slots.append(s)
    return cut(body, slots)


# ------------------------------------------------------------------------------------------------------------
# templates
# ------------------------------------------------------------------------------------------------------------
def build_legs_tab(p: Params, d: dict) -> list[tuple[str, str, list]]:
    """End-launch jack with a flat tab centre contact and legs that straddle the board (Johnson 142-0701-8xx,
    Molex 73251-1150, Amphenol 132322, SV 1521-00002 / 1621-60050 style)."""
    b, a = d["b"], d["a"]
    Lb, Lt = d["body_len_est"], d["thread_len_est"]
    z_if = Lb + Lt                               # reference plane
    W = max(d["body_w_est"], 2 * (d["z_axis"] + p.board_h + d["leg_t_est"]) + 1.0)
    yb = -d["z_axis"]                            # board top (local y)
    parts = []
    # outer body: square block + threaded barrel, bore = outer conductor
    body = box(-W / 2, -W / 2, 0, W / 2, W / 2, Lb)
    barrel_core = cyl(Lb - 0.01, z_if - 0.25, (d["thr_D"] / 2) - 0.6)
    thr = external_thread(d["thr_D"] / 2, d["thr_P"], Lb + 0.3, z_if, p.thread_segments)
    outer = fuse([body, barrel_core] + thr)
    # legs
    legs = []
    xs = [-d["leg_x_est"], d["leg_x_est"]] if p.n_legs == 4 else [0.0]
    for x in xs:
        if p.n_legs == 2 and x == 0.0:
            x = d["leg_x_est"]
        legs.append(box(x - d["leg_w_est"] / 2, yb, -d["leg_len_est"], x + d["leg_w_est"] / 2, yb + d["leg_t_est"],
                        0.01))
        legs.append(box(x - d["leg_w_est"] / 2, yb - p.board_h - d["leg_t_est"], -d["leg_len_est"],
                        x + d["leg_w_est"] / 2, yb - p.board_h, 0.01))
    outer = fuse(outer + legs)
    bore = cyl(-0.1, z_if + 0.1, b / 2)
    outer = cut(outer, [bore])
    parts.append(("outer_body", "gold-plated brass", outer))
    # centre contact: socket at the front, tab at the back
    cc = socket(z_if - 0.05, a, d["socket_bore"], d["socket_depth_est"], Lb * 0.3)
    tab = box(-d["tab_w_est"] / 2, -d["tab_t_est"] / 2, -d["tab_len_est"], d["tab_w_est"] / 2, d["tab_t_est"] / 2,
              Lb * 0.3 + 0.01)
    neck = []
    if p.interface != "SMA":
        z0 = d["bead_z_est"] - d["bead_len_est"] / 2
        cc = cut(cc, [cyl(z_if - d["bead_z_est"] - d["bead_len_est"] / 2, z_if - d["bead_z_est"] + d["bead_len_est"] / 2,
                          a, )], keep_tools=False)
        neck = [cyl(z_if - d["bead_z_est"] - d["bead_len_est"] / 2 - 0.01,
                    z_if - d["bead_z_est"] + d["bead_len_est"] / 2 + 0.01, d["neck"] / 2)]
    cc = fuse(cc + [tab] + neck)
    parts.append(("center_contact", "gold-plated BeCu", cc))
    # dielectrics
    if p.interface == "SMA":
        diel = cut(cyl(0, z_if, b / 2), [cyl(-1, z_if + 1, a / 2)])
        parts.append(("dielectric", "PTFE", diel))
    else:
        zb0, zb1 = z_if - d["bead_z_est"] - d["bead_len_est"] / 2, z_if - d["bead_z_est"] + d["bead_len_est"] / 2
        bead = cut(cyl(zb0, zb1, b / 2), [cyl(zb0 - 1, zb1 + 1, d["neck"] / 2)])
        parts.append(("bead", p.bead, bead))
        rear = cut(cyl(0, Lb * 0.3 + 0.5, b / 2), [cyl(-1, Lb + 1, a / 2)])
        parts.append(("rear_insulator", "PTFE", rear))
    return parts


def build_flange_block(p: Params, d: dict) -> list[tuple[str, str, list]]:
    """Two-hole flange connector screwed into a transition block that clamps the board (Southwest Microwave
    end launch; Rosenberger, SV solderless, Samtec and HUBER+SUHNER edge mounts work the same way)."""
    b, a = d["b"], d["a"]
    yb = -d["z_axis"]
    Wb, Lb, Ht, Hb, Lw = d["block_w_est"], d["block_len_est"], d["block_top_est"], d["block_bot_est"], d["wall_est"]
    Lf, Lt = d["flange_t_est"], d["thread_len_est"]
    z_fl = Lw                          # flange rear face
    z_if = Lw + Lf + Lt                # reference plane
    ys = yb + Ht / 2                   # height of the two #0-80 connector screws (into the top part)
    parts = []
    # --- transition block (top part + bottom clamp plate), slot for the board
    top = box(-Wb / 2, yb + 0.0, -Lb, Wb / 2, yb + Ht, Lw)
    wall_low = box(-Wb / 2, yb - p.board_h - Hb, 0.0, Wb / 2, yb, Lw)     # front wall below the axis (ahead of edge)
    bot = box(-Wb / 2, yb - p.board_h - Hb, -Lb, Wb / 2, yb - p.board_h, 0.0)
    blk = fuse([top, wall_low, bot])
    chan = box(-d["channel_w_est"] / 2, yb - 0.01, -d["channel_len_est"], d["channel_w_est"] / 2,
               yb + d["channel_h_est"], 0.0)
    hole = cyl(-0.1, Lw + 0.1, d["launch_diel"] / 2)
    screws = []
    for x in (-d["screw_dx_est"], d["screw_dx_est"]):
        s = (3, o.addCylinder(x, yb - p.board_h - Hb - 1, -Lb / 2, 0, Ht + p.board_h + Hb + 2, 0, 0.93))  # #1-72
        screws.append(s)
        screws.append(cyl(Lw - 2.6, Lw + 0.1, 0.75, x, ys))                                           # #0-80 tap
    blk = cut(blk, [chan, hole] + screws)
    parts.append(("transition_block", "nickel-plated brass", blk))
    # --- launch dielectric and pin
    ld = cut(cyl(0.0, Lw, d["launch_diel"] / 2), [cyl(-1, Lw + 1, d["launch_pin"] / 2)])
    parts.append(("launch_dielectric", "PTFE", ld))
    # --- connector: flange + threaded barrel, bore = outer conductor
    fl = box(-Wb / 2, -d["flange_h_est"] / 2, z_fl, Wb / 2, d["flange_h_est"] / 2, z_fl + Lf)
    core = cyl(z_fl + Lf - 0.01, z_if - 0.25, d["thr_D"] / 2 - 0.6)
    thr = external_thread(d["thr_D"] / 2, d["thr_P"], z_fl + Lf + 0.8, z_if, p.thread_segments)
    housing = fuse([fl, core] + thr)
    holes = [cyl(z_fl - 0.1, z_fl + Lf + 0.1, 0.95, x, ys) for x in (-d["screw_dx_est"], d["screw_dx_est"])]
    cbore = [cyl(z_fl + Lf - 0.9, z_fl + Lf + 0.1, 1.3, x, ys) for x in (-d["screw_dx_est"], d["screw_dx_est"])]
    housing = cut(housing, [cyl(z_fl - 0.1, z_if + 0.1, b / 2)] + holes + cbore)
    parts.append(("connector_housing", "passivated stainless steel (CRES)", housing))
    # --- centre contact: socket, neck in the bead, step down to the launch pin
    zb0, zb1 = z_if - d["bead_z_est"] - d["bead_len_est"] / 2, z_if - d["bead_z_est"] + d["bead_len_est"] / 2
    cc = socket(z_if - 0.05, a, d["socket_bore"], d["socket_depth_est"], zb1)
    rear = cyl(z_fl + 0.3, zb0, a / 2)
    neck = cyl(zb0 - 0.01, zb1 + 0.01, d["neck"] / 2)
    taper = cone(z_fl, z_fl + 0.3 + 0.001, d["launch_pin"] / 2, a / 2)
    pin = cyl(-d["launch_over_est"], z_fl + 0.01, d["launch_pin"] / 2)
    cc = fuse(cc + [rear, neck, taper, pin])
    parts.append(("center_contact_and_launch_pin", "gold-plated BeCu", cc))
    # --- interface dielectric
    if p.interface == "SMA":
        parts.append(("dielectric", "PTFE", cut(cyl(zb1, z_if, b / 2), [cyl(zb1 - 1, z_if + 1, a / 2)])))
    parts.append(("bead", p.bead, cut(cyl(zb0, zb1, b / 2), [cyl(zb0 - 1, zb1 + 1, d["neck"] / 2)])))
    # --- screws (simplified heads + shank, no thread): 2x #0-80 connector screws, 2x #1-72 board screws
    sc = []
    for x in (-d["screw_dx_est"], d["screw_dx_est"]):
        sc += [cyl(z_fl + Lf - 0.85, z_fl + Lf + 0.65, 1.22, x, ys), cyl(z_fl - 2.4, z_fl + Lf - 0.85, 0.72, x, ys)]
        sc += [(3, o.addCylinder(x, yb + Ht, -Lb / 2, 0, 1.5, 0, 1.5)),
               (3, o.addCylinder(x, yb - p.board_h - Hb + 0.4, -Lb / 2, 0, Ht + p.board_h + Hb - 0.4, 0, 0.9))]
    parts.append(("screws", "stainless steel", fuse(sc)))
    return parts


TEMPLATES = {"legs_tab": build_legs_tab, "flange_block": build_flange_block}

MODELS = [
    Params("sma_legs_tab_jack_h1.57", "SMA", "legs_tab", board_h=1.575, n_legs=4,
           origin={"like": "Johnson 142-0701-801, Molex 73251-1150, Amphenol 132322"}),
    Params("k292_legs_tab_jack_h1.57", "2.92mm", "legs_tab", board_h=1.575,
           origin={"like": "SV 1521-00002 / 1521-60051, Cinch 145-0701-251"}),
    Params("k24_legs_tab_jack_h1.57", "2.4mm", "legs_tab", board_h=1.575,
           origin={"like": "SV 1621-60050 / 1621-60008"}),
    Params("sma_flange_block_jack_010_h0.508", "SMA", "flange_block", board_h=0.508, launch_option=".010",
           origin={"like": "Southwest Microwave 292-04A-12"}),
    Params("k292_flange_block_jack_010_h0.203", "2.92mm", "flange_block", board_h=0.203, launch_option=".010",
           origin={"like": "Southwest Microwave 1092-03A-12"}),
    Params("k292_flange_block_jack_005_h0.203", "2.92mm", "flange_block", board_h=0.203, launch_option=".005",
           origin={"like": "Southwest Microwave 1092-01A-12"}),
    Params("k24_flange_block_jack_007_h0.203", "2.4mm", "flange_block", board_h=0.203, launch_option=".007",
           origin={"like": "Southwest Microwave 1492-03A-12"}),
    Params("k24_flange_block_jack_005_h0.203", "2.4mm", "flange_block", board_h=0.203, launch_option=".005",
           origin={"like": "Southwest Microwave 1492-04A-12"}),
]


def build(p: Params, outdir: Path) -> dict:
    global o
    gmsh.model.add(p.name)
    o = gmsh.model.occ
    d = derived(p)
    parts = TEMPLATES[p.template](p, d)
    o.synchronize()
    # local -> Hertzon frame: rotate +90 deg about X (local z -> -y, local y -> z), lift to z_axis
    allv = [dt for _, _, dts in parts for dt in dts]
    o.rotate(allv, 0, 0, 0, 1, 0, 0, math.pi / 2)
    o.translate(allv, 0, 0, d["z_axis"])
    o.synchronize()
    bodies = []
    for name, mat, dts in parts:
        for k, (dim, tag) in enumerate(dts):
            nm = name if len(dts) == 1 else f"{name.rstrip('s')}_{k + 1}"
            gmsh.model.setEntityName(dim, tag, nm)
            hm = (DIELECTRICS.get(mat) or {}).get("hertzon_material")
            if hm is None:
                from connectors_db import METALS
                hm = (METALS.get(mat) or {}).get("hertzon_material", mat)
            bodies.append({"index": len(bodies), "name": nm, "material": mat, "hertzon_material": hm,
                           "volume_mm3": round(o.getMass(dim, tag), 4)})
    overlaps = []
    flat = [(nm, dt) for nm, _, dts in parts for dt in dts]
    for i in range(len(flat)):
        for j in range(i + 1, len(flat)):
            c1, c2 = o.copy([flat[i][1]]), o.copy([flat[j][1]])
            r, _ = o.intersect(c1, c2)
            v = sum(o.getMass(*q) for q in r if q[0] == 3)
            if r:
                o.remove(r, recursive=True)
            if v > 1e-6:
                overlaps.append((flat[i][0], flat[j][0], round(v, 6)))
    o.synchronize()
    step = outdir / f"{p.name}.step"
    gmsh.write(str(step))
    meta = {"model": p.name, "params": asdict(p), "dimensions_mm": {k: round(v, 4) for k, v in d.items()},
            "frame": "board top copper z=0 at t_cu, board in y>0, edge at y=0, axis along -Y at z=z_axis",
            "bodies": bodies, "step": step.name, "overlaps_mm3": overlaps,
            "estimated": sorted(k for k in d if k.endswith("_est"))}
    (outdir / f"{p.name}.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    gmsh.model.remove()
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default=str(ROOT / "step"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    for p in MODELS:
        if args.only and args.only not in p.name:
            continue
        m = build(p, out)
        print(p.name, "overlaps:", m["overlaps_mm3"], len(m["bodies"]), "bodies:", ", ".join(f"{b['name']}[{b['hertzon_material']}]" for b in m["bodies"]))
    gmsh.finalize()


if __name__ == "__main__":
    main()
