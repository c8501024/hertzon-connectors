"""Preview images of the STEP models: full view and a half section through the axis (x < 0 removed).

python gen/render.py  ->  preview/<model>.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import gmsh
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
COLORS = {"gold": "#e5c100", "Teflon (tm)": "#f4f4f4", "Ultem_1000": "#c8732a", "Kel-F_PCTFE": "#d9e6f2",
          "stainless_steel": "#9aa3aa", "stainless steel": "#8a939a", "nickel": "#b9b9b0", "PEEK": "#d8c8a8",
          "Rogers RO4003C (tm)": "#5f7f3f", "copper": "#d98c4a", "air": "#bfe0ff", "vacuum": "#bfe0ff"}


def body_meshes(step: Path, h: float = 0.08):
    gmsh.model.add(step.stem)
    gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
    out = gmsh.model.occ.importShapes(str(step), highestDimOnly=True)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.model.mesh.generate(2)
    meshes = []
    for dim, tag in out:
        tris = []
        for _, s in gmsh.model.getBoundary([(dim, tag)], oriented=False):
            types, _, nodes = gmsh.model.mesh.getElements(2, abs(s))
            for t, nd in zip(types, nodes):
                if t == 2:
                    tris.append(nd.reshape(-1, 3))
        ntags, xyz, _ = gmsh.model.mesh.getNodes()
        idx = {int(t): i for i, t in enumerate(ntags)}
        P = xyz.reshape(-1, 3)
        T = np.vectorize(idx.get)(np.vstack(tris).astype(int))
        faces = np.hstack([np.full((len(T), 1), 3), T]).ravel()
        meshes.append(pv.PolyData(P, faces).clean())
    gmsh.model.remove()
    return meshes


def render(meta_file: Path, out: Path):
    meta = json.loads(meta_file.read_text())
    meshes = body_meshes(meta_file.with_suffix(".step"))
    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1800, 800), border=False)
    for col, section in ((0, False), (1, True)):
        pl.subplot(0, col)
        for b, m in zip(meta["bodies"], meshes):
            if section:
                m = m.clip(normal=(-1, 0, 0), origin=(0, 0, 0), invert=False)
                if m.n_points == 0:
                    continue
            pl.add_mesh(m, color=COLORS.get(b["hertzon_material"], "#cccccc"), smooth_shading=True,
                        specular=0.4, show_edges=False)
        pl.add_text(("half section (x>0)" if section else meta["model"]), font_size=10)
        pl.camera_position = [(28, -34, 22), (0, -6, 0), (0, 0, 1)] if not section else [(48, -4, 4), (0, -4, 0), (0, 0, 1)]
        pl.camera.zoom(1.2 if not section else 1.25)
    pl.set_background("white")
    pl.screenshot(str(out))
    pl.close()


def main():
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    outdir = ROOT / "preview"
    outdir.mkdir(exist_ok=True)
    for mf in sorted((ROOT / "step").glob("*.json")):
        if len(sys.argv) > 1 and sys.argv[1] not in mf.stem:
            continue
        render(mf, outdir / f"{mf.stem}.png")
        print(mf.stem)
    gmsh.finalize()


if __name__ == "__main__":
    main()
