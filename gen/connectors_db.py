"""Edge-launch (end-launch) PCB connector database: SMA, 2.92 mm and 2.4 mm.

python gen/connectors_db.py  ->  data/interfaces.json, data/connectors.json, data/connectors.csv
python gen/connectors_db.py --hertzon <Hertzon checkout>  ->  also its hertzon/model/data/edge_launch_connectors.json

Every value carries its origin:
  "source"    : read on the manufacturer's page / datasheet / catalogue named in the part's `sources`
  "standard"  : interface dimension from IEEE 287 / MIL-STD-348 / IEC 61169 (same for every vendor)
  "estimate"  : not published (or only in a drawing that was not readable); typical value, check before use
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = 25.4  # mm per inch
C0 = 299792458.0


def z_coax(d_outer: float, d_inner: float, er: float) -> float:
    return 59.9585 / math.sqrt(er) * math.log(d_outer / d_inner)


def te11_cutoff_ghz(d_outer: float, d_inner: float, er: float) -> float:
    """First higher-order mode (TE11) of a coaxial line, approximation fc = c / (pi * (a + b) * sqrt(er))."""
    a, b = d_inner / 2e3, d_outer / 2e3
    return C0 / (math.pi * (a + b) * math.sqrt(er)) / 1e9


def d_inner_for_z(d_outer: float, er: float, z: float = 50.0) -> float:
    return d_outer / math.exp(z * math.sqrt(er) / 59.9585)


def ellipk(k: float) -> float:
    a, b = 1.0, math.sqrt(1 - k * k)
    for _ in range(40):
        a, b = (a + b) / 2, math.sqrt(a * b)
    return math.pi / (2 * a)


def kk(k: float) -> float:
    """K(k) / K(k')"""
    return ellipk(k) / ellipk(math.sqrt(1 - k * k))


def cpwg(w: float, g: float, h: float, er: float) -> tuple[float, float]:
    k = w / (w + 2 * g)
    k3 = math.tanh(math.pi * w / (4 * h)) / math.tanh(math.pi * (w + 2 * g) / (4 * h))
    r = kk(k3) / kk(k)        # = K(k3)/K(k3') * K(k')/K(k)
    eeff = (1 + er * r) / (1 + r)
    z = 60 * math.pi / math.sqrt(eeff) / (kk(k) + kk(k3))
    return z, eeff


def cpwg_width(z0: float, g: float, h: float, er: float) -> float:
    lo, hi = 0.02, 10.0
    for _ in range(80):
        m = (lo + hi) / 2
        if cpwg(m, g, h, er)[0] > z0:
            lo = m
        else:
            hi = m
    return (lo + hi) / 2



# Calibrated GCPW widths for the shielded (flange_block) Hertzon models: there the line runs in a closed box
# (via-fence walls, lid at Hair), and the analytic width gave 45.4 ohm (0.203 mm) / 45.35 ohm (0.508 mm) wave-port
# impedance. Two solves per board and linear interpolation give these widths for 50 ohm. The open board of the
# legs_tab models gives 50.0 ohm with the analytic width (checked on 1.575 mm), so it keeps cpwg_width.
LINE_CAL_SHIELDED: dict[float, float] = {0.203: 0.330, 0.508: 0.755}


def line_for_board(h: float, er: float = 3.55, shielded: bool = False, narrow: bool = False) -> tuple[float, float]:
    """50 ohm GCPW. narrow: smallest practical gap (0.1 mm) so the line stays narrow on thick boards at 40-50 GHz."""
    g = 0.1 if narrow else (0.15 if h <= 0.3 else 0.25)
    w = (LINE_CAL_SHIELDED.get(round(h, 3)) if shielded else None) or cpwg_width(50.0, g, h, er)
    return w, g


# ------------------------------------------------------------------------------------------------------------
# Dielectrics used in these connectors (values at microwave frequencies, typical; vendors rarely publish them)
# hertzon_material: name in Hertzon's System library (hertzon/model/materials.py) or a new one to add.
# ------------------------------------------------------------------------------------------------------------
DIELECTRICS = {
    "PTFE": {"er": 2.05, "tand": 0.0003, "hertzon_material": "Teflon (tm)",
             "note": "Hertzon 'Teflon (tm)' has er 2.1 / tand 0.001; solid PTFE in connectors is 2.02-2.08"},
    "Ultem 1000 (PEI)": {"er": 3.15, "tand": 0.0025, "hertzon_material": "Ultem_1000",
                         "note": "new material, er 3.15 tand 0.0025 at 1-10 GHz (SABIC data, typical)"},
    "Kel-F (PCTFE)": {"er": 2.35, "tand": 0.005, "hertzon_material": "Kel-F_PCTFE",
                      "note": "new material, er 2.3-2.5, tand ~0.005 at GHz (typical)"},
    "PEEK": {"er": 3.2, "tand": 0.003, "hertzon_material": "PEEK"},
    "glass (Corning 7070 type)": {"er": 4.1, "tand": 0.0012, "hertzon_material": "Glass_7070",
                                  "note": "hermetic glass beads, er 4.1 (7070) to 6.7 (Kovar seal glass)"},
    "air": {"er": 1.0006, "tand": 0.0, "hertzon_material": "air"},
}
METALS = {
    "gold-plated BeCu": {"sigma": 4.1e7, "hertzon_material": "gold",
                         "note": "skin depth at 10 GHz 0.79 um; 1.27 um Au is > 1.6 skin depths -> use gold"},
    "gold-plated brass": {"sigma": 4.1e7, "hertzon_material": "gold"},
    "passivated stainless steel (CRES)": {"sigma": 1.1e6, "hertzon_material": "stainless_steel"},
    "nickel-plated brass": {"sigma": 1.45e7, "hertzon_material": "nickel",
                            "note": "Ni is ferromagnetic; for HF use sigma_eff ~ 2-4e6 S/m (mu_r ~ 1 at >10 GHz but"
                                    " Ni losses are high); 'brass' underestimates the loss"},
    "brass": {"sigma": 1.5e7, "hertzon_material": "brass"},
}

# ------------------------------------------------------------------------------------------------------------
# Interfaces (mating side), dimensions in mm
# ------------------------------------------------------------------------------------------------------------
def _iface(name, std, fmax, b, a, er_name, thread, thread_major, thread_pitch, male_pin, notes):
    er = DIELECTRICS[er_name]["er"]
    return {
        "interface": name, "standard": std, "f_max_GHz": fmax,
        "outer_conductor_id_mm": b, "center_conductor_od_mm": a, "dielectric": er_name,
        "z0_ohm": round(z_coax(b, a, er), 2), "te11_cutoff_GHz": round(te11_cutoff_ghz(b, a, er), 1),
        "coupling_thread": thread, "thread_major_mm": thread_major, "thread_pitch_mm": round(thread_pitch, 4),
        "male_pin_od_mm": male_pin, "origin": "standard", "notes": notes,
    }


INTERFACES = {
    "SMA": _iface("SMA", "MIL-STD-348 Fig. 310, IEC 61169-15", 18.0, 4.13, 1.27, "PTFE", "1/4-36 UNS-2A (jack)",
                  6.35, IN / 36, 0.92,
                  "PTFE-filled at the interface. 18 GHz standard, precision/'Super SMA' 26.5-27 GHz. Dielectric "
                  "OD 0.161-0.165 in (4.09-4.19 mm); pin 0.036-0.037 in."),
    "2.92mm": _iface("2.92mm", "IEEE 287 / IEC 61169-35 (K, RPC-2.92, SK)", 40.0, 2.92, 1.27, "air",
                     "1/4-36 UNS-2A (jack)", 6.35, IN / 36, 0.91,
                     "Air at the interface; mechanically compatible with SMA and 3.5 mm. Support bead (Ultem, "
                     "Kel-F or PEEK) a few mm behind the reference plane, compensated by an undercut."),
    "2.4mm": _iface("2.4mm", "IEEE 287 / IEC 61169-32 (RPC-2.40, PC 2.4)", 50.0, 2.40, 1.042, "air",
                    "M7x0.75-6g (jack)", 7.0, 0.75, 0.51,
                    "Air at the interface; intermates with 1.85 mm. Not compatible with SMA/2.92 mm."),
}
NOTE_25 = ("There is no standard 2.5 mm coaxial connector. The user's '2.5 mm' is read as 2.4 mm (50 GHz). "
           "Nearby sizes: 2.92 mm (40 GHz), 2.4 mm (50 GHz), 1.85 mm (67 GHz, intermates with 2.4 mm).")

# Bead compensation estimates: the pin is necked inside the bead so the bead section stays 50 ohm.
for key in ("2.92mm", "2.4mm"):
    f = INTERFACES[key]
    f["bead_compensation_mm"] = {
        bead: {"neck_od_mm": round(d_inner_for_z(f["outer_conductor_id_mm"], DIELECTRICS[bead]["er"]), 3),
               "origin": "estimate (50 ohm with the outer conductor unchanged)"}
        for bead in ("Ultem 1000 (PEI)", "Kel-F (PCTFE)", "PEEK")}

# Southwest Microwave launch pin / dielectric options (inches): A = pin dia, C = dielectric dia
SWM_LAUNCH = {".010": (0.010, 0.0635), ".007": (0.007, 0.039), ".005": (0.005, 0.029)}
SWM_BLOCKS = {
    "-12": ("standard block, thread-in", 0.0, 0.065),
    "-7": ("standard block", 0.090, 0.210),
    "-8": ("standard block", 0.190, 0.300),
    "-9": ("narrow block (.350 in), thread-in", 0.0, 0.065),
    "-14": ("narrow block, thread-in", 0.040, 0.090),
}
SWM_SERIES = {
    # interface: (female base, male base, fmax, {pin option: (female dash, male dash)})
    "SMA": ("292", "293", 27.0, {".010": ("04A", "01A"), ".007": ("06A", "03A"), ".005": ("07A", "04A")}),
    "2.92mm": ("1092", "1093", 40.0, {".010": ("03A", "01A"), ".007": ("04A", "03A"), ".005": ("01A", "04A")}),
    "2.4mm": ("1492", "1493", 50.0, {".010": ("02A", "01A"), ".007": ("03A", "03A"), ".005": ("04A", "04A")}),
}
SWM_DS = "https://mpd.southwestmicrowave.com/wp-content/uploads/End-Launch-Connectors-2024.pdf"


def swm_parts() -> list[dict]:
    parts = []
    for iface, (fb, mb, fmax, pins) in SWM_SERIES.items():
        bead = {"SMA": "Ultem 1000 (PEI)", "2.92mm": "Ultem 1000 (PEI) / Kel-F (PCTFE)",
                "2.4mm": "Ultem 1000 (PEI) / Kel-F (PCTFE)"}[iface]
        for pin, (fd, md) in pins.items():
            dp, dd = SWM_LAUNCH[pin]
            for suf, (style, tmin, tmax) in SWM_BLOCKS.items():
                for gender, base, dash in (("jack (female)", fb, fd), ("plug (male)", mb, md)):
                    pn = f"{base}-{dash}{suf}"
                    if iface == "2.4mm" and gender.startswith("plug") and pin == ".010" and suf in ("-9", "-14"):
                        pn = f"{base}-02A{suf}"   # catalogue lists 1493-02A-9/-14 for the .010 pin
                    parts.append({
                        "manufacturer": "Southwest Microwave", "part_number": pn, "series": "End Launch",
                        "interface": iface, "gender": gender, "f_max_GHz": fmax,
                        "mounting": f"edge clamp, solderless ({style})", "template": "flange_block",
                        "pcb_thickness_mm": [round(tmin * IN, 3), round(tmax * IN, 3)],
                        "launch_pin_od_mm": round(dp * IN, 4), "launch_dielectric_od_mm": round(dd * IN, 4),
                        "materials": {
                            "housing": "passivated stainless steel (CRES)", "center_contact": "gold-plated BeCu",
                            "bead": bead, "interface_dielectric": "PTFE" if iface == "SMA" else "air",
                            "block": "nickel-plated brass", "launch_pin": "gold-plated BeCu",
                            "launch_dielectric": "PTFE", "screws": "#0-80 SHCS (connector), #1-72 SHCS (board)"},
                        "datasheet_url": SWM_DS, "product_url": "https://mpd.southwestmicrowave.com/",
                        "cad": {"step": "on request (catalogue: '3D models for mechanical layouts')",
                                "em_model": "encrypted HFSS model on request (HFSS 18+)", "step_url": None},
                        "origin": "source (catalogue 2024); body dimensions only in drawings (images) -> estimate",
                        "notes": "Super SMA (27 GHz)" if iface == "SMA" else "",
                    })
    return parts


def P(manufacturer, pn, iface, gender, fmax, mounting, template, pcb=None, materials=None, ds=None, url=None,
      cad=None, notes="", origin="source"):
    return {"manufacturer": manufacturer, "part_number": pn, "series": "", "interface": iface, "gender": gender,
            "f_max_GHz": fmax, "mounting": mounting, "template": template, "pcb_thickness_mm": pcb,
            "launch_pin_od_mm": None, "launch_dielectric_od_mm": None, "materials": materials or {},
            "datasheet_url": ds, "product_url": url,
            "cad": cad or {"step": "unknown", "step_url": None, "em_model": None}, "origin": origin, "notes": notes}


SV = "https://www.svmicrowave.com/Attachment/DownloadFile?downloadId="
JOHNSON_DS = "https://www.eg.bucknell.edu/~dkelley/eceg497/lab/Datasheet_SMA-BoardEdge_Johnson.pdf"
JOHNSON_MAT = {"body": "gold-plated brass (or nickel)", "center_contact": "gold-plated BeCu (female)",
               "interface_dielectric": "PTFE", "nut_spring": "BeCu"}
HS_DS = "https://shop.richardsonrfpd.com/docs/rfpd/Precision_compression_mount_PCB_connectors.pdf"
HS_MAT = {"body": "passivated stainless steel (CRES)", "outer_contact": "passivated stainless steel (CRES)",
          "center_contact": "gold-plated BeCu", "bead": "Ultem 1000 (PEI)"}
SAMTEC_MAT = {"body": "passivated stainless steel (CRES)", "center_contact": "gold-plated BeCu (50 uin Au over Ni)",
              "bead": "Ultem 1000 (PEI)", "launch_dielectric": "PTFE", "screws": "stainless steel"}

OTHER_PARTS = [
    # --- Amphenol SV Microwave ---------------------------------------------------------------------------
    P("Amphenol SV Microwave", "1521-00002", "2.92mm", "jack (female)", 40, "edge launch, 4 legs, solder",
      "legs_tab", materials={"finish": "gold"}, ds=SV + "8355",
      url="https://www.svmicrowave.com/292mm-female-pcb-edge-launch-connector-4-leg",
      cad={"step": "yes (direct download)", "step_url": SV + "8356",
           "em_model": f"HFSS {SV}11843 ; CST {SV}14754"}),
    P("Amphenol SV Microwave", "1621-60050", "2.4mm", "jack (female)", 50, "edge launch, 4 legs, solder",
      "legs_tab", materials={"finish": "gold"}, ds=SV + "8374",
      url="https://www.svmicrowave.com/24mm-female-pcb-edge-launch-connector-4-legs",
      cad={"step": "yes (direct download)", "step_url": SV + "8375",
           "em_model": f"HFSS {SV}22688 ; CST {SV}14759"}),
    P("Amphenol SV Microwave", "1521-60051", "2.92mm", "jack (female)", 40, "edge launch, solder", "legs_tab",
      pcb=[1.575, 1.575], url="https://www.svmicrowave.com/292mm-female-pcb-edge-launch-connector-062-pcb-thickness",
      cad={"step": "likely (SV product pages offer CAD download)", "step_url": None, "em_model": "HFSS/CST likely"}),
    P("Amphenol SV Microwave", "1521-60125", "2.92mm", "jack (female)", 40, "edge launch, solder", "legs_tab",
      pcb=[2.362, 2.362], url="https://www.svmicrowave.com/292mm-female-pcb-edge-launch-connector-093-pcb-thickness",
      cad={"step": "likely (SV product pages offer CAD download)", "step_url": None, "em_model": "HFSS/CST likely"}),
    P("Amphenol SV Microwave", "1621-60008", "2.4mm", "jack (female)", 50, "edge launch, solder", "legs_tab",
      pcb=[1.575, 1.575], url="https://www.svmicrowave.com/24mm-female-pcb-edge-launch-connector-062-pcb-thickness",
      cad={"step": "likely (SV product pages offer CAD download)", "step_url": None, "em_model": "HFSS/CST likely"}),
    P("Amphenol SV Microwave", "1621-60051", "2.4mm", "jack (female)", 50, "edge launch, solder", "legs_tab",
      pcb=[2.362, 2.362], url="https://www.svmicrowave.com/24mm-female-pcb-edge-launch-connector-093-pcb-thickness",
      cad={"step": "likely (SV product pages offer CAD download)", "step_url": None, "em_model": "HFSS/CST likely"}),
    P("Amphenol SV Microwave", "1621-60092", "2.4mm", "jack (female)", 50,
      "edge launch, 2-hole flange, torque resistant (TraceGuard), solder", "flange_block",
      materials={"finish": "gold", "seal": "glass"},
      ds="https://www.svmicrowave.com/Attachment/DownloadFile?downloadId=16929&productId=6513",
      url="https://svmicrowave.com/24mm-female-pcb-edge-launch-connector-2-hole-torque-resistant",
      cad={"step": "yes (direct download)",
           "step_url": "https://www.svmicrowave.com/Attachment/DownloadFile?downloadId=22711&productId=6513",
           "em_model": "HFSS/CST on request form"}),
    P("Amphenol SV Microwave", "2921-61778", "SMA", "jack (female)", None, "edge launch, TraceGuard, solder",
      "flange_block", ds=SV + "16927", notes="datasheet not opened; f_max not verified"),
    P("Amphenol SV Microwave", "1521-60235", "2.92mm", "jack (female)", 40, "edge launch, TraceGuard, solder",
      "flange_block", ds=SV + "16925", notes="datasheet not opened"),
    P("Amphenol SV Microwave", "SF2921-61450", "SMA", "jack (female)", 18, "edge, solderless compression, 2-hole",
      "flange_block", ds="https://www.mouser.com/datasheet/3/59/1/SF2921-61450_PDD.pdf",
      url="https://www.mouser.in/new/amphenol/amphenol-sv-solderless-edge-launch"),
    P("Amphenol SV Microwave", "SF1521-60107", "2.92mm", "jack (female)", 40,
      "edge, solderless compression, 2-hole", "flange_block", ds="https://svmicrowave.com/Attachment/DownloadFile?downloadId=8682",
      url="https://www.mouser.in/new/amphenol/amphenol-sv-solderless-edge-launch", notes="usable to 46 GHz (Mouser)"),
    P("Amphenol SV Microwave", "SF1621-60036", "2.4mm", "jack (female)", 50,
      "edge, solderless compression, 2-hole", "flange_block",
      ds="https://svmicrowave.com/Attachment/DownloadFile?downloadId=8725",
      url="https://www.mouser.in/new/amphenol/amphenol-sv-solderless-edge-launch"),
    # --- Rosenberger -------------------------------------------------------------------------------------
    P("Rosenberger", "02K243-40ME3", "2.92mm", "jack (female)", 40,
      "edge launch, clamp (1 clamp, 2 screws), center contact soldered", "flange_block",
      materials={"outer_contact": "brass, >=0.8 um Au over chemical Ni", "center_contact": "BeCu, >=1.27 um Au over Ni",
                 "dielectric": "PEEK (dielectric 1) + PTFE (dielectric 2)", "screws": "stainless steel"},
      ds="https://products.rosenberger.com/_ocassets/db/02K243-40ME3.pdf",
      url="https://www.rfmw.com/products/detail/02k24340me3-rosenberger/241019/",
      cad={"step": "not found on public pages (ask Rosenberger)", "step_url": None, "em_model": None},
      notes="RL >= 14 dB DC-40 GHz; PCB layout B 208"),
    P("Rosenberger", "32K243-40ML5", "SMA", "jack (female)", 18, "edge mount, clamp style", "flange_block",
      materials={"outer_contact": "brass, AuroDur gold", "center_contact": "BeCu, AuroDur gold",
                 "dielectric": "PTFE"},
      ds="https://products.rosenberger.com/_ocassets/db/32K243-40ML5.pdf",
      url="https://www.rfmw.com/products/detail/32k24340ml5-rosenberger/306965/",
      cad={"step": "not found on public pages (ask Rosenberger)", "step_url": None, "em_model": None},
      notes="VSWR 1.1 + 0.02 f[GHz]; PCB layout B 208"),
    # --- Samtec ------------------------------------------------------------------------------------------
    P("Samtec", "292-J-P-EP-ST-EL(W)-01", "2.92mm", "jack (female)", 40, "edge launch, solderless compression",
      "flange_block", materials=SAMTEC_MAT, ds="https://suddendocs.samtec.com/productspecs/292-elw.pdf",
      url="https://www.samtec.com/products/292-EL",
      cad={"step": "yes (Samtec configurator 3D/CAD download, login)", "step_url": None,
           "em_model": "Samtec offers simulation models on request"},
      notes="VSWR 1.10 typ to 18 GHz, 1.15 typ 18-40 GHz; IL 0.035*sqrt(f[GHz]) dB max"),
    P("Samtec", "240-EL series", "2.4mm", "jack (female)", 50, "edge launch, solderless compression",
      "flange_block", url="https://www.samtec.com/products/240-EL",
      cad={"step": "yes (Samtec configurator 3D/CAD download, login)", "step_url": None, "em_model": None}),
    # --- HUBER+SUHNER -----------------------------------------------------------------------------------
    P("HUBER+SUHNER", "85221629 (92_SMA-S50-0-200/122_NE)", "SMA", "jack (female)", 26.5,
      "edge mount, microstrip", "flange_block", materials=HS_MAT, ds=HS_DS,
      cad={"step": "3D libraries mentioned; HFSS model on request", "step_url": None, "em_model": "HFSS on request"},
      notes="RL >= 17 dB to 26.5 GHz"),
    P("HUBER+SUHNER", "85221630 (92_SK-S50-0-200/199_NE)", "2.92mm", "jack (female)", 40,
      "edge mount, microstrip", "flange_block", materials=HS_MAT, ds=HS_DS,
      cad={"step": "3D libraries mentioned; HFSS model on request", "step_url": None, "em_model": "HFSS on request"},
      notes="RL >= 17 dB to 40 GHz"),
    P("HUBER+SUHNER", "85221631 (92_PC24-S50-0-200/199_NE)", "2.4mm", "jack (female)", 50,
      "edge mount, microstrip", "flange_block", materials=HS_MAT, ds=HS_DS,
      cad={"step": "3D libraries mentioned; HFSS model on request", "step_url": None, "em_model": "HFSS on request"},
      notes="RL >= 17 dB to 50 GHz"),
    # --- Cinch / Johnson ---------------------------------------------------------------------------------
    *[P("Cinch Connectivity (Johnson)", pn, "SMA", "jack (female)", 18, "end launch, tab contact, 2 legs, solder",
        "legs_tab", pcb=[t, t], materials=JOHNSON_MAT, ds=JOHNSON_DS,
        url=f"https://www.cinch.com/products/rf-microwave/connectors/{pn}",
        cad={"step": "Cinch/Bel product page (CAD tab), not verified", "step_url": None, "em_model": None},
        notes="VSWR typ < 1.5 to 18 GHz")
      for pn, t in [("142-0701-801", 1.575), ("142-0701-806", 1.575), ("142-0701-851", 1.575),
                    ("142-0701-861", 1.575), ("142-0701-871", 1.575), ("142-0701-876", 1.575),
                    ("142-0711-821", 1.575), ("142-0711-826", 1.575), ("142-0711-871", 1.194),
                    ("142-0711-876", 1.194), ("142-0711-881", 1.194), ("142-0711-886", 1.194),
                    ("142-0701-881", 0.787), ("142-0701-886", 0.787)]],
    P("Cinch Connectivity (Johnson)", "145-0701-251", "2.92mm", "jack (female)", None, "end launch, solder",
      "legs_tab", url="https://www.digikey.at/en/products/detail/cinch-connectivity-solutions-johnson/145-0701-251/17278652",
      notes="f_max not read (2.92 mm: typ. 40 GHz)"),
    # --- Amphenol RF --------------------------------------------------------------------------------------
    P("Amphenol RF", "132322", "SMA", "jack (female)", 18, "end launch, solder", "legs_tab", pcb=[0, 1.575],
      materials={"body": "gold-plated brass", "center_contact": "gold-plated BeCu", "interface_dielectric": "PTFE"},
      url="https://www.amphenolrf.com/132322.html", cad={"step": "Amphenol RF product page (not readable here)",
                                                         "step_url": None, "em_model": None}),
    P("Amphenol RF", "132357", "SMA", "jack (female)", 18, "end launch, flat tab, ribbed, solder", "legs_tab",
      pcb=[1.575, 1.575], materials={"center_contact": "gold-plated BeCu", "body_finish": "gold"},
      ds="https://docs.rs-online.com/9216/A700000006757212.pdf", url="https://www.amphenolrf.com/132357.html"),
    P("Amphenol RF", "132255", "SMA", "jack (female)", None, "end launch, high frequency, 0.062 in", "legs_tab",
      pcb=[1.575, 1.575], url="https://www.amphenolrf.com/132255.html", notes="'High Frequency' variant; f_max not read"),
    P("Amphenol RF", "132289-14", "SMA", "jack (female)", None, "end launch, bulkhead, O-ring, 0.062 in", "legs_tab",
      pcb=[1.575, 1.575], url="https://www.amphenolrf.com/132289-14.html"),
    P("Amphenol RF", "132366-10", "SMA", "plug (male)", None, "end launch, 0.062 in", "legs_tab", pcb=[1.575, 1.575],
      url="https://www.amphenolrf.com/132366-10.html"),
    # --- Molex --------------------------------------------------------------------------------------------
    P("Molex", "73251-1150 (0732511150)", "SMA", "jack (female)", 18, "edge mount, solder", "legs_tab",
      pcb=[1.6, 1.6], materials={"body": "brass, 0.25 um Au", "center_contact": "gold-plated BeCu",
                                 "interface_dielectric": "PTFE"},
      ds="https://www.molex.com/en-us/products/part-detail/732511150?display=pdf",
      url="https://www.molex.com/en-us/products/part-detail/732511150",
      cad={"step": "yes (direct download)",
           "step_url": "https://www.molex.com/content/dam/molex/molex-dot-com/products/automated/en-us/3dcadmodels/732/73251/732511150_stp.zip",
           "em_model": None}),
    P("Molex", "73251-4432 (0732514432)", "SMA", "jack (female)", None, "board edge, end launch, solder", "legs_tab",
      url="https://www.digikey.com/en/products/detail/molex/0732514432/11305622"),
]


HERTZON_FILE = Path("hertzon", "model", "data", "edge_launch_connectors.json")


def hertzon_library(parts: list[dict]) -> dict:
    """The library as Hertzon ships it (hertzon/model/connectors.py reads it): interfaces, materials and the
    parts in one file."""
    return {"format": "hertzon-connectors", "version": 1,
            "source": "https://github.com/c8501024/hertzon-connectors (generated by gen/connectors_db.py)",
            "note_2_5mm": NOTE_25, "interfaces": INTERFACES, "dielectrics": DIELECTRICS, "metals": METALS,
            "swm_launch_options_in": SWM_LAUNCH,
            "parts": [{"id": p["id"], **{k: v for k, v in p.items() if k not in ("id", "series")}} for p in parts]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--hertzon", help="path of a Hertzon checkout: also write its connector library")
    args = ap.parse_args()
    parts = swm_parts() + OTHER_PARTS
    for i, p in enumerate(parts):
        p["id"] = f"{p['manufacturer'].split()[0].lower().replace('+', '')}-{p['part_number'].split()[0].lower()}"
    out = ROOT / "data"
    out.mkdir(exist_ok=True)
    (out / "interfaces.json").write_text(json.dumps(
        {"interfaces": INTERFACES, "dielectrics": DIELECTRICS, "metals": METALS, "note_2_5mm": NOTE_25,
         "swm_launch_options_in": SWM_LAUNCH}, indent=1, ensure_ascii=False), encoding="utf-8")
    (out / "connectors.json").write_text(json.dumps({"version": 1, "parts": parts}, indent=1, ensure_ascii=False),
                                         encoding="utf-8")
    cols = ["id", "manufacturer", "part_number", "interface", "gender", "f_max_GHz", "mounting", "template",
            "pcb_thickness_mm", "launch_pin_od_mm", "launch_dielectric_od_mm", "datasheet_url", "product_url",
            "step", "step_url", "em_model", "materials", "origin", "notes"]
    with open(out / "connectors.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cols)
        for p in parts:
            row = dict(p, **{k: p["cad"].get(k) for k in ("step", "step_url", "em_model")})
            row["materials"] = "; ".join(f"{k}: {v}" for k, v in p["materials"].items())
            row["pcb_thickness_mm"] = "" if not p["pcb_thickness_mm"] else "-".join(map(str, p["pcb_thickness_mm"]))
            w.writerow(["" if row.get(c) is None else row.get(c) for c in cols])
    if args.hertzon:
        f = Path(args.hertzon) / HERTZON_FILE
        f.write_text(json.dumps(hertzon_library(parts), indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"wrote {f}")
    by = {}
    for p in parts:
        by[p["interface"]] = by.get(p["interface"], 0) + 1
    print(len(parts), "parts", by)
    for k, f in INTERFACES.items():
        print(k, f["z0_ohm"], "ohm, TE11", f["te11_cutoff_GHz"], "GHz", f.get("bead_compensation_mm", ""))
    for pin, (dp, dd) in SWM_LAUNCH.items():
        print("SWM launch", pin, "Z(PTFE) =", round(z_coax(dd, dp, 2.05), 1), "ohm; er for 50 ohm =",
              round((59.9585 * math.log(dd / dp) / 50) ** 2, 2))


if __name__ == "__main__":
    main()
