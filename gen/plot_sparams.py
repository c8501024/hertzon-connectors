"""Plot |S11| and |S21| of the solved HF models: sparams/*.s2p -> sparams/overview.png, summary table on stdout."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def read_s2p(path: Path):
    f, s11, s21 = [], [], []
    for line in path.read_text().splitlines():
        if not line.strip() or line[0] in "!#":
            continue
        v = list(map(float, line.split()))
        f.append(v[0] / 1e9)
        s11.append(complex(v[1], v[2]))
        s21.append(complex(v[3], v[4]))
    return f, s11, s21


def db(x):
    return 20 * math.log10(max(abs(x), 1e-12))


def main():
    files = sorted((ROOT / "sparams").glob("*.s2p"))
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=False)
    rows = []
    for p in files:
        f, s11, s21 = read_s2p(p)
        lab = p.stem
        axes[0].plot(f, [db(x) for x in s11], label=lab)
        axes[1].plot(f, [db(x) for x in s21], label=lab)
        w = max(range(len(f)), key=lambda i: db(s11[i]))
        rows.append((lab, f[-1], db(s11[w]), f[w], min(db(x) for x in s21)))
    axes[0].set_ylabel("|S11| (dB)")
    axes[1].set_ylabel("|S21| (dB)")
    for a in axes:
        a.set_xlabel("Frequency (GHz)")
        a.grid(True, alpha=0.3)
    axes[0].axhline(-20, color="k", lw=0.8, ls="--")
    axes[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Hertzon (EMerge FEM): simplified edge-launch models, connector + launch + 50 ohm GCPW")
    fig.tight_layout()
    fig.savefig(ROOT / "sparams" / "overview.png", dpi=130)
    print("| Modell | f_max | schlechtestes S11 | bei | min S21 |\n|---|---|---|---|---|")
    for lab, fm, w11, fw, m21 in rows:
        print(f"| {lab} | {fm:g} GHz | {w11:.1f} dB | {fw:.1f} GHz | {m21:.2f} dB |")


if __name__ == "__main__":
    sys.exit(main())
