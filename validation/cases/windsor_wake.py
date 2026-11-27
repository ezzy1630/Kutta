"""Compare Kutta at the measured no-wheel Set B horizontal PIV vectors."""

import csv
import hashlib
import math

PIV_FILE = "SquareBack_NW_2_5yaw_Z=0194m_FlowField.csv"


def compare_plane(directory, drag, side, speed, export, view, image=None):
    path = directory / PIV_FILE
    samples = []
    missing = 0
    with path.open(newline="") as src:
        for row in csv.DictReader(src):
            u, v = (float(row[f"{a}_mean (m/s)"]) for a in "uv")
            if -9999 in (u, v):
                missing += 1
                continue
            x, y, z = (float(row[f"{a} (m)"]) for a in "xyz")
            samples.append(
                ([x * drag[0] + y * side[0], x * drag[1] + y * side[1], z], u, v)
            )
    errors = [[], [], []]
    chart = []
    measured_reverse = predicted_reverse = missed_reverse = outside = 0
    for i in range(0, len(samples), 12000):
        rows = samples[i : i + 12000]
        coords = [[round(c, 7) for c in xyz] for xyz, _, _ in rows]
        export(
            {
                "probes": [
                    {
                        "type": "points",
                        "name": "Planar PIV",
                        "points": coords,
                        "fields": ["velocity_x", "velocity_y"],
                    }
                ]
            }
        )
        with (view / "probe-0.csv").open(newline="") as src:
            for actual, (xyz, u, v) in zip(csv.DictReader(src), rows, strict=True):
                if actual["valid"] != "1":
                    outside += 1
                    continue
                if any(
                    abs(float(actual[k]) - xyz[j]) > 1e-5 for j, k in enumerate("xyz")
                ):
                    raise ValueError("PIV and CFD coordinates are misaligned")
                vx, vy = (float(actual[f"velocity_{a}"]) / speed for a in "xy")
                predicted_u = vx * drag[0] + vy * drag[1]
                predicted_v = vx * side[0] + vy * side[1]
                errors[0].append(predicted_u - u / 40)
                errors[1].append(predicted_v - v / 40)
                errors[2].append(math.hypot(vx, vy) - math.hypot(u, v) / 40)
                if image:
                    x = xyz[0] * drag[0] + xyz[1] * drag[1]
                    y = xyz[0] * side[0] + xyz[1] * side[1]
                    chart.append((x, y, u / 40, predicted_u))
                measured_reverse += u < 0
                predicted_reverse += predicted_u < 0
                missed_reverse += (u < 0) != (predicted_u < 0)
    if not errors[0]:
        raise ValueError("No matched horizontal PIV samples")
    if image:
        plot_plane(chart, image, samples[0][0][2])
    n = len(errors[0])
    return {
        "source_file": PIV_FILE,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "piv_nominal_speed_m_s": 40,
        "CFD_reference_m_s": speed,
        "matched": n,
        "invalid_experimental": missing,
        "valid_experimental_outside_CFD": outside,
        "axes": "body",
        "normalized_components": {
            key: {"bias": sum(e) / n, "rms_error": math.sqrt(sum(v * v for v in e) / n)}
            for key, e in zip(("u", "v", "horizontal_speed"), errors)
        },
        "reverse_streamwise_counts": {
            "measured": measured_reverse,
            "computed": predicted_reverse,
            "sign_disagreement": missed_reverse,
        },
    }


def plot_plane(points, destination, height):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize

    x, y, measured, predicted = np.asarray(points).T
    fig, axes = plt.subplots(
        1, 3, figsize=(12, 4.25), sharex=True, sharey=True, layout="constrained"
    )
    for ax, values, label, norm, cmap in zip(
        axes,
        (measured, predicted, predicted - measured),
        ("Wind tunnel", "Kutta · steady SST", "Kutta minus tunnel"),
        (Normalize(-0.5, 1.2), Normalize(-0.5, 1.2), Normalize(-0.45, 0.45)),
        ("viridis", "viridis", "RdBu_r"),
        strict=True,
    ):
        dots = ax.scatter(
            x,
            y,
            c=values,
            cmap=cmap,
            norm=norm,
            s=1,
            marker="s",
            linewidths=0,
            rasterized=True,
        )
        ax.set(aspect="equal", xlabel="Body x (m)", title=label)
        fig.colorbar(
            dots,
            ax=ax,
            orientation="horizontal",
            shrink=0.85,
            pad=0.12,
            label="Δu/U∞" if ax is axes[-1] else "u/U∞",
            extend="both",
        )
    axes[0].set_ylabel("Body y (m)")
    fig.suptitle(
        f"Windsor squareback wake at z = {height:.3f} m", fontsize=16, weight="bold"
    )
    fig.savefig(destination, dpi=160)
    plt.close(fig)
