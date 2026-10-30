#!/usr/bin/env python3
"""Spectrum of the frozen-velocity pressure iteration, not the full flow Jacobian."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from analyze import read


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--faces", type=Path, help="Optional captured provisional faces for a fixed-state arithmetic comparison")
    args = parser.parse_args()
    root = args.directory
    J = np.loadtxt(root / "pressure_response.csv", delimiter=",")
    L = np.loadtxt(root / "pressure_correction.csv", delimiter=",")
    alpha = np.loadtxt(root / "pressure_alpha.csv")
    scale = np.sqrt(np.diag(L))
    normalized_L = L / scale[:, None] / scale[None, :]
    normalized_J = J / scale[:, None] / scale[None, :]
    response = alpha[:, None] * np.linalg.solve(normalized_L, normalized_J)
    eigenvalues = np.linalg.eigvals(response)
    iteration = 1 - eigenvalues
    result = {"rows": len(L), "symmetric_correction_error": float(np.max(abs(L - L.T))),
              "normalized_correction_condition": float(np.linalg.cond(normalized_L)),
              "frozen_velocity_spectral_radius": float(np.max(abs(iteration))),
              "unstable_modes_above_1_plus_1e_6": int(np.count_nonzero(abs(iteration) > 1 + 1e-6)),
              "minimum_response_eigenvalue_real": float(np.min(eigenvalues.real)),
              "maximum_response_eigenvalue_real": float(np.max(eigenvalues.real)),
              "scope": "Frozen momentum diagonal, velocity and turbulence at captured coarse-plate state; existing pressure reconstruction evaluated in host double with stored float coefficients. This is not the spectrum of the full pressure-velocity iteration."}
    (root / "spectrum.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if args.faces:
        pressure = np.loadtxt(root / "pressure_state.csv")
        physical = np.loadtxt(root / "rhie_chow_physical.csv")
        faces = list(csv.DictReader(read(args.faces).splitlines()))
        if len(faces) != len(physical):
            raise ValueError("Captured face count mismatch")
        mass, stored_jump, evaluated_flux = np.zeros(len(L)), np.zeros(len(L)), np.zeros(len(L))
        for face in faces:
            f, owner, neighbour = int(face["face"]), int(face["owner"]), int(face["neighbour"])
            advected = sum(float(np.float32(float(face[area]))) * float(face[velocity])
                           for area, velocity in (("sx", "uf"), ("sy", "vf"), ("sz", "wf")))
            impermeable = face["patch"] in ("front", "back") or face["patch"].startswith(("symmetry_", "viscous_solid_"))
            jump = physical[f] * float(face["jump"])
            evaluated = 0 if impermeable else advected + jump
            for cell, sign in ((owner, 1), (neighbour, -1)):
                if cell < 0:
                    continue
                mass[cell] += sign * float(face["flux"])
                stored_jump[cell] += sign * jump
                evaluated_flux[cell] += sign * evaluated
        double_pressure = J @ pressure
        delta = double_pressure - stored_jump
        arithmetic = {
            "stored_mass_l1": float(abs(mass).sum()),
            "pressure_arithmetic_change_l1": float(abs(delta).sum()),
            "pressure_only_fixed_state_counterfactual_l1": float(abs(mass + delta).sum()),
            "face_flux_arithmetic_change_l1": float(abs(evaluated_flux - mass).sum()),
            "double_evaluated_face_flux_l1": float(abs(evaluated_flux).sum()),
            "both_fixed_state_counterfactual_l1": float(abs(evaluated_flux + delta).sum()),
            "scope": "Same stored cell pressures, face velocities and physical Rhie-Chow coefficients, density=1 flat plate. Evaluate final face-flux and pressure-jump arithmetic in double without resolving the flow. This does not test cell-field storage precision or the full iteration Jacobian."
        }
        (root / "pressure-arithmetic.json").write_text(json.dumps(arithmetic, indent=2) + "\n")
        print(json.dumps(arithmetic, indent=2))


if __name__ == "__main__":
    main()
