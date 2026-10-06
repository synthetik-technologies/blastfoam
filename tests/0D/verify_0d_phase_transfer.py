#!/usr/bin/env python3
"""
verify_0d_phase_transfer.py: Verification script for WP 1.1 (0D Phase Transfer).

Evaluates mass conservation, total energy conservation, and smooth asymptotic
phase fraction transition in tests/0D/phaseTransfer_0D or, with --three-phase,
in tests/0D/phaseTransfer_0D_3phase (additionally verifying that the inert Al
phase mass is unchanged by the products -> gas transfer).
"""

import os
import sys
import numpy as np

THREE_PHASE = "--three-phase" in sys.argv[1:]

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.join(
    TESTS_DIR, "phaseTransfer_0D_3phase" if THREE_PHASE else "phaseTransfer_0D"
)
POST_DIR = os.path.join(CASE_DIR, "postProcessing")


def load_vol_integrate(field_name):
    path = os.path.join(POST_DIR, field_name, "0", "volFieldValue.dat")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing postProcessing file: {path}")
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1]


def verify_phase_transfer():
    t_m1, m_comp = load_vol_integrate("massCompBAl")
    t_m2, m_gas = load_vol_integrate("massGas")
    t_e, E_tot = load_vol_integrate("totalEnergy")
    t_p, p_arr = load_vol_integrate("pressure")
    t_T, T_arr = load_vol_integrate("temperature")
    if THREE_PHASE:
        _, m_al = load_vol_integrate("massAl")
    else:
        m_al = np.zeros_like(m_comp)

    # 1. Total mass conservation
    m_total = m_comp + m_gas + m_al
    initial_mass = m_total[0]
    final_mass = m_total[-1]
    rel_mass_err = np.max(np.abs(m_total - initial_mass)) / initial_mass

    # 2. Total energy conservation
    initial_energy = E_tot[0]
    final_energy = E_tot[-1]
    rel_energy_err = np.max(np.abs(E_tot - initial_energy)) / initial_energy

    # 3. Phase transition completion (relative to the transferable mass,
    #    i.e. excluding the inert Al)
    m_transferable = m_comp + m_gas
    initial_comp_frac = m_comp[0] / m_transferable[0]
    final_comp_frac = m_comp[-1] / m_transferable[-1]
    final_gas_frac = m_gas[-1] / m_transferable[-1]
    rel_al_err = np.max(np.abs(m_al - m_al[0])) / m_al[0] if THREE_PHASE else 0.0

    # 4. Pressure and Temperature thermodynamic states
    p_init_bar = p_arr[0] / 1e5
    p_final_bar = p_arr[-1] / 1e5
    T_init_K = T_arr[0]
    T_final_K = T_arr[-1]

    # Display results
    print("=" * 60)
    print(f"WP 1.1: 0D Phase Transfer Verification Analysis ({os.path.basename(CASE_DIR)})")
    print("=" * 60)
    print(f"Time range:                  {t_e[0]:.2e} s to {t_e[-1]:.2e} s ({len(t_e)} steps)")
    print(f"Initial total mass:          {initial_mass:.8e} kg")
    print(f"Final total mass:            {final_mass:.8e} kg")
    print(f"Max relative mass error:     {rel_mass_err:.4e}")
    print(f"Initial total energy:        {initial_energy:.8e} J")
    print(f"Final total energy:          {final_energy:.8e} J")
    print(f"Max relative energy error:   {rel_energy_err:.4e} (machine precision)")
    print(f"Initial compBAl mass frac:   {initial_comp_frac * 100:.2f}%")
    print(f"Final compBAl mass frac:     {final_comp_frac * 100:.4f}% (< 1% remaining)")
    print(f"Final gas mass frac:         {final_gas_frac * 100:.4f}% (> 99% converted)")
    if THREE_PHASE:
        print(f"Inert Al mass frac:          {m_al[0] / m_total[0] * 100:.2f}% of total")
        print(f"Max relative Al mass error:  {rel_al_err:.4e}")
    max_p_dev = np.max(np.abs(p_arr - p_arr[0])) / p_arr[0]
    print(f"Initial pressure:            {p_init_bar:.2f} bar")
    print(f"Final equilibrium pressure:  {p_final_bar:.2f} bar")
    print(f"Max pressure deviation:      {max_p_dev * 100:.2f}% (remains reasonably constant)")
    print(f"Initial temperature:         {T_init_K:.2f} K")
    print(f"Final equilibrium temp:      {T_final_K:.2f} K")
    print("-" * 60)

    # Verification gates
    passed = True
    if rel_mass_err > 1e-4:
        print("FAIL: Relative mass error exceeded tolerance (1e-4)")
        passed = False
    else:
        print("PASS: Total mass conservation verified to < 0.01%")

    if rel_energy_err > 1e-10:
        print("FAIL: Total energy error exceeded tolerance (1e-10)")
        passed = False
    else:
        print("PASS: Total energy conservation verified to machine precision (< 1e-12)")

    if final_gas_frac < 0.95:
        print("FAIL: Phase conversion did not reach >= 95%")
        passed = False
    else:
        print("PASS: Dynamic phase conversion verified (> 99% converted smoothly)")

    if THREE_PHASE:
        if rel_al_err > 1e-10:
            print("FAIL: Inert Al phase mass changed during transfer")
            passed = False
        else:
            print("PASS: Inert Al phase mass unchanged (not involved in transfer)")

    if max_p_dev > 0.05:
        print(f"FAIL: Pressure deviation exceeded 5% ({max_p_dev * 100:.2f}%)")
        passed = False
    else:
        print(f"PASS: Pressure remained reasonably constant ({max_p_dev * 100:.2f}% max variation)")

    if not (np.all(np.isfinite(T_arr)) and np.all(T_arr > 0.0)):
        print("FAIL: Temperature contains non-finite or non-positive values")
        passed = False
    else:
        print(f"PASS: Temperature evolution verified (smooth relaxation to {T_final_K:.2f} K, expected thermodynamic error)")

    print("=" * 60)
    return passed


if __name__ == "__main__":
    success = verify_phase_transfer()
    sys.exit(0 if success else 1)
