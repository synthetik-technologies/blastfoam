#!/usr/bin/env python3
"""
Comparison script between eJWLC and Noble-Abel (Abel-Noble) Equations of State (EOS)
for the ChargeCompBAl_1D blastFoam case.

Features:
  - eJWLC isentrope plotted from Chapman-Jouguet (CJ) state down to 1 bar.
  - Abel-Noble isentrope anchored at the 100 bar transfer point, extended from 500 bar down to 1 bar.
  - Enforces and reports exact mass and energy balance at the 100 bar transfer point.
  - Fully dynamic parsing from constant/physicalProperties, constant/fvModels, and constant/thermo.compressibleGas.
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumulative_trapezoid

RU = 8314.4626  # Universal gas constant [J / (kmol K)]


def parse_physical_properties(filepath):
    """Extract eJWLC products parameters from physicalProperties dictionary."""
    with open(filepath, 'r') as f:
        content = f.read()

    # Match equationOfState subdict inside products block
    m = re.search(r'products\s*\{.*?equationOfState\s*\{([^}]+)\}', content, re.DOTALL)
    block = m.group(1) if m else content

    def get_val(key, default=None):
        m_val = re.search(rf'\b{key}\s+([+-]?[\d.eE+-]+)\s*;', block)
        return float(m_val.group(1)) if m_val else default

    params = {
        'rho0': get_val('rho0', 1624.45),
        'A': get_val('A', 365.113e9),
        'B': get_val('B', 2.656e9),
        'C': get_val('C', 1.043e9),
        'R1': get_val('R1', 3.83395),
        'R2': get_val('R2', 0.8997),
        'omega': get_val('omega', 0.2387),
        'Cv': get_val('Cv', 153.3),
        'T_CJ': get_val('T_CJ', 3432.62),
        'rho_CJ': get_val('rho_CJ', 2098.7),
        'e0': get_val('e0', 72078.0),
    }
    return params


def parse_fv_models(filepath):
    """Extract phase change targetFractions from fvModels dictionary."""
    with open(filepath, 'r') as f:
        content = f.read()

    tf_match = re.search(r'targetFractions\s*\{([^}]+)\}', content)
    fractions = {}
    if tf_match:
        lines = tf_match.group(1).strip().split('\n')
        for line in lines:
            line = re.sub(r'//.*', '', line).strip()
            if not line:
                continue
            parts = re.split(r'\s+', line.rstrip(';').strip())
            if len(parts) >= 2:
                try:
                    fractions[parts[0]] = float(parts[1])
                except ValueError:
                    pass

    t_thresh_m = re.search(r'Tthreshold\s+([+-]?[\d.eE+-]+)\s*;', content)
    p_thresh_m = re.search(r'pThreshold\s+([+-]?[\d.eE+-]+)\s*;', content)

    t_thresh = float(t_thresh_m.group(1)) if t_thresh_m else 1600.0
    p_thresh = float(p_thresh_m.group(1)) if p_thresh_m else 1.0e7

    return fractions, t_thresh, p_thresh


def parse_thermo_compressible_gas(filepath):
    """Extract species properties (molWeight, b, JANAF coeffs) from thermo.compressibleGas."""
    with open(filepath, 'r') as f:
        content = f.read()

    blocks = {}
    lines = content.split('\n')
    current_key = None
    brace_level = 0
    block_lines = []
    for line in lines:
        stripped = line.strip()
        if not current_key:
            m = re.match(r'^([A-Za-z0-9_()+.-]+)\s*\{?$', stripped)
            if m and m.group(1) not in ('species', 'FoamFile', '('):
                current_key = m.group(1)
                brace_level = line.count('{') - line.count('}')
                block_lines = [line]
                if brace_level == 0 and '{' in line and '}' in line:
                    blocks[current_key] = '\n'.join(block_lines)
                    current_key = None
        else:
            block_lines.append(line)
            brace_level += line.count('{') - line.count('}')
            if brace_level == 0 and '{' in ''.join(block_lines):
                blocks[current_key] = '\n'.join(block_lines)
                current_key = None

    species_db = {}
    for sp, sp_block in blocks.items():
        mw_m = re.search(r'molWeight\s+([+-]?[\d.eE+-]+)\s*;', sp_block)
        b_m = re.search(r'\bb\s+([+-]?[\d.eE+-]+)\s*;', sp_block)

        high_m = re.search(r'highCpCoeffs\s*\(\s*([^)]+)\s*\);', sp_block)
        low_m = re.search(r'lowCpCoeffs\s*\(\s*([^)]+)\s*\);', sp_block)
        tlow_m = re.search(r'Tlow\s+([+-]?[\d.eE+-]+)\s*;', sp_block)
        thigh_m = re.search(r'Thigh\s+([+-]?[\d.eE+-]+)\s*;', sp_block)
        tcom_m = re.search(r'Tcommon\s+([+-]?[\d.eE+-]+)\s*;', sp_block)

        mw = float(mw_m.group(1)) if mw_m else 28.0
        b_val = float(b_m.group(1)) if b_m else 1.0e-3
        high_coeffs = [float(x) for x in high_m.group(1).split()] if high_m else [3.5, 0, 0, 0, 0, 0, 0]
        low_coeffs = [float(x) for x in low_m.group(1).split()] if low_m else [3.5, 0, 0, 0, 0, 0, 0]

        species_db[sp] = {
            'molWeight': mw,
            'b': b_val,
            'highCpCoeffs': high_coeffs,
            'lowCpCoeffs': low_coeffs,
            'Tlow': float(tlow_m.group(1)) if tlow_m else 200.0,
            'Thigh': float(thigh_m.group(1)) if thigh_m else 5000.0,
            'Tcommon': float(tcom_m.group(1)) if tcom_m else 1000.0,
        }

    return species_db


def calc_cp_specie(T, sp_data):
    """Calculate cp [J/(kg K)] for a single species using NASA 7-coefficient polynomials."""
    R_sp = RU / sp_data['molWeight']
    T = np.clip(T, sp_data['Tlow'], sp_data['Thigh'])
    a = sp_data['highCpCoeffs'] if T >= sp_data['Tcommon'] else sp_data['lowCpCoeffs']
    cp_over_R = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4
    return cp_over_R * R_sp


def calc_h_specie(T, sp_data):
    """Calculate sensible+formation enthalpy h [J/kg] using NASA 7-coefficient polynomials."""
    R_sp = RU / sp_data['molWeight']
    T = np.clip(T, sp_data['Tlow'], sp_data['Thigh'])
    a = sp_data['highCpCoeffs'] if T >= sp_data['Tcommon'] else sp_data['lowCpCoeffs']
    h_over_R = a[0]*T + a[1]/2.0*T**2 + a[2]/3.0*T**3 + a[3]/4.0*T**4 + a[4]/5.0*T**5 + a[5]
    return h_over_R * R_sp


def compute_mixture_properties(fractions, species_db):
    """Compute mixture-averaged molecular weight, gas constant, co-volume, and thermal functions."""
    tot_frac = sum(fractions.values())
    norm_fracs = {k: v / tot_frac for k, v in fractions.items()}

    inv_W = sum(norm_fracs[k] / species_db[k]['molWeight'] for k in norm_fracs if k in species_db)
    W_mix = 1.0 / inv_W
    R_mix = RU / W_mix
    b_mix = sum(norm_fracs[k] * species_db[k]['b'] for k in norm_fracs if k in species_db)

    def cp_mix(T):
        return sum(norm_fracs[k] * calc_cp_specie(T, species_db[k]) for k in norm_fracs if k in species_db)

    def cv_mix(T):
        return cp_mix(T) - R_mix

    def gamma_mix(T):
        return cp_mix(T) / cv_mix(T)

    def h_mix(T):
        return sum(norm_fracs[k] * calc_h_specie(T, species_db[k]) for k in norm_fracs if k in species_db)

    def e_mix(T):
        return h_mix(T) - R_mix * T

    return {
        'W_mix': W_mix,
        'R_mix': R_mix,
        'b_mix': b_mix,
        'norm_fracs': norm_fracs,
        'cp_mix': cp_mix,
        'cv_mix': cv_mix,
        'gamma_mix': gamma_mix,
        'h_mix': h_mix,
        'e_mix': e_mix,
    }


def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    phys_prop_file = os.path.join(case_dir, 'constant', 'physicalProperties')
    fv_models_file = os.path.join(case_dir, 'constant', 'fvModels')
    thermo_gas_file = os.path.join(case_dir, 'constant', 'thermo.compressibleGas')

    # 1. Parse parameters from OpenFOAM case files
    jwl = parse_physical_properties(phys_prop_file)
    fractions, t_thresh, p_thresh = parse_fv_models(fv_models_file)
    species_db = parse_thermo_compressible_gas(thermo_gas_file)
    mix = compute_mixture_properties(fractions, species_db)

    print("=" * 78)
    print("  eJWLC (CJ to 1 bar) vs. Abel-Noble (500 bar to 1 bar) with 100 bar Transfer")
    print("=" * 78)
    print(f"eJWLC Detonation Products Parameters (constant/physicalProperties):")
    print(f"  rho0    = {jwl['rho0']:.2f} kg/m^3  |  rho_CJ = {jwl['rho_CJ']:.2f} kg/m^3  |  T_CJ = {jwl['T_CJ']:.2f} K")
    print(f"  A       = {jwl['A']/1e9:.3f} GPa    |  B      = {jwl['B']/1e9:.3f} GPa      |  C    = {jwl['C']/1e9:.4f} GPa")
    print(f"  R1      = {jwl['R1']:.5f}      |  R2     = {jwl['R2']:.4f}       |  omega = {jwl['omega']:.4f}")
    print(f"  Cv      = {jwl['Cv']:.2f} J/(kg K) |  e0     = {jwl['e0']:.1f} J/kg")
    print("-" * 78)
    print(f"Abel-Noble Reacting Gas Mixture (constant/fvModels + thermo.compressibleGas):")
    print(f"  Species:                  {', '.join(mix['norm_fracs'].keys())}")
    print(f"  Mixture Molecular Weight: {mix['W_mix']:.3f} g/mol")
    print(f"  Mixture Gas Constant R:   {mix['R_mix']:.2f} J/(kg K)")
    print(f"  Mixture Co-volume b_mix:  {mix['b_mix']:.5e} m^3/kg ({mix['b_mix']*1e3:.4f} cm^3/g)")
    print("=" * 78)

    # 2. Master expansion grid from CJ state down to 1 bar
    v0 = 1.0 / jwl['rho0']
    v_CJ = 1.0 / jwl['rho_CJ']
    V_CJ = v_CJ / v0

    V_max_search = 3000.0
    V_grid = np.geomspace(V_CJ, V_max_search, 6000)
    v_grid = V_grid * v0

    A, B, C = jwl['A'], jwl['B'], jwl['C']
    R1, R2, om = jwl['R1'], jwl['R2'], jwl['omega']
    e0 = jwl['e0']

    # eJWLC state along isentrope
    p_jwl_grid = A * np.exp(-R1 * V_grid) + B * np.exp(-R2 * V_grid) + C * (V_grid ** (-(om + 1.0)))
    T_jwl_grid = jwl['T_CJ'] * (V_CJ / V_grid) ** om
    e_jwl_grid = v0 * ((A/R1)*np.exp(-R1*V_grid) + (B/R2)*np.exp(-R2*V_grid) + (C/om)*(V_grid**(-om))) + e0

    # Locate key pressure points along the expansion
    P_500 = 500.0 * 1e5   # 500 bar
    P_TRANS = 100.0 * 1e5 # 100 bar (Transfer point)
    P_1BAR = 1.0 * 1e5    # 1 bar

    idx_500bar = np.argmin(np.abs(p_jwl_grid - P_500))
    idx_trans = np.argmin(np.abs(p_jwl_grid - P_TRANS))
    idx_1bar = np.argmin(np.abs(p_jwl_grid - P_1BAR))

    # Truncate eJWLC array from CJ state down to 1 bar
    V_jwl = V_grid[:idx_1bar + 1]
    v_jwl = v_grid[:idx_1bar + 1]
    p_jwl = p_jwl_grid[:idx_1bar + 1]
    T_jwl = T_jwl_grid[:idx_1bar + 1]
    e_jwl = e_jwl_grid[:idx_1bar + 1]

    # Isentropic expansion velocity for eJWLC: u(v) = int sqrt(-dp/dv) dv
    dp_dV_jwl = -A * R1 * np.exp(-R1 * V_jwl) - B * R2 * np.exp(-R2 * V_jwl) - C * (om + 1.0) * (V_jwl ** (-(om + 2.0)))
    dp_dv_jwl = dp_dV_jwl / v0
    integrand_u_jwl = np.sqrt(np.maximum(-dp_dv_jwl, 0.0))
    u_jwl = cumulative_trapezoid(integrand_u_jwl, v_jwl, initial=0.0)

    # State parameters at the 100 bar phase transfer point from eJWLC
    v_trans = v_grid[idx_trans]
    p_trans_jwl = p_jwl_grid[idx_trans]
    T_trans_jwl = T_jwl_grid[idx_trans]
    e_trans_jwl = e_jwl_grid[idx_trans]
    u_trans = u_jwl[idx_trans]

    # Energy-matching temperature determination for target reactive gas:
    # In compressible solvers, mass (v) and internal energy (e) are strictly conserved.
    # We solve for T_trans such that e_sensible,gas(T_trans) = e_trans,JWL:
    from scipy.optimize import root_scalar
    e_ref_298 = mix['e_mix'](298.15)
    def energy_residual(T):
        return (mix['e_mix'](T) - e_ref_298) - e_trans_jwl

    res_T = root_scalar(energy_residual, bracket=[300.0, 4500.0])
    T_trans_na = res_T.root
    p_trans_na = (mix['R_mix'] * T_trans_na) / np.maximum(v_trans - mix['b_mix'], 1e-8)

    # Gruneisen parameter comparison:
    cv_trans = mix['cv_mix'](T_trans_na)
    gamma_trans = mix['gamma_mix'](T_trans_na)
    Gruneisen_na = (gamma_trans - 1.0) / (1.0 - mix['b_mix']/v_trans)
    Gruneisen_jwl = om

    print("\n--- Phase Transfer State & Balances @ 100 bar ---")
    print(f"  [MASS BALANCE]              rho_gas = rho_compBAl = {1.0/v_trans:.3f} kg/m^3  (v_trans = {v_trans:.5e} m^3/kg)")
    print(f"  [ENERGY BALANCE]            e_sensible,gas = e_JWL(v_trans) = {e_trans_jwl/1e6:.4f} MJ/kg")
    print(f"                              T_gas (from energy match) = {T_trans_na:.1f} K  (vs. T_JWL = {T_trans_jwl:.1f} K)")
    print(f"  [PRESSURE MISMATCH]         p_NA = {p_trans_na/1e5:.2f} bar  vs.  p_JWL = {p_trans_jwl/1e5:.2f} bar")
    print(f"  [GRUNEISEN EXPLANATION]     Gamma_NA = {Gruneisen_na:.4f}  vs.  omega_JWL = {Gruneisen_jwl:.4f}  (Ratio = {Gruneisen_na/Gruneisen_jwl:.2f}x)")
    print(f"                              Since p_thermal ~ Gamma * rho * e, Noble-Abel naturally yields {Gruneisen_na/Gruneisen_jwl:.2f}x higher pressure")
    print(f"                              for the EXACT SAME internal energy and mass density.")
    print("-" * 78)

    # 3. Integrate Single Continuous Abel-Noble Isentrope anchored at 100 bar transfer point:
    #    - Backward from v_trans down to v_500bar (compressing to 500 bar)
    #    - Forward from v_trans up to v_1bar (expanding to 1 bar)
    b_mix = mix['b_mix']
    R_mix = mix['R_mix']

    v_500 = v_grid[idx_500bar]
    v_1bar = v_grid[idx_1bar]

    def ode_isentrope(v, T_state):
        T = max(float(T_state[0]), 50.0)
        cv = mix['cv_mix'](T)
        denom = max(v - b_mix, 1e-8)
        dT_dv = - (R_mix * T) / (cv * denom)
        return [dT_dv]

    # Grid for Abel-Noble: spans [v_500, v_1bar]
    v_na = v_grid[idx_500bar:idx_1bar + 1]

    # Partition into backward (500 bar to 100 bar) and forward (100 bar to 1 bar)
    v_na_back = v_grid[idx_500bar:idx_trans + 1][::-1]  # from v_trans down to v_500
    v_na_fwd = v_grid[idx_trans:idx_1bar + 1]            # from v_trans up to v_1bar

    # Solve backward from v_trans to v_500 using T_trans_na
    sol_back = solve_ivp(ode_isentrope, [v_trans, v_500], [T_trans_na], t_eval=v_na_back, rtol=1e-9, atol=1e-9)
    T_na_back = sol_back.y[0][::-1]

    # Solve forward from v_trans to v_1bar using T_trans_na
    sol_fwd = solve_ivp(ode_isentrope, [v_trans, v_1bar], [T_trans_na], t_eval=v_na_fwd, rtol=1e-9, atol=1e-9)
    T_na_fwd = sol_fwd.y[0]

    # Concatenate to form a seamless single curve from 500 bar to 1 bar
    T_na = np.concatenate([T_na_back[:-1], T_na_fwd])
    p_na = (R_mix * T_na) / np.maximum(v_na - b_mix, 1e-8)

    # Velocity along Abel-Noble isentrope anchored at u_trans
    gamma_na = np.array([mix['gamma_mix'](t) for t in T_na])
    dp_dv_na = - (p_na * gamma_na) / np.maximum(v_na - b_mix, 1e-8)
    integrand_u_na = np.sqrt(np.maximum(-dp_dv_na, 0.0))

    idx_trans_in_na = idx_trans - idx_500bar
    int_u_na_from_start = cumulative_trapezoid(integrand_u_na, v_na, initial=0.0)
    u_na = u_trans + (int_u_na_from_start - int_u_na_from_start[idx_trans_in_na])

    # Print summary states
    print(f"Abel-Noble Isentrope Boundaries:")
    print(f"  At 500 bar region: v = {v_na[0]:.5e} m^3/kg, p = {p_na[0]/1e5:.1f} bar, T = {T_na[0]:.1f} K, u = {u_na[0]:.1f} m/s")
    print(f"  At 100 bar transfer: v = {v_trans:.5e} m^3/kg, p = {p_na[idx_trans_in_na]/1e5:.1f} bar, T = {T_na[idx_trans_in_na]:.1f} K, u = {u_na[idx_trans_in_na]:.1f} m/s")
    print(f"  At 1 bar endpoint:  v = {v_na[-1]:.5e} m^3/kg, p = {p_na[-1]/1e5:.2f} bar, T = {T_na[-1]:.1f} K, u = {u_na[-1]:.1f} m/s")
    print("-" * 78)

    # Load EuDetCode products isentrope if available
    eudet_prod_file = os.path.join(case_dir, 'products_isentrope_output.dat')
    if not os.path.exists(eudet_prod_file):
        eudet_prod_file = os.path.join(case_dir, 'isentrope_output.dat')
    
    eudet_data = None
    if os.path.exists(eudet_prod_file):
        raw_dat = np.loadtxt(eudet_prod_file)
        if raw_dat.shape[1] >= 6:
            v_eu = raw_dat[:, 0]
            p_eu = raw_dat[:, 1]
            T_eu = raw_dat[:, 2]
            dp_dv_eu = np.gradient(p_eu, v_eu)
            integrand_u_eu = np.sqrt(np.maximum(-dp_dv_eu, 0.0))
            u_eu = cumulative_trapezoid(integrand_u_eu, v_eu, initial=0.0)
            idx_eu_100 = np.argmin(np.abs(p_eu - P_TRANS))
            eudet_data = {
                'v': v_eu,
                'p': p_eu,
                'T': T_eu,
                'u': u_eu,
                'idx_100': idx_eu_100,
                'label': r'EuDetCode Isentrope ($v_{\mathrm{prod}}$)'
            }
            print(f"Loaded EuDetCode isentrope from: {eudet_prod_file}")
            print(f"  EuDetCode @ 100 bar: v = {v_eu[idx_eu_100]*1e3:.3f} cm^3/g, T = {T_eu[idx_eu_100]:.1f} K, u = {u_eu[idx_eu_100]:.1f} m/s")

    # 4. Generate High-Quality Comparison Plots
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- Subplot 1: p-v Expansion Curve (Log-Log) ---
    ax = axes[0, 0]
    if eudet_data is not None:
        ax.loglog(eudet_data['v'] * 1e3, eudet_data['p'] / 1e5, color='forestgreen', linestyle='-.', lw=2.2, label=eudet_data['label'])
    ax.loglog(v_jwl * 1e3, p_jwl / 1e5, 'b-', lw=2.2, label='eJWLC Products Fit (CJ to 1 bar)')
    ax.loglog(v_na * 1e3, p_na / 1e5, 'r--', lw=2.2, label='Abel-Noble Isentrope (500 bar to 1 bar)')
    ax.scatter([v_trans * 1e3], [p_trans_jwl / 1e5], color='blue', s=70, zorder=5, label=f'eJWLC @ 100 bar')
    ax.scatter([v_trans * 1e3], [p_trans_na / 1e5], color='red', s=70, zorder=5, label=f'Abel-Noble @ 100 bar')
    if eudet_data is not None:
        idx_eu = eudet_data['idx_100']
        ax.scatter([eudet_data['v'][idx_eu] * 1e3], [eudet_data['p'][idx_eu] / 1e5], color='forestgreen', s=70, zorder=5, label='EuDetCode @ 100 bar')
    ax.axhline(500, color='gray', linestyle=':', alpha=0.6, label='500 bar Bound')
    ax.axhline(100, color='darkorange', linestyle=':', alpha=0.6)
    ax.axhline(1.0, color='gray', linestyle=':', alpha=0.6)
    ax.set_ylim(0.8, p_jwl[0] / 1e5 * 1.5)
    ax.set_xlabel(r'Specific Volume $v$ [$\mathrm{cm}^3/\mathrm{g}$ or $10^{-3}\,\mathrm{m}^3/\mathrm{kg}$]', fontsize=11)
    ax.set_ylabel(r'Pressure $p$ [bar]', fontsize=11)
    ax.set_title(r'Pressure vs. Specific Volume ($p$--$v$ Isentrope)', fontsize=12, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=8.5, loc='upper right')

    # --- Subplot 2: p-u Curve (Pressure vs Expansion Velocity) ---
    ax = axes[0, 1]
    if eudet_data is not None:
        ax.plot(eudet_data['u'] / 1e3, eudet_data['p'] / 1e5, color='forestgreen', linestyle='-.', lw=2.2, label=r'EuDetCode Expansion ($u_{\mathrm{prod}}$)')
    ax.plot(u_jwl / 1e3, p_jwl / 1e5, 'b-', lw=2.2, label='eJWLC Expansion (CJ to 1 bar)')
    ax.plot(u_na / 1e3, p_na / 1e5, 'r--', lw=2.2, label='Abel-Noble Expansion (500 bar to 1 bar)')
    ax.scatter([u_trans / 1e3], [p_trans_jwl / 1e5], color='blue', s=70, zorder=5, label='eJWLC @ 100 bar')
    ax.scatter([u_trans / 1e3], [p_trans_na / 1e5], color='red', s=70, zorder=5, label='Abel-Noble @ Transfer')
    ax.set_yscale('log')
    ax.set_ylim(0.8, p_jwl[0] / 1e5 * 1.5)
    ax.set_xlabel(r'Expansion Particle Velocity $u$ [$\mathrm{km/s}$]', fontsize=11)
    ax.set_ylabel(r'Pressure $p$ [bar]', fontsize=11)
    ax.set_title(r'Pressure vs. Expansion Velocity ($p$--$u$ Curve)', fontsize=12, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=8.5, loc='upper right')

    # --- Subplot 3: Temperature vs Specific Volume ---
    ax = axes[1, 0]
    if eudet_data is not None:
        ax.plot(eudet_data['v'] * 1e3, eudet_data['T'], color='forestgreen', linestyle='-.', lw=2.2, label=r'EuDetCode Equil. ($T_{\mathrm{true}}$)')
    ax.plot(v_jwl * 1e3, T_jwl, 'b-', lw=2.2, label=r'eJWLC ($T_s = T_{\mathrm{CJ}}(V_{\mathrm{CJ}}/V)^\omega$)')
    ax.plot(v_na * 1e3, T_na, 'r--', lw=2.2, label='Abel-Noble Isentrope (500 bar to 1 bar)')
    ax.scatter([v_trans * 1e3], [T_trans_jwl], color='blue', s=70, zorder=5, label=f'T_JWL ({T_trans_jwl:.0f} K)')
    ax.scatter([v_trans * 1e3], [T_trans_na], color='red', s=70, zorder=5, label=f'T_NA,matched ({T_trans_na:.0f} K)')
    if eudet_data is not None:
        idx_eu = eudet_data['idx_100']
        ax.scatter([eudet_data['v'][idx_eu] * 1e3], [eudet_data['T'][idx_eu]], color='forestgreen', s=70, zorder=5, label=f'T_EuDet ({eudet_data["T"][idx_eu]:.0f} K)')
    ax.set_xscale('log')
    ax.set_xlabel(r'Specific Volume $v$ [$\mathrm{cm}^3/\mathrm{g}$]', fontsize=11)
    ax.set_ylabel(r'Temperature $T$ [K]', fontsize=11)
    ax.set_title(r'Temperature along Isentrope ($T$--$v$ Curve)', fontsize=12, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=8.5, loc='upper right')

    # --- Subplot 4: Relative Deviation & Energy Balance Details ---
    ax = axes[1, 1]
    p_jwl_overlap = p_jwl[idx_500bar:]
    rel_diff = (p_na - p_jwl_overlap) / p_jwl_overlap * 100.0
    ax.plot(v_na * 1e3, rel_diff, 'm-', lw=2.2, label=r'$(p_{\mathrm{NA}} - p_{\mathrm{JWL}}) / p_{\mathrm{JWL}} \times 100\%$')
    if eudet_data is not None:
        # Interpolate eJWLC onto EuDetCode volume points
        p_jwl_on_eu = np.interp(eudet_data['v'], v_jwl, p_jwl)
        rel_diff_jwl_eu = (p_jwl_on_eu - eudet_data['p']) / eudet_data['p'] * 100.0
        ax.plot(eudet_data['v'] * 1e3, rel_diff_jwl_eu, color='forestgreen', linestyle='-.', lw=2.0, label=r'$(p_{\mathrm{JWL}} - p_{\mathrm{EuDet}}) / p_{\mathrm{EuDet}} \times 100\%$')
    ax.axvline(v_trans * 1e3, color='darkorange', linestyle='--', label=f'Transfer $v = {v_trans*1e3:.2f}' + r'\,\mathrm{cm}^3/\mathrm{g}$')
    ax.axhline(0, color='black', lw=1, linestyle='-')
    ax.set_xscale('log')
    ax.set_xlabel(r'Specific Volume $v$ [$\mathrm{cm}^3/\mathrm{g}$]', fontsize=11)
    ax.set_ylabel(r'Relative Pressure Deviation [$\%$]', fontsize=11)
    ax.set_title(r'EOS Deviation & Balance @ 100 bar Transfer', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=8.5, loc='lower right')

    # Add text annotation summarizing mass and energy balance
    balance_lines = [
        "100 bar Conservation Balances:",
        f"  • Mass: v_trans = {v_trans*1e3:.2f} cm³/g (conserved)",
        f"  • Energy: e_sens = {e_trans_jwl/1e6:.3f} MJ/kg (conserved)",
        f"  • T_gas (energy-matched) = {T_trans_na:.1f} K",
        f"  • p_JWL = {p_trans_jwl/1e5:.1f} bar vs p_NA = {p_trans_na/1e5:.1f} bar",
        f"  • Ratio Gamma_NA/omega_JWL = {Gruneisen_na/Gruneisen_jwl:.2f}x",
    ]
    if eudet_data is not None:
        idx_eu = eudet_data['idx_100']
        balance_lines.append(f"  • EuDetCode T_true @ 100 bar = {eudet_data['T'][idx_eu]:.1f} K")
    balance_text = "\n".join(balance_lines)

    ax.text(
        0.05, 0.95, balance_text,
        transform=ax.transAxes,
        fontsize=8.5,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='linen', edgecolor='darkorange', alpha=0.9)
    )

    plt.tight_layout()
    plot_path = os.path.join(case_dir, 'pv_curve_comparison.png')
    plt.savefig(plot_path, dpi=300)
    print(f"\nSaved updated comparison plot to:\n  {plot_path}")
    print("=" * 78)


if __name__ == '__main__':
    main()
