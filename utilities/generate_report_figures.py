#!/usr/bin/env python3
"""
generate_report_figures.py: figures of report/TBX-blastFoam.tex.

Every figure is made from data: the EuDetCode output stored in
utilities/thermo_fitter/eudet_compBAl, the reference models of the test
scripts, and the results of the test cases. The test results are read from
$BLAST_RESULTS_DIR/tests/... (default: the repository), so the figures can
be made from a complete run of the tests in a scratch copy:

    BLAST_RESULTS_DIR=/path/to/copy python3 utilities/generate_report_figures.py [name ...]

A missing input is an error; no figure is drawn from placeholder data.
"""

import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = "/opt/OpenFOAM/blastfoam"
RESULTS = os.environ.get("BLAST_RESULTS_DIR", CODE)
FIG_DIR = os.path.join(CODE, "report/figures")
EUDET = os.path.join(CODE, "utilities/thermo_fitter/eudet_compBAl")
GAS_CASE = os.path.join(CODE, "tests/0D/phaseTransfer_0D")

for sub in ("utilities/thermo_fitter", "tests/0D", "tests/1D", "tests/1D_spherical",
            "tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m"):
    sys.path.insert(0, os.path.join(CODE, sub))

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "legend.fontsize": 8.5,
    "lines.linewidth": 1.6,
    "grid.alpha": 0.35,
    "grid.linestyle": "--",
})


def result(*parts):
    """Path of a test result; error if it does not exist."""
    path = os.path.join(RESULTS, *parts)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found: run the test first (BLAST_RESULTS_DIR={RESULTS})"
        )
    return path


def save(fig, name):
    out = os.path.join(FIG_DIR, name + ".pdf")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    fig.savefig(out.replace(".pdf", ".png"), dpi=200)
    plt.close(fig)
    print(f"Generated {out}")


def read_field(path):
    """Internal field of an ascii OpenFOAM field (x component of vectors)."""
    text = open(path).read()
    body = text.split("internalField", 1)[1].split(";", 1)[0].strip()
    if body.startswith("uniform"):
        v = body.split(None, 1)[1].strip("() ")
        return np.array([float(v.split()[0])])
    lst = body[body.index("(") + 1:body.rindex(")")]
    if "(" in lst:
        return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", lst)])
    return np.array([float(x) for x in lst.split()])


# ---------------------------------------------------------------------------
# Detonation products (EuDetCode output)
# ---------------------------------------------------------------------------

def products_fit():
    import products_only_fit as pof
    return pof, pof.compute(os.path.join(EUDET, "products"), GAS_CASE)


def load_dat(path):
    return np.loadtxt(path, comments="#")


def fig_isentropes():
    """CJ isentropes of the mixture, the products and the inert Al, with fits."""
    mix = load_dat(os.path.join(EUDET, "products/isentrope_output.dat"))
    prod = load_dat(os.path.join(EUDET, "products/products_isentrope_output.dat"))
    fit_prod = load_dat(os.path.join(EUDET, "products/jwl_isentrope_output.dat"))
    fit_mix = load_dat(os.path.join(EUDET, "full/jwl_isentrope_output.dat"))
    _, r = products_fit()
    y_act = r["cj"]["Mass Fraction in Mixture (Y_eq)"]
    v_inert = (mix[:, 0] - y_act*prod[:, 0])/(1.0 - y_act)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.loglog(mix[:, 0]*1e3, mix[:, 1]/1e9, "k-", label=r"Mixture $v_\mathrm{total}$ (EuDetCode)")
    ax1.loglog(fit_mix[:, 0]*1e3, fit_mix[:, 1]/1e9, "k:", label="Full-mixture JWL fit")
    ax1.loglog(prod[:, 0]*1e3, prod[:, 1]/1e9, "b-", label=r"Products $v_\mathrm{prod}$ (EuDetCode)")
    ax1.loglog(fit_prod[:, 0]*1e3, fit_prod[:, 1]/1e9, "b:", label="Products JWL fit (used)")
    ax1.loglog(v_inert*1e3, mix[:, 1]/1e9, "g-", label=r"Inert Al $v_\mathrm{inert}$")
    ax1.set_xlabel(r"Specific volume $v$ [dm$^3$/kg]")
    ax1.set_ylabel(r"Pressure $p$ [GPa]")
    ax1.set_title("(a) CJ isentropes, 70% Comp-B / 30% Al (90% of Al inert)")
    ax1.grid(True, which="both")
    ax1.legend()

    ax2.semilogx(prod[:, 1]/1e5, 100*(fit_prod[:, 1]/prod[:, 1] - 1), "b-", label="Products fit")
    ax2.semilogx(mix[:, 1]/1e5, 100*(fit_mix[:, 1]/mix[:, 1] - 1), "k:", label="Full-mixture fit")
    ax2.axhline(0, color="gray", lw=0.8)
    ax2.set_xlabel("Pressure on the isentrope [bar]")
    ax2.set_ylabel("JWL fit / EuDetCode $-$ 1 [%]")
    ax2.set_title("(b) Fit error along the isentrope")
    ax2.invert_xaxis()
    ax2.grid(True, which="both")
    ax2.legend()
    save(fig, "isentrope_comparison")


def fig_products_temperature():
    """Isentrope temperature: equilibrium vs eJWLC power law."""
    _, r = products_fit()
    prod = load_dat(os.path.join(EUDET, "products/products_isentrope_output.dat"))
    v = prod[:, 0]
    V = v*r["rho0"]
    p_ejwlc = r["A"]*np.exp(-r["R1"]*V) + r["B"]*np.exp(-r["R2"]*V) + r["C"]*V**(-(r["w"] + 1))
    T_w = r["T_CJ"]*(r["V_CJ"]/V)**r["w"]
    T_wT = r["T_CJ"]*(r["V_CJ"]/V)**r["omegaT"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, x, xl in ((ax1, prod[:, 1]/1e5, "Pressure [bar]"),
                      (ax2, v*1e3, r"Specific volume [dm$^3$/kg]")):
        xe = p_ejwlc/1e5 if ax is ax1 else v*1e3
        ax.semilogx(x, prod[:, 2], "k-", label="Equilibrium (EuDetCode)")
        ax.semilogx(xe, T_w, "r--", label=rf"eJWLC, $\omega_T = \omega = {r['w']:.4f}$")
        ax.semilogx(xe, T_wT, "b-.", label=rf"eJWLC, $\omega_T = {r['omegaT']:.4f}$ (used)")
        xt = r["ptrans"]/1e5 if ax is ax1 else r["v_t"]*1e3
        ax.plot([xt], [r["T_t"]], "ko", mfc="none",
                label=rf"Gas EOS at transfer, {r['T_t']:.0f} K")
        ax.set_xlabel(xl)
        ax.set_ylabel("Temperature [K]")
        ax.grid(True, which="both")
    ax1.invert_xaxis()
    ax1.set_title("(a) Temperature against pressure")
    ax2.set_title("(b) Temperature against specific volume")
    ax1.legend()
    save(fig, "temperature_comparison_products")


def gas_isentrope(pof, r, v_end, n=400):
    """Abel-Nobel isentrope of the frozen product gas through the transfer state."""
    gas, Y = r["gas"], r["Y"]
    R, b = r["R_mix"], r["b_mix"]

    def cv(T):
        # c_p - R of the mixture (c_v of the Abel-Nobel gas), from the
        # species enthalpies (JANAF or nasa9)
        c = 0.0
        dT = 1e-3*T
        for sp, y in Y.items():
            ha, Rs, _, _ = pof.species_data(gas, sp, T)
            c += y*((ha(T + dT) - ha(T - dT))/(2*dT) - Rs)
        return c

    vs = np.geomspace(r["v_t"], v_end, n)
    T = r["T_t"]
    Ts = [T]
    for v0, v1 in zip(vs[:-1], vs[1:]):
        # Midpoint step of dT/dv = -p/cv on the isentrope (de = -p dv)
        f0 = -R*T/(v0 - b)/cv(T)
        Tm = T + 0.5*(v1 - v0)*f0
        vm = 0.5*(v0 + v1)
        T = T + (v1 - v0)*(-R*Tm/(vm - b)/cv(Tm))
        Ts.append(T)
    Ts = np.array(Ts)
    return vs, R*Ts/(vs - b), Ts


def fig_transfer_matching():
    """eJWLC, equilibrium and Abel-Nobel gas around the transfer state."""
    pof, r = products_fit()
    prod = load_dat(os.path.join(EUDET, "products/products_isentrope_output.dat"))
    v = prod[:, 0]
    V = v*r["rho0"]
    p_ejwlc = r["A"]*np.exp(-r["R1"]*V) + r["B"]*np.exp(-r["R2"]*V) + r["C"]*V**(-(r["w"] + 1))
    T_wT = r["T_CJ"]*(r["V_CJ"]/V)**r["omegaT"]
    vg_hi, pg_hi, Tg_hi = gas_isentrope(pof, r, v[-1])
    vg_lo, pg_lo, Tg_lo = gas_isentrope(pof, r, r["v_t"]/4)
    vg = np.concatenate((vg_lo[::-1], vg_hi[1:]))
    pg = np.concatenate((pg_lo[::-1], pg_hi[1:]))
    Tg = np.concatenate((Tg_lo[::-1], Tg_hi[1:]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.loglog(v*1e3, prod[:, 1]/1e5, "k-", label="Equilibrium products (EuDetCode)")
    ax1.loglog(v*1e3, p_ejwlc/1e5, "b--", label="eJWLC isentrope")
    ax1.loglog(vg*1e3, pg/1e5, "r-.", label="Abel–Nobel frozen gas, through transfer state")
    ax1.plot([r["v_t"]*1e3], [r["ptrans"]/1e5], "ko", mfc="none", label="Transfer state")
    ax1.set_xlim(5, v[-1]*1e3*1.1)
    ax1.set_ylim(15, 3000)
    ax1.set_xlabel(r"Specific volume [dm$^3$/kg]")
    ax1.set_ylabel("Pressure [bar]")
    ax1.set_title("(a) Isentropes near the transfer pressure")
    ax1.grid(True, which="both")
    ax1.legend()

    ax2.semilogx(prod[:, 1]/1e5, prod[:, 2], "k-", label="Equilibrium (EuDetCode)")
    ax2.semilogx(p_ejwlc/1e5, T_wT, "b--", label=rf"eJWLC, $\omega_T = {r['omegaT']:.4f}$")
    ax2.semilogx(pg/1e5, Tg, "r-.", label="Abel–Nobel frozen gas")
    ax2.plot([r["ptrans"]/1e5], [r["T_t"]], "ko", mfc="none",
             label=rf"$T_\mathrm{{tr}} = {r['T_t']:.1f}$ K")
    ax2.set_xlim(3000, 15)
    ax2.set_ylim(400, 2000)
    ax2.set_xlabel("Pressure [bar]")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_title("(b) Temperature along the expansion")
    ax2.grid(True, which="both")
    ax2.legend()

    save(fig, "pv_curve_comparison")


# ---------------------------------------------------------------------------
# Stage 1: 0D
# ---------------------------------------------------------------------------

def fig_stage1_phase_transfer():
    """0D products -> gas transfer in blastEulerFoam."""
    post = result("tests/0D/phaseTransfer_0D/postProcessing")

    def load(name):
        return np.loadtxt(os.path.join(post, name, "0/volFieldValue.dat"),
                          comments="#", ndmin=2)

    d1, d2 = load("massCompBAl"), load("massGas")
    de, dp, dT = load("totalEnergy"), load("pressure"), load("temperature")
    t = d1[:, 0]*1e3
    m_tot = d1[:, 1] + d2[:, 1]
    y_comp, y_gas = d1[:, 1]/m_tot, d2[:, 1]/m_tot
    E = de[:, 1:].sum(axis=1)
    T_comp = np.where(y_comp > 1e-6, dT[:, -2], np.nan)
    T_gas = np.where(y_gas > 1e-6, dT[:, -1], np.nan)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.5, 7.6))
    ax1.plot(t, y_comp, "b-", label="Products (eJWLC)")
    ax1.plot(t, y_gas, "r--", label="Gas (JANAF / Abel–Nobel)")
    ax1.set_xlabel("Time [ms]")
    ax1.set_ylabel("Mass fraction [-]")
    ax1.set_title(rf"(a) Phase mass fractions (converted: {100*y_gas[-1]:.2f}%)")
    ax1.grid(True)
    ax1.legend()

    err = (E - E[0])/E[0]
    ax2.plot(t, err, "g-")
    ax2.set_xlabel("Time [ms]")
    ax2.set_ylabel(r"$(E - E_0)/E_0$ [-]")
    ax2.set_title(rf"(b) Total energy, max $|\Delta E/E_0| = {np.max(np.abs(err)):.1e}$")
    ax2.grid(True)

    p = dp[:, 1]/1e5
    ax3.plot(t, p, color="darkorange")
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("Pressure [bar]")
    ax3.set_title(rf"(c) Pressure, max $|\Delta p/p_0| = {np.max(np.abs(p/p[0] - 1)):.1e}$")
    ax3.ticklabel_format(useOffset=False)
    ax3.grid(True)

    ax4.plot(t, T_comp, color="purple", label="Products")
    ax4.plot(t, T_gas, color="teal", ls="--", label="Gas")
    ax4.set_xlabel("Time [ms]")
    ax4.set_ylabel("Temperature [K]")
    ax4.set_title("(d) Phase temperatures")
    ax4.ticklabel_format(useOffset=False)
    ax4.grid(True)
    ax4.legend()
    save(fig, "stage1_0d_phase_transfer")


def fig_stage1_al_kinetics():
    """Reference particle model: heating, profiles, burning (Beckstead)."""
    import verify_al_particle_kinetics as al
    d0 = 20e-6
    T_gas = 2800.0
    g = al.simulate_goodman_particle_heating(d0, T_gas=T_gas)
    lu = al.simulate_lumped_particle_heating(d0, T_gas=T_gas)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.3))
    ax1.plot(g["times"]*1e6, g["Ts"], "crimson", label="Surface $T_s$")
    ax1.plot(g["times"]*1e6, g["Tavg"], "forestgreen", ls="-.", label="Mean")
    ax1.plot(g["times"]*1e6, g["Tc"], "royalblue", ls="--", label="Centre $T_c$")
    ax1.plot(lu["times"]*1e6, lu["T_lump"], "dimgray", ls=":", label="Uniform (lumped)")
    ax1.axhline(al.T_MELT_AL, color="navy", ls="--", alpha=0.5)
    ax1.axhline(al.T_IGN_AL, color="darkorange", ls="--", alpha=0.6)
    ax1.set_xlabel(r"Time [$\mu$s]")
    ax1.set_ylabel("Temperature [K]")
    ax1.set_title(rf"(a) {d0*1e6:.0f} $\mu$m particle, ignition at "
                  rf"{g['t_ign_surf']*1e6:.1f} / {lu['t_ign']*1e6:.1f} $\mu$s")
    ax1.grid(True)
    ax1.legend(loc="lower right")

    for key, colour in (("delta_04", "darkviolet"), ("delta_08", "royalblue"),
                        ("penetration", "darkcyan"), ("melt", "darkorange"),
                        ("ign", "crimson")):
        if key in g["snaps"]:
            snap = g["snaps"][key]
            rn, T = al.reconstruct_radial_profile(snap, g["R"])
            ax2.plot(rn, T, color=colour, label=rf"$t = {snap[0]*1e6:.1f}\,\mu$s ({key.replace('_', ' ')})")
    ax2.set_xlabel("$r/R$ [-]")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_title("(b) Internal temperature profiles")
    ax2.grid(True)
    ax2.legend(loc="upper left")

    t_b = al.eval_beckstead_burn_time(d0, 0.21, 0.0, 0.0, T0=T_gas)
    t = np.linspace(0, t_b, 200)
    d = d0*np.maximum(0.0, 1 - t/t_b)**(1/al.BECKSTEAD_N)
    ax3.plot(t*1e3, d*1e6, "navy", label="Diameter ($d^{1.8}$ law)")
    ax3r = ax3.twinx()
    ax3r.plot(t*1e3, 31.0*(1 - (d/d0)**3), "forestgreen", ls="--",
              label="Heat released (31 MJ/kg Al to Al$_2$O$_3$(s))")
    ax3.set_xlabel("Time after ignition [ms]")
    ax3.set_ylabel(r"Diameter [$\mu$m]")
    ax3r.set_ylabel("Heat released [MJ/kg Al]")
    ax3.set_title(rf"(c) Beckstead burn time in air: {t_b*1e3:.2f} ms")
    ax3.grid(True)
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3r.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, loc="center left")
    save(fig, "stage1_al_particle_kinetics")


# ---------------------------------------------------------------------------
# Stage 2: 1D canonical
# ---------------------------------------------------------------------------

def fig_stage2_high_mach():
    """Strong shock in air: blastFoam (NASA-9 air in local equilibrium)
    against the frozen and equilibrium jumps for the same species data."""
    import verify_shock_dissociating_air as vs
    case = result("tests/1D/shock_dissociating_air")
    cfd = np.load(result(case, "shock_fields.npz"))
    vs.CASE_DIR = case
    S = vs.air_species(case)
    Y1 = vs.air_mass_fractions(case)
    machs = np.linspace(2.0, 16.0, 29)
    a1 = vs.sound_speed(S, Y1, 300.0, 101325.0)
    p_ratio = 1 + 2.8/2.4*(machs**2 - 1)
    rho_id = 2.4*machs**2/(0.4*machs**2 + 2)
    T_id = 300.0*p_ratio/rho_id
    fr = [vs.jump(M*a1, 101325.0, 300.0, S, Y1, False) for M in machs]
    eq = [vs.jump(M*a1, 101325.0, 300.0, S, Y1, True) for M in machs]
    Ma, T2 = float(cfd["Ma_s"]), float(cfd["T2"])
    rr = float(cfd["rho2"]/cfd["rho1"])

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.5, 7.6))
    ax1.plot(machs, T_id, "b:", label=r"Ideal gas, $\gamma = 1.4$")
    ax1.plot(machs, [s["T2"] for s in fr], "b--", label="NASA-9 air, frozen")
    ax1.plot(machs, [s["T2"] for s in eq], "r-", label="NASA-9 air (N$_2$ O$_2$ NO N O), equilibrium")
    ax1.plot([Ma], [T2], "k*", ms=12, label=rf"blastFoam, Ma {Ma:.1f}")
    ax1.set_xlabel("Shock Mach number")
    ax1.set_ylabel("Post-shock temperature [K]")
    ax1.set_title("(a) Post-shock temperature")
    ax1.grid(True)
    ax1.legend(loc="upper left")

    ax2.plot(machs, rho_id, "b:", label=r"Ideal gas, $\gamma = 1.4$")
    ax2.plot(machs, [s["rho2"]/s["rho1"] for s in fr], "b--", label="NASA-9 air, frozen")
    ax2.plot(machs, [s["rho2"]/s["rho1"] for s in eq], "r-", label="NASA-9 air, equilibrium")
    ax2.plot([Ma], [rr], "k*", ms=12, label="blastFoam")
    ax2.set_xlabel("Shock Mach number")
    ax2.set_ylabel(r"$\rho_2/\rho_1$ [-]")
    ax2.set_title("(b) Density ratio")
    ax2.grid(True)
    ax2.legend(loc="upper left")

    t_us = float(cfd["t"])*1e6
    for ax, f, lab in ((ax3, "T", "Temperature [K]"), (ax4, "rho", r"Density [kg/m$^3$]")):
        ax.plot(cfd["x"], cfd[f], "k-", label="blastFoam")
        ax.axvline(float(cfd["x_shock"]), color="gray", ls=":")
        ax.set_xlabel("$x$ [m]")
        ax.set_ylabel(lab)
        ax.grid(True)
    xs = float(cfd["x_shock"])
    x_c = xs - (float(cfd["Us"]) - float(cfd["u2"]))*float(cfd["t"])
    ax3.hlines(float(cfd["T2_eq"]), x_c, xs, color="r", lw=2, label="Equilibrium $T_2$")
    ax4.hlines(float(cfd["rho2_eq"]), x_c, xs, color="r", lw=2, label=r"Equilibrium $\rho_2$")
    ax4.set_yscale("log")
    ax3.legend(loc="lower left")
    ax4.legend(loc="lower left")
    ax3.set_title(rf"(c) Temperature at {t_us:.0f} $\mu$s")
    ax4.set_title(rf"(d) Density at {t_us:.0f} $\mu$s")
    save(fig, "stage2_high_mach_shock")


def fig_stage2_multiphase_advection():
    """Contact advection: four flux schemes; van Leer vs Minmod."""
    import verify_multiphase_advection as va
    d = np.load(result("tests/1D/multiphase_advection/advection_fields.npz"))
    n = len(d["p_HLLC_Minmod"])
    x = (np.arange(n) + 0.5)*va.DX
    p0 = 100.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.plot(x, np.where(x < va.X_EXACT, 1.0, 0.0), "k--", lw=1.2, label="Exact")
    for key, lab in (("HLLC_Minmod", "HLLC"), ("AUSMplus_Minmod", "AUSM+"),
                     ("Kurganov_Minmod", "Kurganov"), ("Tadmor_vanLeer", "Tadmor (van Leer)")):
        if f"alpha_{key}" in d:
            ax1.plot(x, d[f"alpha_{key}"], label=lab)
    ax1.set_xlim(va.X_EXACT - 0.15, va.X_EXACT + 0.15)
    ax1.set_xlabel("$x$ [m]")
    ax1.set_ylabel(r"$\alpha_1$ [-]")
    ax1.set_title(rf"(a) Volume fraction at {va.T_END*1e3:.1f} ms (Minmod unless stated)")
    ax1.grid(True)
    ax1.legend()

    for key, lab, c in (("HLLC_vanLeer", "HLLC, van Leer", "r"), ("HLLC_Minmod", "HLLC, Minmod", "g")):
        p = d[f"p_{key}"]/1e5
        ax2.plot(x, p, c, label=rf"{lab} (max {100*np.max(np.abs(p/p0 - 1)):.2f}%)")
    ax2.axhline(p0, color="k", ls="--", lw=1.0)
    ax2.set_xlim(va.X_EXACT - 0.15, va.X_EXACT + 0.15)
    ax2.set_xlabel("$x$ [m]")
    ax2.set_ylabel("Pressure [bar]")
    ax2.set_title("(b) Pressure at the contact")
    ax2.grid(True)
    ax2.legend()
    save(fig, "stage2_multiphase_advection")


def fig_stage2_gunib1d():
    """Houim 2 granular shock tube: blastEulerFoam vs gunIB1D."""
    d = np.load(result("tests/1D/gunib1d_houim2/houim2_profiles.npz"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.plot(d["x_gun"]*1e3, d["p_gun"]/1e5, "k-", label="gunIB1D")
    ax1.plot(d["x_bf"]*1e3, d["p_bf"]/1e5, "r--", label="blastEulerFoam")
    ax1.set_xlabel("$x$ [mm]")
    ax1.set_ylabel("Gas pressure [bar]")
    ax1.set_title(r"(a) Houim 2 at 0.184 ms: pressure")
    ax1.grid(True)
    ax1.legend()
    ax2.plot(d["x_gun"]*1e3, d["u_g_gun"], "k-", label="Gas, gunIB1D")
    ax2.plot(d["x_bf"]*1e3, d["u_g_bf"], "r--", label="Gas, blastEulerFoam")
    # Particle velocity only where the bed is: gunIB1D keeps a residual
    # alpha_s ~ 1e-5 (2% of the bed) with the gas velocity where the particles
    # have left, blastEulerFoam sets the velocity of a residual phase to zero
    alpha_min = 0.05*d["alpha_s_gun"][-1]
    for k, style, code in (("gun", "b-", "gunIB1D"), ("bf", "c--", "blastEulerFoam")):
        u_s = np.where(d["alpha_s_" + k] > alpha_min, d["u_s_" + k], np.nan)
        ax2.plot(d["x_" + k]*1e3, u_s, style, label=f"Particles, {code}")
    ax2.set_xlabel("$x$ [mm]")
    ax2.set_ylabel("Velocity [m/s]")
    ax2.set_title(r"(b) Gas and particle velocities (particles where $\alpha_s > 0.05\,\alpha_{s,0}$)")
    ax2.grid(True)
    ax2.legend()
    save(fig, "stage2_gunib1d_benchmarks")


def fig_stage2_planar_znd():
    """Planar Comp-B detonation."""
    from verify_planar_znd import get_eudetcode_reference
    case = result("tests/1D/planar_znd")
    ref = get_eudetcode_reference()
    times = sorted((d for d in os.listdir(case) if d.replace(".", "").replace("e-", "").isdigit()
                    and float(d) > 0), key=float)
    sel = [min(times, key=lambda d: abs(float(d) - t)) for t in (5e-6, 10e-6, 15e-6)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    for d, c in zip(sel, ("darkgreen", "navy", "crimson")):
        p = read_field(os.path.join(case, d, "p"))
        x = (np.arange(len(p)) + 0.5)*0.2/len(p)*1e3
        ax1.plot(x, p/1e9, color=c, label=rf"$t = {float(d)*1e6:.1f}\,\mu$s")
    ax1.axhline(ref["P_CJ"]/1e9, color="k", ls="--", lw=1.0,
                label=rf"$p_\mathrm{{CJ}}$ (EuDetCode) = {ref['P_CJ']/1e9:.1f} GPa")
    ax1.set_xlabel("Distance from the closed end [mm]")
    ax1.set_ylabel("Pressure [GPa]")
    ax1.set_title("(a) Pressure")
    ax1.grid(True)
    ax1.legend()

    d = sel[1]
    u = read_field(os.path.join(case, d, "U"))
    rho = read_field(os.path.join(case, d, "rho"))
    x = (np.arange(len(u)) + 0.5)*0.2/len(u)*1e3
    ax2.plot(x, u, "navy", label="Velocity")
    ax2.set_xlabel("Distance from the closed end [mm]")
    ax2.set_ylabel("Velocity [m/s]")
    ax2r = ax2.twinx()
    ax2r.plot(x, rho/ref["rho0"], "crimson", ls="--", label=r"$\rho/\rho_0$")
    ax2r.set_ylabel(r"$\rho/\rho_0$ [-]")
    ax2.set_title(rf"(b) Velocity and density at {float(d)*1e6:.1f} $\mu$s")
    ax2.grid(True)
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc="upper right")
    save(fig, "stage2_planar_znd_detonation")


# ---------------------------------------------------------------------------
# Stage 3: 1D spherical
# ---------------------------------------------------------------------------

def fig_stage3_spherical_blast():
    """1 kg C4 spherical blast against the Kinney-Graham TNT correlations."""
    from verify_kingery_bulmash import (
        kinney_graham_overpressure, kinney_graham_positive_impulse, P0_ATM)
    case = result("tests/1D_spherical/spherical_blast")
    pk = np.loadtxt(result(case, "postProcessing/radialPeakPressure/0.005/p_temporalMax.dat"))
    im = np.loadtxt(result(case, "postProcessing/radialImpulse/0.005/p_temporalIntegral.dat"))
    probes = np.loadtxt(result(case, "postProcessing/pressureProbes/0/p"), comments="#")

    Z = np.geomspace(0.4, 5.0, 100)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.3))
    m = (pk[:, 1] >= 0.4) & (pk[:, 5] - P0_ATM > 5e3)
    ax1.loglog(Z, kinney_graham_overpressure(Z)/1e5, "k--", label="Kinney–Graham (TNT)")
    ax1.loglog(pk[m, 1], (pk[m, 5] - P0_ATM)/1e5, "crimson", label="blastFoam (300 cells)")
    ax1.set_xlabel(r"$Z = R/W^{1/3}$ [m/kg$^{1/3}$]")
    ax1.set_ylabel("Peak overpressure [bar]")
    ax1.set_title("(a) Peak overpressure")
    ax1.grid(True, which="both")
    ax1.legend()

    m = (im[:, 1] >= 0.4) & (im[:, 5] > 1.0)
    ax2.plot(Z, kinney_graham_positive_impulse(Z), "k--", label="Kinney–Graham (TNT)")
    ax2.plot(im[m, 1], im[m, 5], "royalblue", label="blastFoam (300 cells)")
    ax2.set_xlabel(r"$Z = R/W^{1/3}$ [m/kg$^{1/3}$]")
    ax2.set_ylabel(r"Positive impulse [Pa$\,$s]")
    ax2.set_title("(b) Positive impulse")
    ax2.grid(True)
    ax2.legend()

    for i, (r, c) in enumerate(zip((0.5, 1.0, 1.5, 2.0), ("navy", "crimson", "forestgreen", "darkorange"))):
        ax3.plot(probes[:, 0]*1e3, probes[:, i + 1]/1e5, color=c, label=f"$r = {r:.1f}$ m")
    ax3.set_xlim(0, 4)
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel("Pressure [bar]")
    ax3.set_title("(c) Pressure histories")
    ax3.grid(True)
    ax3.legend()
    save(fig, "stage3_spherical_blast_kb")


def fig_stage3_spherical_multiphase():
    """Spherical Comp-B/Al charge: inert, reacting, Al combustion."""
    cases = [("ChargeCompBAl_1D_spherical", "Inert Al", "-"),
             ("ChargeCompBAl_1D_spherical_reacting", "Gas chemistry", "--"),
             ("ChargeCompBAl_1D_spherical_alCombustion", "Al burning", ":")]
    dirs = {c: result("tests/1D_spherical", c) for c, _, _ in cases}
    probes_r = [0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12.5, 9.2))
    pmax = 0
    for c, lab, ls in cases:
        pk = np.loadtxt(result(dirs[c], "postProcessing/radialPeakPressure/4e-05/p_temporalMax.dat"))
        ax1.plot(pk[:, 1]*1e3, pk[:, 5]/1e9, ls=ls, label=lab)
        pmax = max(pmax, pk[:, 5].max()/1e9)
    pk = np.loadtxt(result(dirs[cases[0][0]], "postProcessing/radialPeakPressure/4e-05/p_temporalMax.dat"))
    i_pk = np.argmax(pk[:, 5])
    i_bo = np.searchsorted(pk[:, 1], 0.05)
    ax1.plot(pk[i_pk, 1]*1e3, pk[i_pk, 5]/1e9, "ko", mfc="none",
             label=f"Peak {pk[i_pk, 5]/1e9:.2f} GPa")
    ax1.plot(pk[i_bo, 1]*1e3, pk[i_bo, 5]/1e9, "ks", mfc="none",
             label=f"Breakout {pk[i_bo, 5]/1e9:.2f} GPa")
    ax1.axvline(50, color="gray", ls=":")
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 1.08*pmax)
    ax1.set_xlabel("$r$ [mm]")
    ax1.set_ylabel("Peak pressure [GPa]")
    ax1.set_title("(a) Peak pressure (temporalMax)")
    ax1.grid(True)
    ax1.legend()

    colours = ("navy", "darkgreen", "darkorange", "purple")
    for c, lab, ls in cases:
        pr = np.loadtxt(result(dirs[c], "postProcessing/pressureProbes/0/p"), comments="#")
        for k, idx in enumerate((1, 3, 4, 5)):
            ax2.plot(pr[:, 0]*1e6, pr[:, idx + 1]/1e9, color=colours[k], ls=ls,
                     label=f"$r = {probes_r[idx]*1e3:.0f}$ mm" if c == cases[0][0] else None)
    ax2.set_xlabel(r"Time [$\mu$s]")
    ax2.set_ylabel("Pressure [GPa]")
    ax2.set_title("(b) Pressure histories (line style as in (a))")
    ax2.grid(True)
    ax2.legend()

    for t, col in zip(("5e-06", "1e-05", "2e-05", "4e-05"), colours):
        for c, lab, ls in cases[:2]:
            s = np.loadtxt(result(dirs[c], "postProcessing/sampleLine", t, "lineX.xy"))
            ax3.plot(s[:, 0]*1e3, s[:, 1]/1e9, color=col, ls=ls,
                     label=rf"{float(t)*1e6:.0f} $\mu$s" if c == cases[0][0] else None)
    ax3.axvline(50, color="gray", ls=":")
    ax3.set_xlim(0, 200)
    ax3.set_xlabel("$r$ [mm]")
    ax3.set_ylabel("Pressure [GPa]")
    ax3.set_title("(c) Pressure profiles, inert (solid) and gas chemistry (dashed)")
    ax3.grid(True)
    ax3.legend()

    alc = dirs[cases[2][0]]
    for t, col in zip(("5e-06", "1e-05", "2e-05", "4e-05"), colours):
        a = read_field(result(alc, t, "alpha.Al"))
        d = np.broadcast_to(read_field(result(alc, t, "d.Al")), a.shape)
        r = (np.arange(len(a)) + 0.5)*0.2/len(a)*1e3
        dec = np.where(a >= 1e-4, np.maximum(1 - d/30e-6, 1e-12), np.nan)
        ax4.semilogy(r, dec, color=col, label=rf"{float(t)*1e6:.0f} $\mu$s")
    ax4.set_xlim(0, 200)
    ax4.set_ylim(1e-12, 1e-2)
    ax4.axvline(50, color="gray", ls=":")
    ax4.set_xlabel("$r$ [mm]")
    ax4.set_ylabel(r"Relative diameter decrease $1 - d/d_0$ ($\alpha_\mathrm{Al} \geq 10^{-4}$)")
    ax4r = ax4.twinx()
    a40 = read_field(result(alc, "4e-05", "alpha.Al"))
    g40 = read_field(result(alc, "4e-05", "alpha.gas"))
    r = (np.arange(len(a40)) + 0.5)*0.2/len(a40)*1e3
    ax4r.plot(r, a40, "crimson", lw=1.0, ls="--", label=r"$\alpha_\mathrm{Al}$, 40 $\mu$s")
    ax4r.plot(r, g40, "saddlebrown", lw=1.0, ls="-.", label=r"$\alpha_\mathrm{gas}$, 40 $\mu$s")
    ax4r.set_ylabel("Volume fraction [-]")
    ax4r.set_ylim(-0.02, 1.02)
    ax4.set_title(r"(d) Al burning case: diameter decrease ($d_0 = 30\,\mu$m) and volume fractions")
    ax4.grid(True, which="both")
    h1, l1 = ax4.get_legend_handles_labels()
    h2, l2 = ax4r.get_legend_handles_labels()
    ax4.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=7.5)
    for ax in (ax1, ax2, ax3):
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    save(fig, "stage3_spherical_multiphase")


def fig_stage1_al_solver():
    """Al particle models in blastEulerFoam (Phase D): 0D tests vs references."""
    import verify_al_melting as vm
    import verify_particle_drag as vd
    import verify_particle_heating as vh
    import verify_al_burning as vb

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8))

    # (a) Melting plateau
    data = np.loadtxt(result("tests/0D/alMelting_0D/postProcessing/temperature/0/volFieldValue.dat"))
    al, t_ref, T_al_ref, T_g_ref = vm.reference(data[-1, 0])
    ax1.plot(t_ref*1e3, T_al_ref, "k-", lw=1.0, label="Reference (ODE)")
    ax1.plot(t_ref*1e3, T_g_ref, "k-", lw=1.0)
    k = max(1, len(data)//60)
    ax1.plot(data[::k, 0]*1e3, data[::k, 1], "o", ms=3.5, mfc="none", color="crimson", label="Al, blastEulerFoam")
    ax1.plot(data[::k, 0]*1e3, data[::k, 2], "s", ms=3.5, mfc="none", color="navy", label="Gas, blastEulerFoam")
    ax1.axhline(al.Tm, color="gray", ls=":", lw=1.0)
    ax1.set_xlabel("$t$ [ms]")
    ax1.set_ylabel("Temperature [K]")
    ax1.set_title(r"(a) Melting: Al (30 $\mu$m, $\alpha$ 1e-4) in gas at 2500 K")
    ax1.grid(True)
    ax1.legend(loc="center right")

    # (b) Drag: slip relaxation
    for (name, d, p, T, U0, alpha, dt, t_end), c in zip(vd.RUNS, ("crimson", "navy")):
        path = result("tests/0D/particle_drag", name, "postProcessing/velocity/0/volFieldValue.dat")
        text = open(path).read().replace("(", " ").replace(")", " ")
        dd = np.loadtxt(text.splitlines(), comments="#", ndmin=2)
        t_ref, ur_ref, Re0, Ma0, Kn0 = vd.reference(d, p, T, U0, alpha, t_end)
        ax2.plot(t_ref*1e3, ur_ref/U0, "-", color=c, lw=1.0,
                 label=rf"{name}: $d$ {d*1e6:g} $\mu$m, $Re_0$ {Re0:.0f}, $Ma_0$ {Ma0:.2f}, $Kn_0$ {Kn0:.2g}")
        k = max(1, len(dd)//40)
        ax2.plot(dd[::k, 0]*1e3, (dd[::k, 1] - dd[::k, 4])/U0, "o", ms=3, mfc="none", color=c)
    ax2.set_xscale("log")
    ax2.set_xlabel("$t$ [ms]")
    ax2.set_ylabel(r"Slip velocity / $U_0$")
    ax2.set_title("(b) Loth drag: lines reference, symbols blastEulerFoam")
    ax2.grid(True)
    ax2.legend(fontsize=7.5)

    # (c) Heat transfer with the Kavanau Nusselt number
    for (name, d, p, T, Tp, U0, alpha, dt, t_end), c in zip(vh.RUNS, ("crimson", "navy")):
        path = result("tests/0D/particle_heating", name, "postProcessing/state/0/volFieldValue.dat")
        text = open(path).read().replace("(", " ").replace(")", " ")
        dd = np.loadtxt(text.splitlines(), comments="#", ndmin=2)
        t_ref, T_ref, MaByRe, Nu0 = vh.reference(d, p, T, Tp, U0, alpha, t_end)
        ax3.plot(t_ref*1e3, T_ref, "-", color=c, lw=1.0,
                 label=rf"{name}: $d$ {d*1e6:g} $\mu$m, {p/1e5:g} bar, $U_0$ {U0:g} m/s, $Nu_0$ {Nu0:.2f}")
        k = max(1, len(dd)//40)
        ax3.plot(dd[::k, 0]*1e3, dd[::k, 7], "o", ms=3, mfc="none", color=c)
    ax3.set_xscale("log")
    ax3.set_xlabel("$t$ [ms]")
    ax3.set_ylabel("Al temperature [K]")
    ax3.set_title("(c) Kavanau heat transfer: lines reference, symbols blastEulerFoam")
    ax3.grid(True)
    ax3.legend(fontsize=7.5)

    # (d) Beckstead burning
    x = np.linspace(0, 1, 50)
    ax4.plot(x, 1 - x, "k-", lw=1.0, label="Beckstead $d^n$ law")
    for j, (d0, c, m) in enumerate(zip(vb.D0S, ("crimson", "darkgreen", "navy"), ("o", "s", "^"))):
        tb = vb.burn_time(d0)
        dd = np.loadtxt(result("tests/0D/al_burning", f"d{d0*1e6:g}um", "postProcessing/state/0/volFieldValue.dat"))
        up = np.nonzero(np.diff(dd[:, 1]) > 0)[0]
        n = up[0] + 1 if len(up) else len(dd)
        k = max(1, n//24)
        ax4.plot(dd[j*k//3:n:k, 0]/tb, (dd[j*k//3:n:k, 1]/d0)**vb.N, m, ms=4, mfc="none", color=c,
                 label=rf"$d_0$ = {d0*1e6:g} $\mu$m ($t_b$ = {tb*1e3:.3g} ms)")
    ax4.set_xlabel(r"$t/t_b$")
    ax4.set_ylabel(r"$(d/d_0)^{1.8}$")
    ax4.set_title("(d) Ignited Al in air at 2500 K, 1 atm")
    ax4.grid(True)
    ax4.legend()
    save(fig, "stage1_al_solver")


def fig_stage3_blast_kb():
    """1 kg TNT and Comp-B free-air blast vs Kingery-Bulmash (E1)."""
    import kingery_bulmash as kb
    from verify_kingery_bulmash import kinney_graham_overpressure, kinney_graham_positive_impulse

    def radial(case):
        pp = result("tests/1D_spherical", case, "postProcessing")
        t = sorted(os.listdir(os.path.join(pp, "radialPeakPressure")), key=float)[-1]
        pm = np.loadtxt(os.path.join(pp, "radialPeakPressure", t, "p_temporalMax.dat"))
        im = np.loadtxt(os.path.join(pp, "radialImpulse", t, "p_temporalIntegral.dat"))
        return pm[:, 1], pm[:, 5] - 101325.0, im[:, 5]

    Z = np.geomspace(0.3, 9.5, 200)
    p_kb, i_kb, _ = kb.free_air_arrays(Z, 1.0)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4.4))
    ax1.loglog(Z, p_kb/1e3, "k-", lw=1.2, label="Kingery-Bulmash (TNT)")
    ax1.loglog(Z, kinney_graham_overpressure(Z)/1e3, "k:", lw=1.0, label="Kinney-Graham (TNT)")
    ax2.loglog(Z, i_kb, "k-", lw=1.2, label="Kingery-Bulmash (TNT)")
    ax2.loglog(Z, kinney_graham_positive_impulse(Z), "k:", lw=1.0, label="Kinney-Graham (TNT)")
    runs = [("tnt_blast", "TNT, 4 mm", "crimson", "-"),
            ("tnt_blast_n5000", "TNT, 2 mm", "darkorange", "-"),
            ("tnt_blast_near1mm", "TNT, 1 mm to 0.3 m, 2 mm beyond", "olive", ":"),
            ("tnt_blast_n10000", "TNT, 1 mm", "goldenrod", "-"),
            ("compB_blast", "Comp-B, 4 mm", "navy", "--"),
            ("compB_blast_n5000", "Comp-B, 2 mm", "royalblue", "--")]
    for case, label, c, ls in runs:
        if not os.path.isdir(os.path.join(RESULTS, "tests/1D_spherical", case)):
            continue
        r, dp, imp = radial(case)
        sel = (r > 0.3) & (r < 9.5)
        ax1.loglog(r[sel], dp[sel]/1e3, ls, color=c, lw=1.1, label=label)
        ax2.loglog(r[sel], imp[sel], ls, color=c, lw=1.1, label=label)
        eq_p, eq_i = (1.11, 0.98) if case.startswith("compB") else (1.0, 1.0)
        pk = np.array([kb.free_air(x, eq_p)[0] for x in r[sel]])
        ik = np.array([kb.free_air(x, eq_i)[1] for x in r[sel]])
        ax3.plot(r[sel], 100*(dp[sel]/pk - 1), ls, color=c, lw=1.1, label=label + " (overpressure)")
        if case in ("tnt_blast_n10000", "compB_blast_n5000", "tnt_blast", "compB_blast"):
            ax3.plot(r[sel], 100*(imp[sel]/ik - 1), ls, color=c, lw=0.8, alpha=0.5)
    ax1.set_xlabel(r"$Z$ [m/kg$^{1/3}$] (1 kg charges)")
    ax1.set_ylabel("Peak incident overpressure [kPa]")
    ax1.set_title("(a) Peak overpressure")
    ax2.set_xlabel(r"$Z$ [m/kg$^{1/3}$]")
    ax2.set_ylabel("Positive impulse [Pa s]")
    ax2.set_title("(b) Positive impulse")
    ax3.axhspan(-5, 5, color="green", alpha=0.08)
    ax3.axhline(0, color="k", lw=0.8)
    ax3.set_xscale("log")
    ax3.set_xlabel(r"$Z$ [m/kg$^{1/3}$]")
    ax3.set_ylabel("Deviation from Kingery-Bulmash [%]")
    ax3.set_title("(c) Deviation (thin: impulse; Comp-B with TNT equiv.)")
    ax3.set_ylim(-40, 40)
    for ax in (ax1, ax2, ax3):
        ax.grid(True, which="both")
    ax1.legend(fontsize=7.5)
    ax3.legend(fontsize=7)
    save(fig, "stage3_blast_kb")


def fig_stage3_e2_ignition():
    """Comp-B/Al on the 2 m case (plan step E2): three Al ignition
    temperatures against the inert and products-only runs. The runs are made
    with tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m/setup_e2.py
    (burning to oxide, transported ignited fraction) and read from
    $E2_RESULTS_DIR/{inert,products,Tign2050,Tign932,Tign0}."""
    from compare_e2 import first_phase
    base = os.environ.get("E2_RESULTS_DIR")
    if not base:
        raise FileNotFoundError("E2_RESULTS_DIR not set (directory of the E2 runs)")
    runs = [("Tign2050", r"$T_\mathrm{ign}$ = 2050 K (oxide melting)", "crimson", "-"),
            ("Tign932", r"$T_\mathrm{ign}$ = 933 K (Al melting, ruptured shell)", "darkorange", "--"),
            ("Tign0", r"$T_\mathrm{ign}$ = 0 (burning from the start)", "navy", ":")]
    P0 = 101325.0

    def path(c, *parts):
        q = os.path.join(base, c, *parts)
        if not os.path.exists(q):
            raise FileNotFoundError(q)
        return q

    def last(c, sub):
        return sorted(os.listdir(path(c, "postProcessing", sub)), key=float)[-1]

    def peak(c):
        pm = np.loadtxt(path(c, "postProcessing/radialPeakPressure", last(c, "radialPeakPressure"), "p_temporalMax.dat"))
        return pm[:, 1], pm[:, 5] - P0

    def al_history(c):
        r = peak(c)[0]
        f = np.r_[0, 0.5*(r[1:] + r[:-1]), r[-1] + 0.5*(r[-1] - r[-2])]
        V = np.diff(f**3)
        times = sorted((d for d in os.listdir(path(c)) if re.fullmatch(r"[0-9.e-]+", d) and float(d) > 0), key=float)
        t, M, T, I = [], [], [], []
        for d in times:
            a = read_field(path(c, d, "alpha.Al"))
            rho = np.broadcast_to(read_field(path(c, d, "rho.Al")), a.shape)
            Tp = np.broadcast_to(read_field(path(c, d, "T.Al")), a.shape)
            ig = np.broadcast_to(read_field(path(c, d, "ignited.Al")), a.shape)
            m = a*rho*V
            t.append(float(d)); M.append(m.sum()); T.append(np.sum(m*Tp)/m.sum())
            I.append(np.sum(m*ig)/m.sum())
        return np.array(t), np.array(M), np.array(T), np.array(I)

    def probes(c):
        p = path(c, "postProcessing/pressureProbes/0/p")
        radii = [float(re.search(r"\(([^ ]+)", line).group(1)) for line in open(p) if line.startswith("# Probe")]
        return np.array(radii), np.loadtxt(p)

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    (ax1, ax2), (ax3, ax4) = axs
    M0 = None
    for c, label, col, ls in runs:
        t, M, T, I = al_history(c)
        M0 = M0 or M[0]
        ax1.plot(t*1e3, 100*(1 - M/M0), ls, color=col, label=label)
        ax1.plot(t*1e3, 100*I, ls, color=col, lw=0.8, alpha=0.6)
        ax2.plot(t*1e3, T, ls, color=col, label=label)
    ax1.set_xlabel("t [ms]"); ax1.set_ylabel("Al burnt [% of initial mass]")
    ax1.set_title("(a) Aluminium burnt (thin: ignited, % of remaining mass)")
    ax2.axhline(933.47, color="grey", lw=0.8); ax2.axhline(2050, color="grey", lw=0.8, ls="--")
    ax2.set_xlabel("t [ms]"); ax2.set_ylabel("Mass-mean Al temperature [K]")
    ax2.set_title("(b) Aluminium temperature")

    r_in, dp_in = peak("inert")
    sel = (r_in > 0.06) & (r_in < 1.9)
    R, d_in = probes("inert")
    imp_in = np.array([first_phase(d_in[:, 0], d_in[:, j + 1] - P0)[2] for j in range(len(R))])
    for c, label, col, ls in [("products", "products transfer, no Al burning", "seagreen", "-.")] + runs:
        r, dp = peak(c)
        ax3.plot(r[sel], 100*(dp[sel]/dp_in[sel] - 1), ls, color=col, label=label)
        Rc, d = probes(c)
        res = [first_phase(d[:, 0], d[:, j + 1] - P0) for j in range(len(Rc))]
        imp = np.array([x[2] for x in res]); done = np.array([x[3] for x in res])
        ok = (Rc >= 0.08) & done
        ax4.plot(Rc[ok], 100*(imp[ok]/imp_in[ok] - 1), ls, marker="o", ms=3.5, color=col, label=label)
    ax3.axhline(0, color="k", lw=0.8)
    ax3.set_xlabel("r [m]"); ax3.set_ylabel("Peak overpressure vs inert [%]")
    ax3.set_title("(c) Peak overpressure")
    ax4.axhline(0, color="k", lw=0.8)
    ax4.set_xlabel("r [m]"); ax4.set_ylabel("First positive-phase impulse vs inert [%]")
    ax4.set_title("(d) Impulse (gauges with a completed positive phase)")
    for ax in axs.flat:
        ax.grid(True)
    ax1.legend(); ax3.legend(fontsize=7.5)
    save(fig, "stage3_e2_ignition")


def fig_stage5_g1_axisymmetric():
    """2D axisymmetric repeat of the spherical Comp-B/Al case (plan step G1):
    peak overpressure and arrival against the angle from the symmetry plane at
    r = 0.10 and 0.15 m, HLLC and hybrid HLLC/HLL, with the 1D values of the
    same flux. Read from $G1_RESULTS_DIR/{hllc,hybrid} (runs of
    tests/2D_axi/ChargeCompBAl_2D_axi_spherical, hybridHLL no / yes) and
    $G1_RESULTS_DIR/{1d_hllc,1d_hybrid} (tests/1D_spherical/
    ChargeCompBAl_1D_spherical)."""
    base = os.environ.get("G1_RESULTS_DIR")
    if not base:
        raise FileNotFoundError("G1_RESULTS_DIR not set (directory of the G1 runs)")
    P0 = 101325.0
    angles = np.arange(0, 91, 5)
    radii = (0.10, 0.15)
    i1d = {0.10: 6, 0.15: 7}

    def probes(c, name):
        root = os.path.join(base, c, "postProcessing", name)
        if not os.path.isdir(root):
            raise FileNotFoundError(root)
        a = np.vstack([np.loadtxt(os.path.join(root, t, "p"), comments="#")
                       for t in sorted(os.listdir(root), key=float)])
        return a[:, 0], a[:, 1:]

    def arrival(t, p):
        level = P0 + 0.5*(p.max() - P0)
        i = int(np.argmax(p > level))
        return t[i - 1] + (level - p[i - 1])/(p[i] - p[i - 1])*(t[i] - t[i - 1])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    for c, oned, label, colour in (("hllc", "1d_hllc", "HLLC", "crimson"),
                                   ("hybrid", "1d_hybrid", "hybrid HLLC/HLL", "navy")):
        t, P = probes(c, "angularProbes")
        t1, P1 = probes(oned, "pressureProbes")
        for k, (r, ls) in enumerate(zip(radii, ("-", "--"))):
            pr = P[:, k*len(angles):(k + 1)*len(angles)]
            pk = (pr.max(axis=0) - P0)/1e6
            ta = np.array([arrival(t, pr[:, i]) for i in range(len(angles))])
            pk1 = (P1[:, i1d[r]].max() - P0)/1e6
            ta1 = arrival(t1, P1[:, i1d[r]])
            ax1.plot(angles, pk/pk1, ls, color=colour, marker="o", ms=3,
                     label=f"{label}, r = {r:.2f} m")
            ax2.plot(angles, (ta/ta1 - 1)*100, ls, color=colour, marker="o",
                     ms=3, label=f"{label}, r = {r:.2f} m")
    ax1.axhline(1, color="k", lw=0.8)
    ax1.set_xlabel("Angle from the symmetry plane [deg]")
    ax1.set_ylabel("Peak overpressure / 1D (same flux)")
    ax1.set_title("(a) Peak overpressure")
    ax2.axhline(0, color="k", lw=0.8)
    ax2.set_xlabel("Angle from the symmetry plane [deg]")
    ax2.set_ylabel("Arrival vs 1D (same flux) [%]")
    ax2.set_title("(b) Arrival time")
    for ax in (ax1, ax2):
        ax.set_xticks(range(0, 91, 15))
        ax.grid(True)
    ax1.legend(fontsize=7.5)
    save(fig, "stage5_g1_axisymmetric")


FIGURES = {
    "isentropes": fig_isentropes,
    "products_temperature": fig_products_temperature,
    "transfer_matching": fig_transfer_matching,
    "phase_transfer": fig_stage1_phase_transfer,
    "al_kinetics": fig_stage1_al_kinetics,
    "al_solver": fig_stage1_al_solver,
    "high_mach": fig_stage2_high_mach,
    "advection": fig_stage2_multiphase_advection,
    "gunib1d": fig_stage2_gunib1d,
    "planar_znd": fig_stage2_planar_znd,
    "spherical_blast": fig_stage3_spherical_blast,
    "blast_kb": fig_stage3_blast_kb,
    "spherical_multiphase": fig_stage3_spherical_multiphase,
    "e2_ignition": fig_stage3_e2_ignition,
    "g1_axisymmetric": fig_stage5_g1_axisymmetric,
}


def main():
    names = sys.argv[1:] or list(FIGURES)
    unknown = [n for n in names if n not in FIGURES]
    if unknown:
        sys.exit(f"Unknown figures {unknown}; choose from {list(FIGURES)}")
    os.makedirs(FIG_DIR, exist_ok=True)
    for n in names:
        FIGURES[n]()


if __name__ == "__main__":
    main()
