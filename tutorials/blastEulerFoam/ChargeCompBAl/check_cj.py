import numpy as np
from scipy.optimize import fsolve

def solve_cj(rho0, A, B, R1, R2, omega, Cv, E0_vol, D):
    """
    Calculates the CJ state for a JWL explosive.
    """
    e0 = E0_vol / rho0

    def get_p_jwl(V, T):
        # P = A*exp(-R1*V) + B*exp(-R2*V) + omega*rho0*(Cv*T + e0)/V
        return A * np.exp(-R1*V) + B * np.exp(-R2*V) + omega * rho0 * (Cv * T + e0) / V

    def get_cs_sq(V, T):
        rho = rho0 / V
        p = get_p_jwl(V, T)
        pi_val = A*(omega/(R1*V) - 1.0)*np.exp(-R1*V) + B*(omega/(R2*V) - 1.0)*np.exp(-R2*V)
        
        # Consistent with JWLI.H implementation in blastFoam
        delta = (
            A * np.exp(-R1*V) * (omega*(1.0/(R1*rho0) + 1.0/rho) - R1*rho0/(rho**2))
          + B * np.exp(-R2*V) * (omega*(1.0/(R2*rho0) + 1.0/rho) - R2*rho0/(rho**2))
        ) / omega
        
        h = ((omega + 1.0)*p + pi_val)/(omega*rho)
        return (h - delta)*omega

    def equations(vars):
        V, T = vars
        p = get_p_jwl(V, T)
        # Rayleigh line: P = rho0 * D^2 * (1 - V)
        eq1 = p - rho0 * (D**2) * (1.0 - V)
        # Sonic condition: D * (rho0/rho) = c  => D*V = c
        eq2 = (D**2) * (V**2) - get_cs_sq(V, T)
        return [eq1, eq2]

    try:
        # Initial guess: V ~ 0.75, T ~ 3000K
        V_cj, T_cj = fsolve(equations, [0.75, 3000.0])
        p_cj = get_p_jwl(V_cj, T_cj)
        return V_cj, p_cj, T_cj
    except Exception:
        return None, None, None

if __name__ == "__main__":
    # Parameters from constant/physicalProperties
    rho0 = 1820
    A = 617e9
    B = 16.9e9
    R1 = 4.4
    R2 = 1.2
    e0 = 4.45e6;
    omega = 0.28
    Cv = 1383
    E0_vol = 8.1e9
    D = 7612

    print(f"Detonation Properties for:")
    print(f"rho0: {rho0} kg/m3")
    print(f"D:    {D} m/s")
    print(f"E0:   {E0_vol} J/m3")
    print("-" * 30)

    V_cj, p_cj, T_cj = solve_cj(rho0, A, B, R1, R2, omega, Cv, E0_vol, D)
    
    if V_cj:
        print(f"Theoretical CJ State:")
        print(f"V_cj:   {V_cj:.4f}")
        print(f"p_cj:   {p_cj/1e9:.2f} GPa")
        print(f"T_cj:   {T_cj:.2f} K")
        print(f"rho_cj: {rho0/V_cj:.2f} kg/m3")
    else:
        print("Error: Could not solve for CJ state. Check parameter consistency.")
