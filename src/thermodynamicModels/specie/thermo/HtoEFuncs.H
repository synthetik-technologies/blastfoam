inline scalar Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp(rho, e, T) - EquationOfState::CpMCv(rho, e, T);
}

inline scalar Es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hs(rho, e, T) - EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}

inline scalar Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ha(rho, e, T) - EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}
