inline scalar Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Cv(rho, e, T) + EquationOfState::CpMCv(rho, e, T);
}

inline scalar Hs
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Es(rho, e, T) + EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}

inline scalar Ha
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Ea(rho, e, T) + EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}
