inline scalar Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar cp(Cp(rho, e, T));
    return cp - EquationOfState::CpMCv(rho, e, T, cp);
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


inline scalar Ets
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hts(rho, e, T) - EquationOfState::p(rho, e, T)/max(rho, 1e-10);
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


inline scalar Eta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hta(rho, e, T) - EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}
