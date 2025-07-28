inline scalar Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar cv(Cv(rho, e, T));
    return cv + EquationOfState::CpMCv(rho, e, T, cv);
}

inline scalar Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Es(rho, e, T) + EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}

inline scalar Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ea(rho, e, T) + EquationOfState::p(rho, e, T)/max(rho, 1e-10);
}

inline scalar Cpt(const scalar T) const
{
    return Cvt(T);
}

inline scalar Hst(const scalar T) const
{
    return Est(T);
}

inline scalar Hat(const scalar T) const
{
    return Eat(T);
}
