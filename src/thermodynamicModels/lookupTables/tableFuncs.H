typedef scalar (*modFuncType)(const scalar&);

inline static scalar noneS(const scalar& x)
{
    return x;
}

inline static scalar lnS(const scalar& x)
{
    return Foam::log(Foam::max(x, small));
}

inline static scalar expS(const scalar& x)
{
    return Foam::exp(x);
}

inline static scalar log10S(const scalar& x)
{
    return Foam::log10(Foam::max(x, small));
}

inline static scalar pow10S(const scalar& x)
{
    return Foam::pow(10, x);
}

inline void setMod
(
    const word& mod,
    modFuncType& modF,
    modFuncType& invModF
)
{
    if (mod == "log10")
    {
        modF = &log10S;
        invModF = &pow10S;
    }
    else if (mod == "pow10")
    {
        modF = &pow10S;
        invModF = &log10S;
    }
    else if (mod == "ln")
    {
        modF = &lnS;
        invModF = &expS;
    }
    else if (mod == "exp")
    {
        modF = &expS;
        invModF = &lnS;
    }
    else
    {
        modF = &noneS;
        invModF = &noneS;
    }
}
