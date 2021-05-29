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
    else if (mod == "none")
    {
        modF = &noneS;
        invModF = &noneS;
    }
    else
    {
        FatalErrorInFunction
            << mod << " is not a valid mod scheme" << nl
            << "Options are: " << nl
            << "    none" << nl
            << "    log10" << nl
            << "    pow10" << nl
            << "    ln" << nl
            << "    exp" << nl
            << abort(FatalError);
    }
}


//- Interpolation types
typedef scalar (*interpFuncType)(const scalar&, const scalar&, const scalar&);

inline static scalar linearClampInterp
(
    const scalar& x,
    const scalar& xLow,
    const scalar& xHigh
)
{
    return max(min((x - xLow)/(xHigh - xLow), 1.0), 0.0);
}


inline static scalar linearExtrapolateInterp
(
    const scalar& x,
    const scalar& xLow,
    const scalar& xHigh
)
{
    return (x - xLow)/(xHigh - xLow);
}


inline static scalar ceilInterp
(
    const scalar& x,
    const scalar& xLow,
    const scalar& xHigh
)
{
    return 1.0;
}


inline static scalar floorInterp
(
    const scalar& x,
    const scalar& xLow,
    const scalar& xHigh
)
{
    return 0.0;
}


inline void setInterp
(
    const word& interpolationScheme,
    interpFuncType& interp
)
{
    if (interpolationScheme == "linearClamp")
    {
        interp = &linearClampInterp;
    }
    else if (interpolationScheme == "linearExtrapolated")
    {
        interp = &linearExtrapolateInterp;
    }
    else if (interpolationScheme == "ceil")
    {
        interp = &ceilInterp;
    }
    else if (interpolationScheme == "floor")
    {
        interp = &floorInterp;
    }
    else
    {
        FatalErrorInFunction
            << interpolationScheme << " is not a valid interpolation scheme" << nl
            << "Options are: " << nl
            << "    linearClamp" << nl
            << "    linearExtrapolated" << nl
            << "    ceil" << nl
            << "    floor" << nl
            << abort(FatalError);
    }
}
