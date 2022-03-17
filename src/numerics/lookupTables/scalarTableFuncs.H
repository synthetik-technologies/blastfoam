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

virtual void setMod
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


// inline void setInterp
// (
//     const word& interpolationScheme,
//     interp2DFuncType& interp
// )
// {
//     Info<< "Selecting " << interpolationScheme << " interpolation" << endl;
//     if (interpolationScheme == "linearClamp")
//     {
//         interp = &bilinearClampInterp;
//     }
//     else if (interpolationScheme == "linearExtrapolated")
//     {
//         interp = &bilinearExtrapolatedInterp;
//     }
//     else
//     {
//         FatalErrorInFunction
//             << interpolationScheme << " is not a valid interpolation scheme" << nl
//             << "Options are: " << nl
//             << "    linearClamp" << nl
//             << "    linearExtrapolated" << nl
//             << abort(FatalError);
//     }
// }
