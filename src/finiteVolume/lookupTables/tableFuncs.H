typedef Type (*modFuncType)(const Type&);

inline static Type noneS(const Type& x)
{
    return x;
}

virtual void setMod
(
    const word& mod,
    modFuncType& modF,
    modFuncType& invModF
)
{
    if (mod == "none")
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
            << abort(FatalError);
    }
}


//- Interpolation types
typedef scalar (*interpFuncType)(const scalar, const scalar, const scalar);

inline static scalar linearClampInterp
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
)
{
    return max(min((x - xLow)/(xHigh - xLow), 1.0), 0.0);
}


inline static scalar linearExtrapolateInterp
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
)
{
    return (x - xLow)/(xHigh - xLow);
}


inline static scalar ceilInterp
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
)
{
    return 1.0;
}


inline static scalar floorInterp
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
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
