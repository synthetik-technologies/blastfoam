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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Pointer to function to lookup
typedef label (*findIndexFunc)(const scalar, const List<scalar>&);

//- Lookup based on uniform indexing
inline static label findUniformIndexes
(
    const scalar xy,
    const List<scalar>& XY
)
{
    scalar xyMin = XY[0];
    scalar dxy = XY[1] - XY[0];

    scalar ij = (xy - xyMin)/dxy;
    if (ij <= 0)
    {
        return 0;
    }
    return min(floor(ij), XY.size() - 2);
}

//- Lookup based on non uniform indexing
inline static label findNonuniformIndexes
(
    const scalar xy,
    const List<scalar>& XY
)
{
    if (xy < XY[0])
    {
        return 0;
    }

    for (label ij = 0; ij < XY.size() - 1; ij++)
    {
        if (xy < XY[ij] && xy < XY[ij+1])
        {
            return ij;
        }
    }
    return XY.size() - 2;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline static scalar linearWeight
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
)
{
    return (x - xLow)/(xHigh - xLow);
}

//- Interpolation types
typedef void (*interpFuncType)
(
    const scalar,
    const label,
    const scalarList&,
    labelList&,
    scalarList&
);


inline static void linearExtrapolatedInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    label hi = min(xs.size() - 1, I + 1);
    label lo = hi - 1;

    indices.setSize(2);
    weights.setSize(2);

    indices[0] = lo;
    indices[1] = hi;

    weights[1] = (x - xs[lo])/(xs[hi] - xs[lo]);
    weights[0] = 1.0 - weights[1];
}


inline static void linearClampInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    if (x < xs[0])
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = 0;
        weights[0] = 1.0;

        return;
    }
    else if (x > xs.last())
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = xs.size() - 1;
        weights[0] = 1.0;
        return;
    }
    return linearExtrapolatedInterp(x, I, xs, indices, weights);
}


inline static void quadraticExtrapolatedInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    if (I < 1 || I > xs.size() - 2)
    {
        return linearExtrapolatedInterp(x, I, xs, indices, weights);
    }

    label lo = I - 1;
    label mid = I;
    label hi = I + 1;

    const scalar& x0 = xs[lo];
    const scalar& x1 = xs[mid];
    const scalar& x2 = xs[hi];

    indices.resize(3);
    weights.resize(3);

    indices[0] = lo;
    indices[1] = mid;
    indices[2] = hi;

    weights[0] = (x - x1)*(x - x2)/(x0 - x1)/(x0 - x2);
    weights[1] = (x - x2)*(x - x0)/(x1 - x2)/(x1 - x0);
    weights[2] = (x - x0)*(x - x1)/(x2 - x0)/(x2 - x1);
}


inline static void quadraticClampInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    if (x < xs[0])
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = 0;
        weights[0] = 1.0;

        return;
    }
    else if (x > xs.last())
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = xs.size() - 1;
        weights[0] = 1.0;
        return;
    }
    return quadraticExtrapolatedInterp(x, I, xs, indices, weights);
}


inline static void cubicExtrapolatedInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    if (I < 1)
    {
        return linearExtrapolatedInterp(x, I, xs, indices, weights);
    }
    if (I > xs.size() - 3)
    {
        return quadraticExtrapolatedInterp(x, I, xs, indices, weights);
    }

    indices.resize(4);
    weights.resize(4);

    indices[0] = I - 1;
    indices[1] = I;
    indices[2] = I + 1;
    indices[3] = I + 2;

    const scalar& x0 = xs[I-1];
    const scalar& x1 = xs[I];
    const scalar& x2 = xs[I+1];
    const scalar& x3 = xs[I+2];


    weights[0] = (x - x1)*(x - x2)*(x - x3)/(x0 - x1)/(x0 - x2)/(x0 - x3);
    weights[1] = (x - x2)*(x - x3)*(x - x0)/(x1 - x0)/(x1 - x2)/(x1 - x3);
    weights[2] = (x - x3)*(x - x0)*(x - x1)/(x2 - x0)/(x2 - x1)/(x2 - x3);
    weights[3] = (x - x0)*(x - x1)*(x - x2)/(x3 - x0)/(x3 - x1)/(x3 - x2);

//     const scalar dx20 = xs[I+1] - xs[I-1];
//     const scalar dx31 = xs[I+2] - xs[I];
//
//     const scalar x3 = pow3(x);
//     const scalar x2 = sqr(x);
//
//     weights[0] = (-x3 + x2*2.0 - x)/dx20;
//     weights[1] = x3*(2.0 - 1.0/dx31) + x2*(1.0/dx31 - 3.0) + 1.0;
//     weights[2] = x3*(1.0/dx20 - 2.0) + x2*(3.0 - 2.0/dx20) + 1.0/dx20;
//     weights[3] = (x3 - x2)/dx31;
}


inline static void cubicClampInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    if (x < xs[0])
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = 0;
        weights[0] = 1.0;

        return;
    }
    else if (x > xs.last())
    {
        indices.setSize(1);
        weights.setSize(1);

        indices[0] = xs.size() - 1;
        weights[0] = 1.0;
        return;
    }
    return cubicExtrapolatedInterp(x, I, xs, indices, weights);
}


inline static void ceilInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = x != xs[I] ? I + 1 : I;
    weights[0] = 1.0;
}


inline static void floorInterp
(
    const scalar x,
    const label I,
    const scalarList& xs,
    labelList& indices,
    scalarList& weights
)
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = I;
    weights[0] = 1.0;
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
        interp = &linearExtrapolatedInterp;
    }
    else if (interpolationScheme == "quadraticClamp")
    {
        interp = &quadraticClampInterp;
    }
    else if (interpolationScheme == "quadraticExtrapolated")
    {
        interp = &quadraticExtrapolatedInterp;
    }
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &cubicClampInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &cubicExtrapolatedInterp;
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
            << "    quadraticClamp" << nl
            << "    quadraticExtrapolated" << nl
            << "    cubicClamp" << nl
            << "    cubicExtrapolated" << nl
            << "    ceil" << nl
            << "    floor" << nl
            << abort(FatalError);
    }
}

