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


//- Interpolation types
typedef Type (*interp1DFuncType)
(
    const scalar,
    const label,
    const List<scalar>&,
    const List<Type>&
);

inline static scalar linearWeight
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh
)
{
    return (x - xLow)/(xHigh - xLow);
}


inline static Type linearExtrapolateInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    return ys[i] + linearWeight(x, xs[i], xs[i+1])*(ys[i+1] - ys[i]);
}


inline static Type linearClampInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    if (x < xs[0])
    {
        return ys[0];
    }
    else if (x > xs.last())
    {
        return ys.last();
    }
    return linearExtrapolateInterp(x, i, xs, ys);
}


inline static Type cubicExtrapolateInterp
(
    const scalar X,
    const label I,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    const label i = min(max(I, 1), xs.size() - 3);
    const scalar x0 = xs[i-1];
    const scalar x1 = xs[i];
    const scalar x2 = xs[i+1];
    const scalar x3 = xs[i+2];
    const scalar dx10 = x1 - x0;
    const scalar dx21 = x2 - x1;
    const scalar dx32 = x3 - x2;

//     if (mag(dx10 - dx21) < small && mag(dx21 - dx32) < small)
//     {
//         scalar x = linearWeight(X, x1, x2);
//         Info<<"x: "<<x<<endl;
//         return
//             ys[i-1]*(-0.5*pow3(x) + sqr(x) - 0.5*x)
//           + ys[i]*(1.5*pow3(x) - 2.5*sqr(x) + 1.0)
//           + ys[i+1]*(-1.5*pow3(x) + 2.0*sqr(x) + 0.5*x)
//           + ys[i+2]*(0.5*pow3(x) + 0.5*sqr(x));
//     }

    const scalar x = X;
    const scalar dx20 = dx10 + dx21;
    const scalar dx31 = dx21 + dx32;

    const scalar den0 =
        dx20*(sqr(x1) - 2.0*x1*x2 + sqr(x2));
    const scalar den1 =
        dx31*(pow3(x1) - 3.0*sqr(x1)*x2 + 3.0*x1*sqr(x2) - pow3(x2));
    const scalar den2 =
        dx20*(pow3(x1) - 3.0*sqr(x1)*x2 + 3.0*x1*sqr(x2) - pow3(x2));
    const scalar den3 =
        dx31*(sqr(x1) - 2.0*x1*x2 + sqr(x2));


    const Type a =
      - ys[i-1]/den0
      + ys[i]*(x2 - x1 - 2.0*dx31)/den1
      + ys[i+1]*(2.0*dx20 + x1 - x2)/den2
      + ys[i+2]/den3;

    const Type b =
        ys[i-1]*(x1 + 2.0*x2)/den0
      + ys[i]*(3.0*dx31*(x1 + x2) + 2.0*sqr(x1) - x1*x2 - sqr(x2))/den1
      - ys[i+1]*(3.0*dx20*(x1 + x2) + sqr(x1) + x1*x2 - 2.0*sqr(x2))/den2
      - ys[i+2]*(2.0*x1 + x2)/den3;

    const Type c =
      - ys[i-1]*x2*(2.0*x1 + x2)/den0
      - ys[i]*x1*(6.0*dx31*x2 + sqr(x1) + x1*x2 - 2.0*sqr(x2))/den1
      + ys[i+1]*x2*(6.0*dx20*x1 + 2.0*sqr(x1) - x1*x2 - sqr(x2))/den2
      + ys[i+2]*x1*(x1 + 2.0*x2)/den3;

    const Type d =
        ys[i-1]*x1*sqr(x2)/den0
      + ys[i]*x2*(dx31*(3.0*x1*x2 - sqr(x2)) + sqr(x1)*(x1 - x2))/den1
      + ys[i+1]*x1*(dx20*(sqr(x1) - 3.0*x1*x2) + sqr(x2)*(x2 - x1))/den2
      - ys[i+2]*sqr(x1)*x2/den3;

    return a*pow3(x) + b*sqr(x) + c*x + d;

}


inline static Type cubicClampInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    if (x < xs[0])
    {
        return ys[0];
    }
    else if (x > xs.last())
    {
        return ys.last();
    }
    return cubicExtrapolateInterp(x, i, xs, ys);
}


inline static Type ceilInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    return ys[i+1];
}


inline static Type floorInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    return ys[0];
}


inline void setInterp
(
    const word& interpolationScheme,
    interp1DFuncType& interp
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
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &cubicExtrapolateInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &cubicExtrapolateInterp;
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
            << "    cubicClamp" << nl
            << "    cubicExtrapolated" << nl
            << "    ceil" << nl
            << "    floor" << nl
            << abort(FatalError);
    }
}


//- Interpolation types
typedef Type (*interp2DFuncType)
(
    const scalar,
    const scalar,
    const label,
    const label,
    const List<scalar>&,
    const List<scalar>&,
    const Field<Field<Type>>&
);


inline static Type bilinearClampInterp
(
    const scalar x,
    const scalar y,
    const label i,
    const label j,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const Field<Field<Type>>& zs
)
{
    Type m(linearClampInterp(y, j, ys, zs[i]));
    Type p(linearClampInterp(y, j, ys, zs[i+1]));
    return linearClampInterp(x, 0, {xs[i], xs[i+1]}, {m, p});
}


inline static Type bilinearExtrapolateInterp
(
    const scalar x,
    const scalar y,
    const label i,
    const label j,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const Field<Field<Type>>& zs
)
{
    Type m(linearExtrapolateInterp(x, i, xs, zs[i]));
    Type p(linearExtrapolateInterp(x, i+1, xs, zs[i+1]));
    return linearExtrapolateInterp(y, 0, {ys[j], ys[j+1]}, {m, p});
}


inline static Type bicubicClampInterp
(
    const scalar x,
    const scalar y,
    const label I,
    const label J,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const Field<Field<Type>>& zs
)
{
    const label i = min(max(I, 1), xs.size() - 3);
    const label j = min(max(J, 1), ys.size() - 3);

    Type fim(cubicClampInterp(x, i-1, xs, zs[i-1]));
    Type fi(cubicClampInterp(x, i, xs, zs[i]));
    Type fip(cubicClampInterp(x, i+1, xs, zs[i+1]));
    Type fipp(cubicClampInterp(x, i+2, xs, zs[i+2]));
    return cubicClampInterp
    (
        y,
        1,
        {ys[j-1], ys[j], ys[j+1], ys[j+2]},
        {fim, fi, fip, fipp}
    );
}


inline static Type bicubicExtrapolateInterp
(
    const scalar x,
    const scalar y,
    const label I,
    const label J,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const Field<Field<Type>>& zs
)
{
    const label i = min(max(I, 1), xs.size() - 3);
    const label j = min(max(J, 1), ys.size() - 3);

    Type fim(cubicExtrapolateInterp(x, i-1, xs, zs[i-1]));
    Type fi(cubicExtrapolateInterp(x, i, xs, zs[i]));
    Type fip(cubicExtrapolateInterp(x, i+1, xs, zs[i+1]));
    Type fipp(cubicExtrapolateInterp(x, i+2, xs, zs[i+2]));
    return cubicExtrapolateInterp
    (
        y,
        1,
        {ys[j-1], ys[j], ys[j+1], ys[j+2]},
        {fim, fi, fip, fipp}
    );
}


inline void setInterp
(
    const word& interpolationScheme,
    interp2DFuncType& interp
)
{
    Info<< "Selecting " << interpolationScheme << " interpolation" << endl;
    if (interpolationScheme == "linearClamp")
    {
        interp = &bilinearClampInterp;
    }
    else if (interpolationScheme == "linearExtrapolated")
    {
        interp = &bilinearExtrapolateInterp;
    }
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &bicubicClampInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &bicubicExtrapolateInterp;
    }
    else
    {
        FatalErrorInFunction
            << interpolationScheme << " is not a valid interpolation scheme" << nl
            << "Options are: " << nl
            << "    linearClamp" << nl
            << "    linearExtrapolated" << nl
            << "    cubicClamp" << nl
            << "    cubicExtrapolated" << nl
            << abort(FatalError);
    }
}

