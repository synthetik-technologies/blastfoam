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


inline static Type linearExtrapolatedInterp
(
    const scalar x,
    const label I,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    const label i = min(max(I, 1), xs.size() - 2);
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
    return linearExtrapolatedInterp(x, i, xs, ys);
}


inline static Type quadraticExtrapolatedInterp
(
    const scalar x,
    const label I,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    const label i = min(max(I, 1), xs.size() - 2);
    const scalar x0 = xs[i-1];
    const scalar x1 = xs[i];
    const scalar x2 = xs[i+1];

    return
        (x - x1)*(x - x2)/(x0 - x1)/(x0 - x2)*ys[i-1]
      + (x - x2)*(x - x0)/(x1 - x2)/(x1 - x0)*ys[i]
      + (x - x0)*(x - x1)/(x2 - x0)/(x2 - x1)*ys[i+1];

}


inline static Type quadraticClampInterp
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
    return quadraticExtrapolatedInterp(x, i, xs, ys);
}


inline static Type cubicExtrapolatedInterp
(
    const scalar x,
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

    return
        (x - x1)*(x - x2)*(x - x3)/(x0 - x1)/(x0 - x2)/(x0 - x3)*ys[i-1]
      + (x - x2)*(x - x3)*(x - x0)/(x1 - x2)/(x1 - x3)/(x1 - x0)*ys[i]
      + (x - x3)*(x - x0)*(x - x1)/(x2 - x3)/(x2 - x0)/(x2 - x1)*ys[i+1]
      + (x - x0)*(x - x1)*(x - x2)/(x3 - x0)/(x3 - x1)/(x3 - x2)*ys[i+2];
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
    return cubicExtrapolatedInterp(x, i, xs, ys);
}


inline static Type ceilInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    return ys[x != xs[i] ? i+1 : i];
}


inline static Type floorInterp
(
    const scalar x,
    const label i,
    const List<scalar>& xs,
    const List<Type>& ys
)
{
    return ys[i];
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    const label I,
    const label J,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const Field<Field<Type>>& zs
)
{
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);
    Type m(linearClampInterp(y, j, ys, zs[i]));
    Type p(linearClampInterp(y, j, ys, zs[i+1]));
    return linearClampInterp(x, 0, {xs[i], xs[i+1]}, {m, p});
}


inline static Type bilinearExtrapolatedInterp
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
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);
    Type m(linearExtrapolatedInterp(x, i, xs, zs[i]));
    Type p(linearExtrapolatedInterp(x, i+1, xs, zs[i+1]));

    return linearExtrapolatedInterp(y, 0, {ys[j], ys[j+1]}, {m, p});
}


inline static Type biquadraticClampInterp
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
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);

    Type fim(quadraticClampInterp(x, i-1, xs, zs[i-1]));
    Type fi(quadraticClampInterp(x, i, xs, zs[i]));
    Type fip(quadraticClampInterp(x, i+1, xs, zs[i+1]));

    return quadraticClampInterp
    (
        y,
        1,
        {ys[j-1], ys[j], ys[j+1]},
        {fim, fi, fip}
    );
}


inline static Type biquadraticExtrapolatedInterp
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
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);

    Type fim(quadraticClampInterp(x, i-1, xs, zs[i-1]));
    Type fi(quadraticClampInterp(x, i, xs, zs[i]));
    Type fip(quadraticClampInterp(x, i+1, xs, zs[i+1]));
    return quadraticExtrapolatedInterp
    (
        y,
        1,
        {ys[j-1], ys[j], ys[j+1]},
        {fim, fi, fip}
    );
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


inline static Type bicubicExtrapolatedInterp
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

    Type fim(cubicExtrapolatedInterp(x, i-1, xs, zs[i-1]));
    Type fi(cubicExtrapolatedInterp(x, i, xs, zs[i]));
    Type fip(cubicExtrapolatedInterp(x, i+1, xs, zs[i+1]));
    Type fipp(cubicExtrapolatedInterp(x, i+2, xs, zs[i+2]));
    return cubicExtrapolatedInterp
    (
        y,
        1,
        {ys[j-1], ys[j], ys[j+1], ys[j+2]},
        {fim, fi, fip, fipp}
    );
}


inline static Type bifloor
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
    return zs[i][j];
}


inline static Type biceil
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
    return zs[x != xs[i] ? i+1 : i][y != ys[j] ? j+1 : j];
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
        interp = &bilinearExtrapolatedInterp;
    }
    else if (interpolationScheme == "quadraticClamp")
    {
        interp = &biquadraticClampInterp;
    }
    else if (interpolationScheme == "quadraticExtrapolated")
    {
        interp = &biquadraticExtrapolatedInterp;
    }
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &bicubicClampInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &bicubicExtrapolatedInterp;
    }
    else if (interpolationScheme == "floor")
    {
        interp = &bifloor;
    }
    else if (interpolationScheme == "ceil")
    {
        interp = &biceil;
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
            << "    floor" << nl
            << "    ceil" << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Interpolation types
typedef Type (*interp3DFuncType)
(
    const scalar,
    const scalar,
    const scalar,
    const label,
    const label,
    const label,
    const List<scalar>&,
    const List<scalar>&,
    const List<scalar>&,
    const Field<Field<Field<Type>>>&
);


inline static Type trilinearClampInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label i,
    const label j,
    const label k,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    Field<Field<Type>> tf(2, Field<Type>(2, 0.0));
    tf[0][0] = fs[i][j][k];
    tf[0][1] = fs[i][j][k+1];
    tf[1][0] = fs[i][j+1][k];
    tf[1][1] = fs[i][j+1][k+1];
    Type m
    (
        bilinearClampInterp
        (
            y,
            z,
            0,
            0,
            {ys[j], ys[j+1]},
            {zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j][k];
    tf[0][1] = fs[i+1][j][k+1];
    tf[1][0] = fs[i+1][j+1][k];
    tf[1][1] = fs[i+1][j+1][k+1];
    Type p
    (
        bilinearClampInterp
        (
            y,
            z,
            0,
            0,
            {ys[j], ys[j+1]},
            {zs[k], zs[k+1]},
            tf
        )
    );

    return linearClampInterp(x, 0, {xs[i], xs[i+1]}, {m, p});
}


inline static Type trilinearExtrapolatedInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label i,
    const label j,
    const label k,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    Field<Field<Type>> tf(2, Field<Type>(2));
    tf[0][0] = fs[i][j][k];
    tf[0][1] = fs[i][j][k+1];
    tf[1][0] = fs[i][j+1][k];
    tf[1][1] = fs[i][j+1][k+1];

    Type m
    (
        bilinearExtrapolatedInterp
        (
            y,
            z,
            0,
            0,
            {ys[j], ys[j+1]},
            {zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j][k];
    tf[0][1] = fs[i+1][j][k+1];
    tf[1][0] = fs[i+1][j+1][k];
    tf[1][1] = fs[i+1][j+1][k+1];

    Type p
    (
        bilinearExtrapolatedInterp
        (
            y,
            z,
            0,
            0,
            {ys[j], ys[j+1]},
            {zs[k], zs[k+1]},
            tf
        )
    );

    return linearExtrapolatedInterp(x, 0, {xs[i], xs[i+1]}, {m, p});
}


inline static Type triquadraticClampInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label I,
    const label J,
    const label K,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);
    const label k = min(max(K, 1), zs.size() - 2);

    Field<Field<Type>> tf(3, Field<Type>(3, 0.0));
    tf[0][0] = fs[i-1][j-1][k-1];
    tf[0][1] = fs[i-1][j-1][k];
    tf[0][2] = fs[i-1][j-1][k+1];
    tf[1][0] = fs[i-1][j][k-1];
    tf[1][1] = fs[i-1][j][k];
    tf[1][2] = fs[i-1][j][k+1];
    tf[2][0] = fs[i-1][j+1][k-1];
    tf[2][1] = fs[i-1][j+1][k];
    tf[2][2] = fs[i-1][j+1][k+1];
    Type fim
    (
        biquadraticClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i][j-1][k-1];
    tf[0][1] = fs[i][j-1][k];
    tf[0][2] = fs[i][j-1][k+1];
    tf[1][0] = fs[i][j][k-1];
    tf[1][1] = fs[i][j][k];
    tf[1][2] = fs[i][j][k+1];
    tf[2][0] = fs[i][j+1][k-1];
    tf[2][1] = fs[i][j+1][k];
    tf[2][2] = fs[i][j+1][k+1];
    Type fi
    (
        biquadraticClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j-1][k-1];
    tf[0][1] = fs[i+1][j-1][k];
    tf[0][2] = fs[i+1][j-1][k+1];
    tf[1][0] = fs[i+1][j][k-1];
    tf[1][1] = fs[i+1][j][k];
    tf[1][2] = fs[i+1][j][k+1];
    tf[2][0] = fs[i+1][j+1][k-1];
    tf[2][1] = fs[i+1][j+1][k];
    tf[2][2] = fs[i+1][j+1][k+1];
    Type fip
    (
        biquadraticClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    return quadraticClampInterp
    (
        x,
        1,
        {xs[i-1], xs[i], xs[i+1]},
        {fim, fi, fip}
    );
}


inline static Type triquadraticExtrapolatedInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label I,
    const label J,
    const label K,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    const label i = min(max(I, 1), xs.size() - 2);
    const label j = min(max(J, 1), ys.size() - 2);
    const label k = min(max(K, 1), zs.size() - 2);

    Field<Field<Type>> tf(3, Field<Type>(3, 0.0));
    tf[0][0] = fs[i-1][j-1][k-1];
    tf[0][1] = fs[i-1][j-1][k];
    tf[0][2] = fs[i-1][j-1][k+1];
    tf[1][0] = fs[i-1][j][k-1];
    tf[1][1] = fs[i-1][j][k];
    tf[1][2] = fs[i-1][j][k+1];
    tf[2][0] = fs[i-1][j+1][k-1];
    tf[2][1] = fs[i-1][j+1][k];
    tf[2][2] = fs[i-1][j+1][k+1];
    Type fim
    (
        biquadraticExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i][j-1][k-1];
    tf[0][1] = fs[i][j-1][k];
    tf[0][2] = fs[i][j-1][k+1];
    tf[1][0] = fs[i][j][k-1];
    tf[1][1] = fs[i][j][k];
    tf[1][2] = fs[i][j][k+1];
    tf[2][0] = fs[i][j+1][k-1];
    tf[2][1] = fs[i][j+1][k];
    tf[2][2] = fs[i][j+1][k+1];
    Type fi
    (
        biquadraticExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j-1][k-1];
    tf[0][1] = fs[i+1][j-1][k];
    tf[0][2] = fs[i+1][j-1][k+1];
    tf[1][0] = fs[i+1][j][k-1];
    tf[1][1] = fs[i+1][j][k];
    tf[1][2] = fs[i+1][j][k+1];
    tf[2][0] = fs[i+1][j+1][k-1];
    tf[2][1] = fs[i+1][j+1][k];
    tf[2][2] = fs[i+1][j+1][k+1];
    Type fip
    (
        biquadraticExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1]},
            {zs[k-1], zs[k], zs[k+1]},
            tf
        )
    );

    return quadraticExtrapolatedInterp
    (
        x,
        1,
        {xs[i-1], xs[i], xs[i+1]},
        {fim, fi, fip}
    );
}


inline static Type tricubicClampInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label I,
    const label J,
    const label K,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    const label i = min(max(I, 1), xs.size() - 3);
    const label j = min(max(J, 1), ys.size() - 3);
    const label k = min(max(K, 1), zs.size() - 3);

    Field<Field<Type>> tf(4, Field<Type>(4, 0.0));
    tf[0][0] = fs[i-1][j-1][k-1];
    tf[0][1] = fs[i-1][j-1][k];
    tf[0][2] = fs[i-1][j-1][k+1];
    tf[0][3] = fs[i-1][j-1][k+2];
    tf[1][0] = fs[i-1][j][k-1];
    tf[1][1] = fs[i-1][j][k];
    tf[1][2] = fs[i-1][j][k+1];
    tf[1][3] = fs[i-1][j][k+2];
    tf[2][0] = fs[i-1][j+1][k-1];
    tf[2][1] = fs[i-1][j+1][k];
    tf[2][2] = fs[i-1][j+1][k+1];
    tf[2][3] = fs[i-1][j+1][k+2];
    tf[3][0] = fs[i-1][j+2][k-1];
    tf[3][1] = fs[i-1][j+2][k];
    tf[3][2] = fs[i-1][j+2][k+1];
    tf[3][3] = fs[i-1][j+2][k+2];
    Type fim
    (
        bicubicClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i][j-1][k-1];
    tf[0][1] = fs[i][j-1][k];
    tf[0][2] = fs[i][j-1][k+1];
    tf[0][3] = fs[i][j-1][k+2];
    tf[1][0] = fs[i][j][k-1];
    tf[1][1] = fs[i][j][k];
    tf[1][2] = fs[i][j][k+1];
    tf[1][3] = fs[i][j][k+2];
    tf[2][0] = fs[i][j+1][k-1];
    tf[2][1] = fs[i][j+1][k];
    tf[2][2] = fs[i][j+1][k+1];
    tf[2][3] = fs[i][j+1][k+2];
    tf[3][0] = fs[i][j+2][k-1];
    tf[3][1] = fs[i][j+2][k];
    tf[3][2] = fs[i][j+2][k+1];
    tf[3][3] = fs[i][j+2][k+2];
    Type fi
    (
        bicubicClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j-1][k-1];
    tf[0][1] = fs[i+1][j-1][k];
    tf[0][2] = fs[i+1][j-1][k+1];
    tf[0][3] = fs[i+1][j-1][k+2];
    tf[1][0] = fs[i+1][j][k-1];
    tf[1][1] = fs[i+1][j][k];
    tf[1][2] = fs[i+1][j][k+1];
    tf[1][3] = fs[i+1][j][k+2];
    tf[2][0] = fs[i+1][j+1][k-1];
    tf[2][1] = fs[i+1][j+1][k];
    tf[2][2] = fs[i+1][j+1][k+1];
    tf[2][3] = fs[i+1][j+1][k+2];
    tf[3][0] = fs[i+1][j+2][k-1];
    tf[3][1] = fs[i+1][j+2][k];
    tf[3][2] = fs[i+1][j+2][k+1];
    tf[3][3] = fs[i+1][j+2][k+2];
    Type fip
    (
        bicubicClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i+2][j-1][k-1];
    tf[0][1] = fs[i+2][j-1][k];
    tf[0][2] = fs[i+2][j-1][k+1];
    tf[0][3] = fs[i+2][j-1][k+2];
    tf[1][0] = fs[i+2][j][k-1];
    tf[1][1] = fs[i+2][j][k];
    tf[1][2] = fs[i+2][j][k+1];
    tf[1][3] = fs[i+2][j][k+2];
    tf[2][0] = fs[i+2][j+1][k-1];
    tf[2][1] = fs[i+2][j+1][k];
    tf[2][2] = fs[i+2][j+1][k+1];
    tf[2][3] = fs[i+2][j+1][k+2];
    tf[3][0] = fs[i+2][j+2][k-1];
    tf[3][1] = fs[i+2][j+2][k];
    tf[3][2] = fs[i+2][j+2][k+1];
    tf[3][3] = fs[i+2][j+2][k+2];
    Type fipp
    (
        bicubicClampInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    return cubicClampInterp
    (
        x,
        1,
        {xs[i-1], xs[i], xs[i+1], xs[i+2]},
        {fim, fi, fip, fipp}
    );
}


inline static Type tricubicClampExtrapolatedInterp
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label I,
    const label J,
    const label K,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    const label i = min(max(I, 1), xs.size() - 3);
    const label j = min(max(J, 1), ys.size() - 3);
    const label k = min(max(K, 1), zs.size() - 3);

    Field<Field<Type>> tf(4, Field<Type>(4, 0.0));
    tf[0][0] = fs[i-1][j-1][k-1];
    tf[0][1] = fs[i-1][j-1][k];
    tf[0][2] = fs[i-1][j-1][k+1];
    tf[0][3] = fs[i-1][j-1][k+2];
    tf[1][0] = fs[i-1][j][k-1];
    tf[1][1] = fs[i-1][j][k];
    tf[1][2] = fs[i-1][j][k+1];
    tf[1][3] = fs[i-1][j][k+2];
    tf[2][0] = fs[i-1][j+1][k-1];
    tf[2][1] = fs[i-1][j+1][k];
    tf[2][2] = fs[i-1][j+1][k+1];
    tf[2][3] = fs[i-1][j+1][k+2];
    tf[3][0] = fs[i-1][j+2][k-1];
    tf[3][1] = fs[i-1][j+2][k];
    tf[3][2] = fs[i-1][j+2][k+1];
    tf[3][3] = fs[i-1][j+2][k+2];
    Type fim
    (
        bicubicExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i][j-1][k-1];
    tf[0][1] = fs[i][j-1][k];
    tf[0][2] = fs[i][j-1][k+1];
    tf[0][3] = fs[i][j-1][k+2];
    tf[1][0] = fs[i][j][k-1];
    tf[1][1] = fs[i][j][k];
    tf[1][2] = fs[i][j][k+1];
    tf[1][3] = fs[i][j][k+2];
    tf[2][0] = fs[i][j+1][k-1];
    tf[2][1] = fs[i][j+1][k];
    tf[2][2] = fs[i][j+1][k+1];
    tf[2][3] = fs[i][j+1][k+2];
    tf[3][0] = fs[i][j+2][k-1];
    tf[3][1] = fs[i][j+2][k];
    tf[3][2] = fs[i][j+2][k+1];
    tf[3][3] = fs[i][j+2][k+2];
    Type fi
    (
        bicubicExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i+1][j-1][k-1];
    tf[0][1] = fs[i+1][j-1][k];
    tf[0][2] = fs[i+1][j-1][k+1];
    tf[0][3] = fs[i+1][j-1][k+2];
    tf[1][0] = fs[i+1][j][k-1];
    tf[1][1] = fs[i+1][j][k];
    tf[1][2] = fs[i+1][j][k+1];
    tf[1][3] = fs[i+1][j][k+2];
    tf[2][0] = fs[i+1][j+1][k-1];
    tf[2][1] = fs[i+1][j+1][k];
    tf[2][2] = fs[i+1][j+1][k+1];
    tf[2][3] = fs[i+1][j+1][k+2];
    tf[3][0] = fs[i+1][j+2][k-1];
    tf[3][1] = fs[i+1][j+2][k];
    tf[3][2] = fs[i+1][j+2][k+1];
    tf[3][3] = fs[i+1][j+2][k+2];
    Type fip
    (
        bicubicExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    tf[0][0] = fs[i+2][j-1][k-1];
    tf[0][1] = fs[i+2][j-1][k];
    tf[0][2] = fs[i+2][j-1][k+1];
    tf[0][3] = fs[i+2][j-1][k+2];
    tf[1][0] = fs[i+2][j][k-1];
    tf[1][1] = fs[i+2][j][k];
    tf[1][2] = fs[i+2][j][k+1];
    tf[1][3] = fs[i+2][j][k+2];
    tf[2][0] = fs[i+2][j+1][k-1];
    tf[2][1] = fs[i+2][j+1][k];
    tf[2][2] = fs[i+2][j+1][k+1];
    tf[2][3] = fs[i+2][j+1][k+2];
    tf[3][0] = fs[i+2][j+2][k-1];
    tf[3][1] = fs[i+2][j+2][k];
    tf[3][2] = fs[i+2][j+2][k+1];
    tf[3][3] = fs[i+2][j+2][k+2];
    Type fipp
    (
        bicubicExtrapolatedInterp
        (
            y,
            z,
            1,
            1,
            {ys[j-1], ys[j], ys[j+1], ys[j+2]},
            {zs[k-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    return cubicExtrapolatedInterp
    (
        x,
        1,
        {xs[i-1], xs[i], xs[i+1], xs[i+2]},
        {fim, fi, fip, fipp}
    );
}


inline static Type trifloor
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label i,
    const label j,
    const label k,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    return fs[i][j][k];
}


inline static Type triceil
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label i,
    const label j,
    const label k,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<scalar>& zs,
    const Field<Field<Field<Type>>>& fs
)
{
    return
        fs
        [x != xs[i] ? i+1 : i]
        [y != ys[j] ? j+1 : j]
        [z != zs[k] ? k+1 : k];
}


inline void setInterp
(
    const word& interpolationScheme,
    interp3DFuncType& interp
)
{
    Info<< "Selecting " << interpolationScheme << " interpolation" << endl;
    if (interpolationScheme == "linearClamp")
    {
        interp = &trilinearClampInterp;
    }
    else if (interpolationScheme == "linearExtrapolated")
    {
        interp = &trilinearExtrapolatedInterp;
    }
    else if (interpolationScheme == "quadraticClamp")
    {
        interp = &triquadraticClampInterp;
    }
    else if (interpolationScheme == "quadraticExtrapolated")
    {
        interp = &triquadraticExtrapolatedInterp;
    }
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &tricubicClampInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &tricubicClampExtrapolatedInterp;
    }
    else if (interpolationScheme == "floor")
    {
        interp = &trifloor;
    }
    else if (interpolationScheme == "ceil")
    {
        interp = &triceil;
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
            << "    floor" << nl
            << "    ceil" << nl
            << abort(FatalError);
    }
}
