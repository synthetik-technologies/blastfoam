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


inline static Type trilinearExtrapolateInterp
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
        bilinearExtrapolateInterp
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
        bilinearExtrapolateInterp
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
    return linearExtrapolateInterp(x, 0, {xs[i], xs[i+1]}, {m, p});
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
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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


inline static Type tricubicClampExtrapolated
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
        bicubicExtrapolateInterp
        (
            y,
            z,
            1,
            1,
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
        bicubicExtrapolateInterp
        (
            y,
            z,
            1,
            1,
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
        bicubicExtrapolateInterp
        (
            y,
            z,
            1,
            1,
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
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
        bicubicExtrapolateInterp
        (
            y,
            z,
            1,
            1,
            {ys[i-1], ys[j], ys[j+1], ys[j+2]},
            {zs[i-1], zs[k], zs[k+1], zs[k+2]},
            tf
        )
    );

    return cubicExtrapolateInterp
    (
        x,
        1,
        {xs[i-1], xs[i], xs[i+1], xs[i+2]},
        {fim, fi, fip, fipp}
    );
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
        interp = &trilinearExtrapolateInterp;
    }
    else if (interpolationScheme == "cubicClamp")
    {
        interp = &tricubicClampInterp;
    }
    else if (interpolationScheme == "cubicExtrapolated")
    {
        interp = &tricubicClampExtrapolated;
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
