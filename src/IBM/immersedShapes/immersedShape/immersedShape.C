/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "immersedShape.H"
#include "immersedBoundaryObject.H"
#include "wedgePolyPatch.H"
#include "SortableList.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "orientedSurface.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedShape, 0);
    defineRunTimeSelectionTable(immersedShape, twoD);
    defineRunTimeSelectionTable(immersedShape, threeD);
}


Foam::immersedShape::immersedShape
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    pMesh_(pMesh),
    object_(ibo),
    dx_(dict.lookup<scalar>("dx")),
    points0_(),
    faceCentresOld_(),
    geometricD_(Zero),
    ri_(-1),
    ai_(-1),
    angle_(0.0),
    scale_(1.0),
    ei_(-1),
    centreOfMass_(Zero),
    orientation_(tensor::I),
    name_(dict.dictName()),
    write_(dict.lookupOrDefault("write", false)),
    patchPtr_(nullptr),
    bb_(Zero, Zero)
{
    if (pMesh.nGeometricD() == 3)
    {
        geometricD_ = vector::one;
    }
    else if (pMesh.nGeometricD() == 2)
    {
        // Empty index is used for 2-D and axisymmetric
        label i = 0;
        forAll(pMesh.geometricD(), dimi)
        {
            if (pMesh.geometricD()[dimi] == -1)
            {
                ei_ = dimi;
            }
            else if (i == 0)
            {
                xi_ = dimi;
                i++;
            }
            else
            {
                yi_ = dimi;
            }
        }

        boundBox bb(pMesh_.points());
        scalar minE(returnReduce(bb.min()[ei_], minOp<scalar>()));
        scalar maxE(returnReduce(bb.max()[ei_], maxOp<scalar>()));

        //- Remove any offset added in the empty direction
        centreOfMass_[ei_] = 0.5*(minE + maxE);
        bb_.min()[ei_] = minE;
        bb_.max()[ei_] = maxE;

        if (pMesh.nSolutionD() == 3) // Axisymmetric
        {
            const polyBoundaryMesh& bMesh = pMesh.boundaryMesh();
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];
                if (isA<wedgePolyPatch>(patch))
                {
                    const wedgePolyPatch& wedge =
                        refCast<const wedgePolyPatch>(patch);
                    vector axis(wedge.axis());
                    angle_ = acos(wedge.cosAngle());

                    scale_ = angle_/Foam::constant::mathematical::pi;

                    ai_ = 0;
                    forAll(axis, i)
                    {
                        if (mag(axis[i]) > mag(axis[ai_]))
                        {
                            ai_ = i;
                        }
                    }

                    if (ai_ != 0 && ei_ != 0)
                    {
                        ri_ = 0;
                    }
                    else if (ai_ != 1 && ei_ != 1)
                    {
                        ri_ = 1;
                    }
                    else if (ai_ != 2 && ei_ != 2)
                    {
                        ri_ = 2;
                    }
                    xi_ = ai_;
                    yi_ = ri_;
                    break;
                }
            }
        }

    }
    else
    {
        FatalErrorInFunction
            << "IBM is not supported in 1D"
            << abort(FatalError);
    }
}


Foam::immersedShape::~immersedShape()
{}


void Foam::immersedShape::read(const dictionary& dict)
{
    if (needCentreOfMass())
    {
        vector com(dict.lookup<vector>("centreOfMass"));
        if (ei_ != -1)
        {
            this->centreOfMass_[xi_] = com[xi_];
            this->centreOfMass_[yi_] = com[yi_];
        }
        else
        {
            this->centreOfMass_ = com;
        }
    }

//     if (needOrientation())
//     {
//         orientation_ =
//             tensor
//             (
//                 dict.lookupOrDefault<tensor>("orientation", tensor::I)
//             );
//     }
}


void Foam::immersedShape::writeVTK() const
{
    if (write_ && Pstream::master())
    {
        Info<<"writing "<< name_ << ".vtk" << endl;
        vtkWritePolyData::write
        (
            name_ + ".vtk",
            name_,
            false,
            patchPtr_->points(),
            labelList(),
            edgeList(),
            patchPtr_()
        );
    }
}


Foam::tensor Foam::immersedShape::RotM(const scalar theta) const
{
    if (ai_ == 0)
    {
        return tensor
        (
            1.0, 0.0, 0.0,
            0.0, cos(theta), -sin(theta),
            0.0, sin(theta), cos(theta)
        );
    }
    else if (ai_ == 1)
    {
        return tensor
        (
            cos(theta), 0.0, sin(theta),
            0.0, 1.0, 0.0,
            -sin(theta), 0.0, cos(theta)
        );
    }
    else
    {
        return tensor
        (
            cos(theta), -sin(theta), 0.0,
            sin(theta), cos(theta), 0.0,
            0.0, 0.0, 1.0
        );
    }
}

Foam::tmp<Foam::scalarField>
Foam::immersedShape::sortPointsPolar(pointField& points) const
{
    //- Store current points for mapping
    pointField oldPoints(points);

    vector centre(sum(points)/scalar(points.size()));
    if (ai_ != -1)
    {
        if (mag(min(points)[ri_]) < small)
        {
            centre[ri_] = 0.0;
        }
    }

    //- Sort points counter-clockwise
    vectorField diff(points - centre);
    scalarField thetas(diff.size(), 0.0);
    forAll(thetas, i)
    {
        thetas[i] = atan2(diff[i][yi_], diff[i][xi_]);
    }

    SortableList<scalar> sortTheta(thetas);
    const labelList& map = sortTheta.indices();

    forAll(points, i)
    {
        points[i] = oldPoints[map[i]];
    }
    return tmp<scalarField>(new scalarField(sortTheta));
}


void Foam::immersedShape::insidePoints
(
    pointField& points,
    labelList& indices
) const
{
    label pi = 0;
    forAll(points, i)
    {
        if (bb_.contains(points[i]))
        {
            points[pi] = points[i];
            indices[pi++] = i;
        }
    }
    points.resize(pi);
    indices.resize(pi);
}


//- Make list of points axisymmetric and create face indexing
void Foam::immersedShape::extrudeAxi
(
    pointField& points,
    faceList& faces,
    const bool onAxis
) const
{
    pointField origPoints(points);

    // First and last points are not duplicated since they are on the axis
    // This must be handled in the get2DPoints functions
    label nx = points.size();
    if (onAxis)
    {
        nx -= 2;
    }

    //- Positive, negative, and end points on the axis
    origPoints = points;
    if (onAxis)
    {
        points.resize(nx*2 + 2);
    }
    else
    {
        points.resize(nx*2);
    }

    label pi = 0;
    //- Extrapolate to first point on axis
    if (onAxis)
    {
        points[pi++] = origPoints[0];
    }

    //- Add points at front
    tensor R(RotM(angle_));
    label start = onAxis ? 1 : 0;
    for (label i = start; i < origPoints.size() - start; i++)
    {
        points[pi++] = R & origPoints[i];
    }

    //- Add points at front
    R = RotM(-angle_);
    for (label i = start; i < origPoints.size() - start; i++)
    {
        points[pi++] = R & origPoints[i];
    }

    //- Extrapolate to last point on axis
    if (onAxis)
    {
        points.last() = origPoints.last();
    }

    faces.resize(2*nx);

    label fi = 0;
    if (onAxis)
    {
        faces[fi++] = labelledTri(0, 1, nx+1, 0);
        for (label i = 1; i < nx; i++)
        {
            faces[fi++] = labelledTri(i, i+1, nx+i+1, 0);
            faces[fi++] = labelledTri(nx+i+1, nx+i, i, 0);
        }
        label Nx = points.size() - 1;
        faces[fi++] = labelledTri(Nx - 1, nx, Nx, 0);
    }
    else
    {
        label fi = 0;
        for (label i = 0; i < origPoints.size() - 1; i++)
        {
            faces[fi++] = labelledTri(i, i+1, nx+i+1, 0);
            faces[fi++] = labelledTri(nx+i+1, nx+i, i, 0);
        }

        label Nx = points.size();
        faces[fi++] = labelledTri(nx-1, nx-2, Nx-2, 0);
        faces[fi++] = labelledTri(Nx-2, Nx-1, nx-1, 0);
    }

    faceList origFaces(faces);
    label facei = 0;
    forAll(origFaces, i)
    {
        if (origFaces[i].mag(points) > small)
        {
            faces[facei++] = origFaces[i];
        }
    }
    faces.resize(facei);
}

void Foam::immersedShape::extrude2D
(
    pointField& points,
    faceList& faces
) const
{
    pointField origPoints(points);

    //- Positive, negative, and end points on the axis
    points.resize(points.size()*2);

    label pi = 0;

    //- Get offsets
    vector offset(Zero);
    boundBox bb(pMesh_.points());
    scalar minE(returnReduce(bb.min()[ei_], minOp<scalar>()));
    scalar maxE(returnReduce(bb.max()[ei_], maxOp<scalar>()));

    offset[ei_] = minE;
    forAll(origPoints, i)
    {
        points[pi++] = origPoints[i] + offset;
    }

    //- Add points at front
    offset[ei_] = maxE;
    forAll(origPoints, i)
    {
        points[pi++] = origPoints[i] + offset;
    }


    label nx = origPoints.size();
    faces.resize(2*nx + 2);
    label fi = 0;
    for (label i = 0; i < origPoints.size() - 1; i++)
    {
        faces[fi++] = labelledTri(i+1, i, nx+i+1, 0);
        faces[fi++] = labelledTri(nx+i, nx+i+1, i, 0);
    }

    label Nx = points.size();
    faces[fi++] = labelledTri(nx-2, nx-1, Nx-2, 0);
    faces[fi++] = labelledTri(Nx-1, Nx-2, nx-1, 0);

    faces[fi++] = labelledTri(nx-1, 0, Nx-1, 0);
    faces[fi++] = labelledTri(nx, Nx-1, 0, 0);

    faceList origFaces(faces);
    label facei = 0;
    forAll(origFaces, i)
    {
        if (origFaces[i].mag(points) > small)
        {
            faces[facei++] = origFaces[i];
        }
    }
    faces.resize(facei);
}


void Foam::immersedShape::refine3D(triSurface& triMesh)
{
    bool changing = true;
    label iter = 0;

    pointField points(triMesh.points());
    if (adjustPoints(points))
    {
        triMesh.movePoints(points);
    }

    const pointField& tPoints(triMesh.points());
    while (changing)
    {
        changing = false;

        labelList refineFaces(triMesh.size(), -1);
        label fi = 0;
        forAll(triMesh, facei)
        {
            label nEdges = 0;
            const face& f(triMesh[facei]);
            forAll(f, pi)
            {
                const point& p1(tPoints[f[pi]]);
                const point& p2(tPoints[f[(pi+1)%3]]);
                scalar L(mag(p1 - p2));
                if (L > dx_)
                {
                    nEdges++;
                }
            }
            if (nEdges > 1)
            {
                refineFaces[fi++] = facei;
                changing = true;
            }
        }
        refineFaces.resize(fi);
        if (debug)
        {
            Info<< "Iter " << iter++ << ": Refining " << fi << " faces" << endl;
        }

        triMesh = triSurfaceTools::redGreenRefine(triMesh, refineFaces);

        pointField points(triMesh.points());
        if (adjustPoints(points))
        {
            triMesh.movePoints(points);
        }
    }

    orientedSurface::orient
    (
        triMesh,
        vector::zero,
        true
    );
}


void Foam::immersedShape::movePoints()
{
    faceCentresOld_ = patchPtr_->faceCentres();
    patchPtr_->movePoints(object_.transform(points0_));

    point mp(min(patchPtr_->points()));
    point Mp(max(patchPtr_->points()));
    if (ai_ != -1)
    {
        bb_.min()[xi_] = mp[xi_];
        bb_.min()[yi_] = mp[yi_];
        bb_.max()[xi_] = Mp[xi_];
        bb_.max()[yi_] = Mp[yi_];
    }
    else
    {
        bb_.min() = mp;
        bb_.max() = Mp;
    }
}


Foam::tmp<Foam::pointField> Foam::immersedShape::discretizeLine
(
    const vector& p1,
    const vector& p2,
    const bool addEnd
) const
{
    vector l(p2 - p1);
    label nPoints = mag(l)/dx_ + 1;
    tmp<pointField> tmpPoints
    (
        new pointField(addEnd ? nPoints : nPoints - 1)
    );
    pointField& points(tmpPoints.ref());

    vector dl = l/scalar(nPoints - 1);
    forAll(points, i)
    {
        points[i] = p1 + dl*scalar(i);
    }
    return tmpPoints;
}


bool Foam::immersedShape::calcLine
(
    const vector& p1,
    const vector& p2,
    scalar& m,
    scalar& b
) const
{
    scalar x1 = p1[xi_];
    scalar y1 = p1[yi_];
    scalar x2 = p2[xi_];
    scalar y2 = p2[yi_];

    // Handle undefined lines
    if (mag(x2 - x1) < small)
    {
        m = great;
        b = -great;
        return true;
    }

    // Calculate slope
    m = (y2 - y1)/(x2 - x1);

    // Horizontal line
    if (mag(m) < small)
    {
        b = y1;
        return false;
    }

    // Calculate the intercept
    b = y1 - m*x1;
    return false;
}


bool Foam::immersedShape::intersection
(
    const vector& p1, const vector& p2,
    const vector& p3, const vector& p4
) const
{
    const scalar& x1 = p1[xi_];
    const scalar& x2 = p2[xi_];
    const scalar& x3 = p3[xi_];
    const scalar& x4 = p4[xi_];

    scalar m12 = min(x1, x2);
    scalar m34 = min(x3, x4);

    scalar M12 = max(x1, x2);
    scalar M34 = max(x3, x4);

    scalar mxa = max(m12, m34);
    scalar Mxa = min(M12, M34);

    // Not in bounds
    if (M12 < m34)
    {
        return false;
    }

    scalar slope12, b12;
    bool vert12 = calcLine(p1, p2, slope12, b12);

    scalar slope34, b34;
    bool vert34 = calcLine(p3, p4, slope34, b34);

    // Parallel lines
    if (mag(slope12 - slope34) < small)
    {
        return false;
    }

    // Vertical lines
    if (vert12)
    {
        if (m34 < m12 && m12 < M34)
        {
            return true;
        }
        return false;
    }
    if (vert34)
    {
        if (m12 < m34 && m12 < M12)
        {
            return true;
        }
        return false;
    }
    scalar xa = (b34 - b12)/(slope12 - slope34);
    if (xa < mxa || xa > Mxa)
    {
        return false;
    }
    return true;
}


Foam::vector Foam::immersedShape::calcNormal
(
    const vector& p1,
    const vector& p2
) const
{
    scalar slope, b;
    vector n(Zero);

    bool vert = calcLine(p1, p2, slope, b);

    if (vert)
    {
        n[xi_] = 1.0;
        return n;
    }

    if (mag(slope) < small)
    {
        n[yi_] = 1.0;
        return n;
    }

    scalar xm(0.5*(p1[xi_] + p2[xi_]));
    scalar ym(0.5*(p1[yi_] + p2[yi_]));
    scalar bm = ym + slope*xm;

    n[xi_] = xm;
    n[yi_] = (ym - bm)/xm;
    n /= mag(n);

    return n;
}


void Foam::immersedShape::correctNormal
(
    vector& normal,
    const point& start,
    const point& end,
    const point& p1,
    const point& p2
) const
{
    if(!intersection(start, end, p1, p2))
    {
        normal *= -1.0;
    }
}

// ************************************************************************* //
