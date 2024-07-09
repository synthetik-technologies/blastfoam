#ifndef element_H
#define element_H

#include "finiteElements.H"
#include "vector.H"
#include "tensor.H"
#include "tensor2D.H"
#include "List.H"
#include "SubList.H"
#include "degenerateMatcher.H"
#include "UautoPtr.H"

namespace Foam
{

class element
:
    public List<label>
{
protected:

    const finiteElement* elem_;
    const integrationRule* ir_;
    GeoType type_;

    List<scalarList> shapes_;
    List<scalarRectangularMatrix> dshapes_;
    List<tensor> Js_;
    List<tensor> invJs_;
    List<scalar> Ws_;
    List<scalarRectangularMatrix> Bs_;

public:

    element()
    :
        List<label>(),
        elem_(nullptr),
        ir_(nullptr),
        type_(GeoType::UNKNOWN),

        shapes_(0),
        dshapes_(0),
        Js_(0),
        invJs_(0),
        Ws_(0),
        Bs_(0)
    {}

    element
    (
        const primitiveMesh& mesh,
        const label index,
        const label order,
        const GeoType type
    )
    :
        List<label>(),
        elem_(nullptr),
        ir_(nullptr),
        type_(GeoType::UNKNOWN),

        shapes_(0),
        dshapes_(0),
        Js_(0),
        invJs_(0),
        Ws_(0),
        Bs_(0)
    {
        initialize(mesh, index, order, type);
    }


    inline void setIntegrationRule
    (
        const integrationRule& ir,
        const pointField& nodes
    );

    void initialize
    (
        const primitiveMesh& mesh,
        const label index,
        const label order,
        const GeoType et
    )
    {
        type_ = et;
        switch (type_)
        {
            case GeoType::CELL:
            {
                cellShape cs(degenerateMatcher::match(mesh, index));
                elem_ = getRefFiniteElement(cs.model().name(), order);
                if
                (
                    elem_->elemType() == ElementType::HEX
                )
                {
                    faceList faces(cs.model().faces(cs));
                    const cell c(identity(faces.size()));

                    label bottomFaceI = -1;
                    scalar minz = great;
                    bool flip = true;
                    forAll(faces, fi)
                    {
                        vector fc(faces[fi].centre(mesh.points()));

                        if (fc.z() < minz)
                        {
                            minz = fc.z();
                            bottomFaceI = fi;
                        }
                    }
                    face& bottomFace = faces[bottomFaceI];

                    {
                        vector fa(bottomFace.area(mesh.points()));
                        if (fa.z() < 0) bottomFace.flip();

                        scalarField magPts(mag(pointField(mesh.points(), bottomFace)));
                        inplaceRotateList
                        (
                            bottomFace,
                            bottomFace.size() - findMin(magPts)
                        );
                    }
                    labelList::operator=(bottomFace);

                    edgeList edges = c.edges(faces);
                    label topStart = -1;
                    forAll(edges, ei)
                    {
                        const edge& e = edges[ei];
                        if (e[0] == bottomFace[0])
                        {
                            if (findIndex(bottomFace, e[1]) < 0)
                            {
                                topStart = e[1];
                                break;
                            }
                        }
                        else if (e[1] == bottomFace[0])
                        {
                            if (findIndex(bottomFace, e[0]) < 0)
                            {
                                topStart = e[0];
                                break;
                            }
                        }
                    }

                    face topFace(bottomFace.size());
                    {
                        const label topFaceI =
                            c.opposingFaceLabel(bottomFaceI, faces);
                        label pi = findIndex(faces[topFaceI], topStart);
                        forAll(topFace, i)
                        {
                            topFace[i] = faces[topFaceI][pi];
                            pi = faces[topFaceI].fcIndex(pi);
                        }

                        vector fa(topFace.area(mesh.points()));
                        if (fa.z() < 0) topFace.flip();
                    }
                    labelList::append(topFace);
                }
                else
                {
                    this->transfer(cs);
                }

                break;
            }
            case GeoType::FACE:
            {
                const face& f = mesh.faces()[index];
                if (f.size() > 4 || f.size() < 3)
                {
                    FatalErrorInFunction
                        << "Only triangle or quads are supported for 2D finite elements" << endl
                        << abort(FatalError);
                }
                List<label>::operator=(f);
                elem_ = getRefFiniteElement
                (
                    f.size() == 3 ? "tri" : "quad",
                    order
                );
                break;
            }
            case GeoType::EDGE:
            {
                const edge& e = mesh.edges()[index];
                List<label>::operator=({e[0], e[1]});
                elem_ = getRefFiniteElement("seg", order);
                break;
            }
            case GeoType::POINT:
            {
                this->setSize(1);
                List<label>::operator=(index);
                elem_ = getRefFiniteElement("pt", order);
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Geometric type has not been set" << endl
                    << abort(FatalError);
            }
        }
    }


    //- Return the element type
    inline label type() const {return elem_->type();}
    inline label geoType() const {return toUnderlying(type_);}

    inline label order() const {return elem_->order();}

    inline const integrationRule& ir() const
    {
        return bool(ir_) ? *ir_ : elem_->ir();
    }

//     inline label size() const {return elem_->size();}

    label nDims() const {return elem_->nDims();}

    inline label nPoints() const {return elem_->nPoints();}
    inline label nEdges() const {return elem_->nEdges();}
    inline label nFaces() const {return elem_->nFaces();}
    inline label nNodes() const {return elem_->nNodes();}
    inline label nEdgeNodes() const {return elem_->nEdgeNodes();}
    inline label nFaceNodes(const label fi) const {return elem_->nFaceNodes(fi);}
    inline label nCellNodes() const {return elem_->nCellNodes();}

    inline const UList<vector>& nodes() const {return elem_->nodes();}

    const finiteElement& elem() const
    {
        return *elem_;
    }

    // Return edges
    edgeList edges() const
    {
        const labelList& pts = *this;
        switch (elem_->elemType())
        {
            case ElementType::PT:
                return edgeList();
            case ElementType::SEG:
                return edgeList({edge((*this)[0], (*this)[1])});
            case ElementType::TRI:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[1], pts[2]),
                        edge(pts[2], pts[0])
                    }
                );
            case ElementType::QUAD:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[1], pts[2]),
                        edge(pts[2], pts[3]),
                        edge(pts[3], pts[0])
                    }
                );
            case ElementType::TET:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[0], pts[2]),
                        edge(pts[0], pts[3]),
                        edge(pts[1], pts[2]),
                        edge(pts[1], pts[3]),
                        edge(pts[2], pts[3])
                    }
                );
            case ElementType::HEX:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[1], pts[2]),
                        edge(pts[3], pts[2]),
                        edge(pts[0], pts[3]),
                        edge(pts[4], pts[5]),
                        edge(pts[5], pts[6]),
                        edge(pts[7], pts[6]),
                        edge(pts[4], pts[7]),
                        edge(pts[0], pts[4]),
                        edge(pts[1], pts[5]),
                        edge(pts[2], pts[6]),
                        edge(pts[3], pts[7])
                    }
                );
            case ElementType::WEDGE:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[1], pts[2]),
                        edge(pts[2], pts[0]),
                        edge(pts[3], pts[4]),
                        edge(pts[4], pts[5]),
                        edge(pts[5], pts[3]),
                        edge(pts[0], pts[3]),
                        edge(pts[1], pts[4]),
                        edge(pts[2], pts[5])
                    }
                );
            case ElementType::PYRAMID:
                return edgeList
                (
                    {
                        edge(pts[0], pts[1]),
                        edge(pts[1], pts[2]),
                        edge(pts[2], pts[3]),
                        edge(pts[3], pts[0]),
                        edge(pts[0], pts[4]),
                        edge(pts[1], pts[4]),
                        edge(pts[2], pts[4]),
                        edge(pts[3], pts[4])
                    }
                );
            default:
                NotImplemented;
                return edgeList();
        }
    }
    // Return faces
    faceList faces() const
    {
        const labelList& pts = *this;
        switch (elem_->elemType())
        {
            case ElementType::PT:
            case ElementType::SEG:
                return faceList();
            case ElementType::TRI:
            case ElementType::QUAD:
                return faceList({face(pts)});
            case ElementType::TET:
                return faceList
                (
                    {
                        triFace({pts[1], pts[2], pts[3]}),
                        triFace({pts[0], pts[3], pts[2]}),
                        triFace({pts[0], pts[1], pts[3]}),
                        triFace({pts[0], pts[2], pts[1]})
                    }
                );
            case ElementType::HEX:
                return faceList
                (
                    {
                        face({pts[3], pts[2], pts[1], pts[0]}),
                        face({pts[0], pts[1], pts[5], pts[4]}),
                        face({pts[1], pts[2], pts[6], pts[5]}),
                        face({pts[2], pts[3], pts[7], pts[6]}),
                        face({pts[3], pts[0], pts[4], pts[7]}),
                        face({pts[4], pts[5], pts[6], pts[7]})
                    }
                );
            case ElementType::WEDGE:
                return faceList
                (
                    {
                        triFace({pts[3], pts[2], pts[1]}),
                        triFace({pts[3], pts[4], pts[5]}),
                        face({pts[0], pts[1], pts[4], pts[3]}),
                        face({pts[1], pts[2], pts[5], pts[4]}),
                        face({pts[2], pts[0], pts[3], pts[5]})
                    }
                );
            case ElementType::PYRAMID:
                return faceList
                (
                    {
                        face({pts[3], pts[2], pts[1], pts[0]}),
                        triFace({pts[0], pts[1], pts[4]}),
                        triFace({pts[1], pts[2], pts[4]}),
                        triFace({pts[2], pts[3], pts[4]}),
                        triFace({pts[3], pts[0], pts[4]})
                    }
                );
            default:
                NotImplemented;
                return faceList();
        }
    }

    void setPoints(List<vector>& pts) const
    {
        elem_->setPoints(pts);
    }

    const scalarList& shape(const label i) const
    {
        return shapes_[i];
    }
    void calcShape(const vector& pt, List<scalar>& shape) const
    {
        elem_->calcShape(pt, shape);
    }

    const scalarRectangularMatrix& dshape(const label i) const
    {
        return dshapes_[i];
    }
    void calcDShape
    (
        const vector& pt,
        List<scalar>& shape,
        scalarRectangularMatrix& dshape
    ) const
    {
        elem_->calcDShape(pt, shape, dshape);
    }

    const tensor& J(const label i) const {return Js_[i];}
    tensor J(const vector& ip, const UList<vector>& points) const
    {
        UIndirectList<vector> pts(points, *this);

        List<scalar> shape(this->size());
        scalarRectangularMatrix dshape(this->size(), 3);
        elem_->calcDShape(ip, shape, dshape);

        tensor J(Zero);
        forAll(pts, nodei)
        {
            for (label dimi = 0; dimi < 3; dimi++)
            {
                for (label dimj = 0; dimj < dshape.n(); dimj++)
                {
                    J(dimi, dimj) += pts[nodei][dimi]*dshape(nodei, dimj);
                }
            }
        }
        return J;
    }

    const tensor& Jinv(const label i) const {return invJs_[i];}
    tensor Jinv(const vector& ip, const UList<vector>& points) const
    {
        if (type_ == GeoType::POINT)
        {
            return tensor::zero;
        }

        UIndirectList<vector> pts(points, *this);

        List<scalar> shape(this->size());
        scalarRectangularMatrix dshape(this->size(), 3);
        elem_->calcDShape(ip, shape, dshape);

        tensor J(Zero);
        forAll(pts, nodei)
        {
            for (label dimi = 0; dimi < 3; dimi++)
            {
                for (label dimj = 0; dimj < dshape.n(); dimj++)
                {
                    J(dimi, dimj) += pts[nodei][dimi]*dshape(nodei, dimj);
                }
            }
        }

        if (type_ == GeoType::CELL)
        {
            return J.inv();
        }

        label valid[3] = {-1, -1, -1};
        label vi = 0;
        for (label dimj = 0; dimj < 3; dimj++)
        {
            bool allZero = true;
            for (label dimi = 0; dimi < 3; dimi++)
            {
                if (mag(J(dimi, dimj)) > small)
                {
                    allZero = false;
                    break;
                }
            }
            if (!allZero)
            {
                valid[vi++] = dimj;
            }
        }

        if (vi == 3)
        {
            return J.inv();
        }
        else if (vi == 2)
        {
            const label i = valid[0];
            const label j = valid[1];
            scalar E = sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
            scalar G = sqr(J(0, j)) + sqr(J(1, j)) + sqr(J(2, j));
            scalar F =
                J(0, i)*J(0, j)
              + J(1, i)*J(1, j)
              + J(2, i)*J(2, j);
            const scalar T = 1.0/E*G - F*F;
            E *= T;
            G *= T;
            F *= T;

            tensor invJ(Zero);
            invJ(0, i) = J(0, i)*G - J(0, j)*F;
            invJ(1, i) = J(0, j)*E - J(0, i)*F;
            invJ(2, i) = J(1, i)*G - J(1, j)*F;
            invJ(0, j) = J(1, j)*E - J(1, i)*F;
            invJ(1, j) = J(2, i)*G - J(2, j)*F;
            invJ(2, j) = J(2, j)*E - J(2, i)*F;
            return invJ;

        }
        else if (vi == 1)
        {
            const label i = valid[0];
            const scalar T = 1.0/sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
            tensor invJ(Zero);
            invJ(0, i) = J(0, i)*T;
            invJ(1, i) = J(1, i)*T;
            invJ(2, i) = J(2, i)*T;
            return invJ;
        }
        return tensor::zero;

    }

    const scalar& W(const label i) const {return Ws_[i];}
    scalar W(const vector& ip, const UList<vector>& points) const
    {
        if (type_ == GeoType::POINT)
        {
            return 0.0;
        }

        UIndirectList<vector> pts(points, *this);

        List<scalar> shape(this->size());
        scalarRectangularMatrix dshape(this->size(), 3);
        elem_->calcDShape(ip, shape, dshape);

        tensor J(Zero);
        forAll(pts, nodei)
        {
            for (label dimi = 0; dimi < 3; dimi++)
            {
                for (label dimj = 0; dimj < dshape.n(); dimj++)
                {
                    J(dimi, dimj) += pts[nodei][dimi]*dshape(nodei, dimj);
                }
            }
        }

        if (type_ == GeoType::CELL)
        {
            return det(J);
        }

        label valid[3] = {-1, -1, -1};
        label vi = 0;
        for (label dimj = 0; dimj < 3; dimj++)
        {
            bool allZero = true;
            for (label dimi = 0; dimi < 3; dimi++)
            {
                if (mag(J(dimi, dimj)) > small)
                {
                    allZero = false;
                    break;
                }
            }
            if (!allZero)
            {
                valid[vi++] = dimj;
            }
        }

        if (vi == 3)
        {
            return det(J);
        }
        else if (vi == 2)
        {
            const label i = valid[0];
            const label j = valid[1];
            const scalar E = sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
            const scalar G = sqr(J(0, j)) + sqr(J(1, j)) + sqr(J(2, j));
            const scalar F =
                J(0, i)*J(0, j)
              + J(1, i)*J(1, j)
              + J(2, i)*J(2, j);
            return sqrt(E*G - F*F);
        }
        else if (vi == 1)
        {
            const label i = valid[0];
            return sqrt(sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i)));
        }
        return 0.0;

    }

    //- Compute the gradient at an integration point
    template<class T>
    inline typename outerProduct<T, vector>::type grad
    (
        const UList<T>& vals,
        const scalarRectangularMatrix& dshape
    ) const
    {
        const label nCmpts = pTraits<T>::nComponents;
        typename outerProduct<T, vector>::type g(Zero);
        vector gshape(Zero);
        for (label j = 0; j < dshape.m(); j++)
        {
            for (label dimi = 0; dimi < dshape.n(); dimi++)
            {
                gshape[dimi] = dshape(j, dimi);
            }
            g += vals[(*this)[j]]*gshape;
        }
        return g;
    }

    const scalarRectangularMatrix& B(const label i) const
    {
        return Bs_[i];
    }

//     //- Compute the divergence at an integration point
//     template<class T>
//     inline typename innerProduct<T, vector>::type div
//     (
//         const UList<T>& vals,
//         const scalarRectangularMatrix& dshape
//     ) const
//     {
//         const label nCmpts = pTraits<T>::nComponents;
//         typename innerProduct<T, vector>::type d(Zero);
//         vector gshape(Zero);
//         for (label j = 0; j < dshape.m(); j++)
//         {
//             for (label dimi = 0; dimi < dshape.n(); dimi++)
//             {
//                 gshape[dimi] = dshape(j, dimi);
//             }
//             d += vals[(*this)[j]] & (gshape & invJ_);
//         }
//         return d;
//     }

    using List<label>::operator=;
};


inline void element::setIntegrationRule
(
    const integrationRule& ir,
    const pointField& nodes
)
{
    ir_ = &ir;

    shapes_.setSize(ir.size());
    dshapes_.setSize(ir.size());
    Js_.setSize(ir.size());
    invJs_.setSize(ir.size());
    Ws_.setSize(ir.size());
    Bs_.setSize(ir.size());

    forAll(ir, i)
    {
        const integrationPoint& ip = ir[i];
        elem().calcDShape(ip, shapes_[i], dshapes_[i]);
        Js_[i] = this->J(ip, nodes);
        invJs_[i] = this->Jinv(ip, nodes);
        Ws_[i] = this->W(ip, nodes);

        const tensor& invJ = invJs_[i];
        const scalarList& shape = shapes_[i];
        const scalarRectangularMatrix& dshape = dshapes_[i];

        scalarRectangularMatrix& B = Bs_[i];
        B.setSize(shape.size(), 3);
        B = Zero;

        forAll(shape, si)
        {
            for (label dimi = 0; dimi < dshape.n(); dimi++)
            {
                for (label dimj = 0; dimj < 3; dimj++)
                {
                    B(si, dimj) += dshape(si, dimi)*invJ(dimi, dimj);
                }
            }
        }
    }
}


} // End namespace blast

#endif // element_H
