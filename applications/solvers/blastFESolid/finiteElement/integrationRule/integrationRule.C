/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "integrationRule.H"
#include "GaussianQuadrature.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::FixedList
<
    Foam::List<Foam::integrationRule>,
    Foam::ElementType::SIZE
> Foam::integrationRule::integrationRules;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::integrationRule& Foam::integrationRule::ptRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::PT];
    if (!irs.size())
    {
        irs.append
        (
            integrationRule
            (
                {integrationPoint(0.5, 0.5, 0.5, 1.0)}
            )
        );
    }
    return irs[0];
}


const Foam::integrationRule& Foam::integrationRule::segRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::SEG];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {

        ir.setSize(order+1);
        List<scalar> x, w;
        GaussianQuadrature::calcLegendre(order+1, x, w);
        forAll(ir, i)
        {
            ir[i].set1(x[i], w[i]);
        }
    }
    return ir;
}


const Foam::integrationRule& Foam::integrationRule::triRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::TRI];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        switch (order)
        {
            case 0:
            case 1:
            {
                ir.setSize(1);
                ir.setTriMid(0, 0.5);
                break;
            }
            case 2:
            {
                ir.setSize(3);
                ir.setTri3(0, 1.0/6.0, 1.0/6.0);
                break;
            }
            case 3:
            {
                ir.setSize(4);
                ir.setTriMid(0, -9.0/32.0);
                ir.setTri3(1, 0.2, 25.0/96.0);
                break;
            }
            case 4:
            {
                ir.setSize(6);
                ir.setTri3
                (
                    0,
                    0.091576213509770743460,
                    0.054975871827660933819
                );
                ir.setTri3
                (
                    3,
                    0.44594849091596488632,
                    0.11169079483900573285
                );
                break;
            }
        }
    }
    return ir;
}


const Foam::integrationRule& Foam::integrationRule::quadRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::QUAD];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        const integrationRule& irs = segRule(order);
        ir.setSize((order+1)*(order+1));
        label I = 0;
        forAll(irs, j)
        {
            const integrationPoint& ipj = irs[j];
            forAll(irs, i)
            {
                const integrationPoint& ipi = irs[j];
                ir[I++].set2(ipi.x(), ipj.x(), ipi.w()*ipj.w());
            }
        }
    }
    return ir;
}


const Foam::integrationRule& Foam::integrationRule::tetRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::TET];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        switch (order)
        {
            case 0:
            case 1:
            {
                ir.setSize(1);
                ir.setTetMid(0, 1.0/6.0);
                break;
            }
            case 2:
            {
                ir.setSize(4);
                ir.setTet4b(0, 0.58541019662496845446, 1.0/24.0);
                break;
            }
            case 3:
            {
                ir.setSize(5);
                ir.setTetMid(0, -2.0/15.0);
                ir.setTet4b(1, 0.5, 3.0/40.0);
                break;
            }
            case 4:
            {
                ir.setSize(11);
                ir.setTet4(0, 1.0/14.0, 343.0/45000.0);
                ir.setTetMid(4, -74.0/5625.0);
                ir.setTet6(5, 0.10059642383320079500, 28.0/1125.0);
                break;
            }
        }
    }
    return ir;
}


const Foam::integrationRule& Foam::integrationRule::hexRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::HEX];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        const integrationRule& irs = segRule(order);
        ir.setSize((order+1)*(order+1)*(order+1));

        label I = 0;
        forAll(irs, k)
        {
            const integrationPoint& ipk = irs[k];
            forAll(irs, j)
            {
                const integrationPoint& ipj = irs[j];
                forAll(irs, i)
                {
                    const integrationPoint& ipi = irs[i];
                    ir[I++].set3
                    (
                        ipi.x(),
                        ipj.x(),
                        ipk.x(),
                        ipi.w()*ipj.w()*ipk.w()
                    );
                }
            }
        }
    }
    return ir;
}


const Foam::integrationRule& Foam::integrationRule::prismRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::PRISM];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        const integrationRule& irt = triRule(order);
        const integrationRule& irs = segRule(order);

        ir.setSize(irt.size()*irs.size());

        label I = 0;
        forAll(irs, is)
        {
            const integrationPoint& ips = irs[is];
            forAll(irt, it)
            {
                const integrationPoint& ipt = irt[it];
                integrationPoint& ip = ir[I++];
                ip.x() = ipt.x();
                ip.y() = ipt.y();
                ip.z() = ips.x();
                ip.w() = ipt.w()*ips.w();
            }
        }
    }
    return ir;
}

const Foam::integrationRule& Foam::integrationRule::pyrRule
(
    const label order
)
{
    List<integrationRule>& irs = integrationRules[ElementType::PYR];
    if (irs.size() <= order) irs.setSize(order+1);

    integrationRule& ir = irs[order];
    if (!ir.size())
    {
        ir = hexRule(order);

        forAll(ir, i)
        {
            integrationPoint& ip = ir[i];
            const scalar ipz = ip.z();
            ip.x() = ip.x()*(1.0 - ipz);
            ip.x() = ip.y()*(1.0 - ipz);
            ip.x() = ipz;
            ip.w() = ip.w()/3.0;
        }
    }
    return ir;
}

const Foam::integrationRule& Foam::integrationRule::getRule
(
    const ElementType::Type type,
    const label order
)
{
    switch (type)
    {
        case ElementType::PT: return ptRule(order);
        case ElementType::SEG: return segRule(order);
        case ElementType::TRI: return triRule(order);
        case ElementType::QUAD: return quadRule(order);
        case ElementType::TET: return tetRule(order);
        case ElementType::HEX: return hexRule(order);
        case ElementType::PRISM: return prismRule(order);
        case ElementType::PYR: return pyrRule(order);
        default:
        {
            FatalErrorInFunction
                << "Unknown element type" << endl
                << abort(FatalError);
        }
    }
    return ptRule(order);
}


// ************************************************************************* //
