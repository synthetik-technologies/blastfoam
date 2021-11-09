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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedMeshMapper::mapImmersedToBoundary(const Field<Type>& fDist) const
{
    Field<Type> f;
    if (Pstream::parRun())
    {
        List<Field<Type>> sendData(Pstream::nProcs());
        sendData[Pstream::myProcNo()] = fDist;
        Pstream::gatherList(sendData);
        Pstream::scatterList(sendData);

        //- build global field
        forAll(sendData, i)
        {
            f.append(sendData[i]);
        }
    }
    else
    {
        f = fDist;
    }

    tmp<Field<Type>> tmpF
    (
        new Field<Type>(immersedObjectPtr_->size(), Zero)
    );
    Field<Type>& F(tmpF.ref());

    forAll(immersedMeshToObjectPoints_, i)
    {
        forAll(immersedMeshToObjectPoints_[i], j)
        {
            if
            (
                globalImmersedMeshFaceMapPtr_->isLocal
                (
                    immersedMeshToObjectPoints_[i][j]
                )
            )
            {
                label facei =
                    globalImmersedMeshFaceMapPtr_->toLocal
                    (
                        immersedMeshToObjectPoints_[i][j]
                    );
                F[i] += f[facei]*immersedMeshToObjectWeights_[i][j];
            }
        }
    }
    Pstream::listCombineGather(F, plusEqOp<Type>());
    Pstream::listCombineScatter(F);
    return tmpF;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedMeshMapper::mapObjectToBoundary(const Field<Type>& fDist) const
{
    Field<Type> f(fDist);
    if (Pstream::parRun())
    {
        reduce(f, sumOp<List<Type>>());
    }
    tmp<Field<Type>> tmpF
    (
        new Field<Type>(immersedMesh_.boundaryMesh()[interfaceIndex_].size(), Zero)
    );
    Field<Type>& F(tmpF.ref());

    forAll(immersedObjectToMeshPoints_, i)
    {
        forAll(immersedObjectToMeshPoints_[i], j)
        {
            label facei = immersedObjectToMeshPoints_[i][j];
            F[i] += f[facei]*immersedObjectToMeshWeights_[i][j];
        }
    }
    return tmpF;

}

// ************************************************************************* //
