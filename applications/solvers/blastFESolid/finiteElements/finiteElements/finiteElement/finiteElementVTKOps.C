/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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

#include "finiteElement.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::finiteElement::barycentricToVTKTriangle
(
    label* b,
    label ref
)
{
   // Cf. https://git.io/JvW8f
   label maxI = ref;
   label minI = 0;
   label bmin = min(min(b[0], b[1]), b[2]);
   label idx = 0;

   // scope into the correct triangle
   while (bmin > minI)
   {
      idx += 3*ref;
      maxI -= 2;
      ++minI;
      ref -= 3;
   }
   for (label d = 0; d < 3; d++)
   {
      if (b[(d + 2) % 3] == maxI)
      {
         // we are on a vertex
         return idx;
      }
      ++idx;
   }
   for (label d = 0; d < 3; d++)
   {
      if (b[(d + 1) % 3] == minI)
      {
         // we are on an edge
         return idx + b[d] - (minI + 1);
      }
      idx += maxI - (minI + 1);
   }
   return idx;
}


Foam::label Foam::finiteElement::barycentricToVTKTetra(label* b, label ref)
{
    // Cf. https://git.io/JvW8c
    label idx = 0;

    label maxI = ref;
    label minI = 0;

    label bmin = min(min(min(b[0], b[1]), b[2]), b[3]);

    // scope into the correct tetra
    while (bmin > minI)
    {
        idx += 2*(ref*ref + 1);
        maxI -= 3;
        minI++;
        ref -= 4;
    }

    // When a linearized tetra vertex is cast into barycentric coordinates,
    // one of its coordinates is maximal and the other three are minimal.
    // These are the indices of the maximal barycentric coordinate for each
    // vertex.
    static const label vertexMaxCoords[4] = {3, 0, 1, 2};

    // Each linearized tetra edge holds two barycentric tetra coordinates
    // constant and varies the other two. These are the coordinates that are
    // held constant for each edge.
    static const label edgeMinCoords[6][2] =
    {
        {1, 2},
        {2, 3},
        {0, 2},
        {0, 1},
        {1, 3},
        {0, 3}
    };

    // The coordinate that increments when traversing an edge (i.e. the
    // coordinate of the nonzero component of the second vertex of the edge).
    static const label edgeCountingCoord[6] = {0, 1, 3, 2, 2, 2};

    // When describing a linearized tetra face, there is a mapping between the
    // four-component barycentric tetra system and the three-component
    // barycentric triangle system. These are the constant indices within the
    // four-component system for each face (e.g. face 0 holds barycentric
    // tetra coordinate 1 constant).
    static const label faceMinCoord[4] = {1, 3, 0, 2};

    // When describing a linearized tetra face, there is a mapping between the
    // four-component barycentric tetra system and the three-component
    // barycentric triangle system. These are the relevant indices within the
    // four-component system for each face (e.g. face 0 varies across the
    // barycentric tetra coordinates 0, 2 and 3).
    static const label faceBCoords[4][3] =
    {
        {0, 2, 3},
        {2, 0, 1},
        {2, 1, 3},
        {1, 0, 3}
    };


    for (label vertex = 0; vertex < 4; vertex++)
    {
        if (b[vertexMaxCoords[vertex]] == maxI)
        {
            // we are on a vertex
            return idx;
        }
        idx++;
    }

    for (label edge = 0; edge < 6; edge++)
    {
        if
        (
            b[edgeMinCoords[edge][0]] == minI
        && b[edgeMinCoords[edge][1]] == minI
        )
        {
            // we are on an edge
            return idx + b[edgeCountingCoord[edge]] - (minI + 1);
        }
        idx += maxI - (minI + 1);
    }

    for (label face = 0; face < 4; face++)
    {
        if (b[faceMinCoord[face]] == minI)
        {
            // we are on a face
            label projectedb[3];
            for (label i = 0; i < 3; i++)
            {
                projectedb[i] = b[faceBCoords[face][i]] - minI;
            }
            // we must subtract the indices of the face's vertices and edges, which
            // total to 3*ref
            return (idx + barycentricToVTKTriangle(projectedb, ref) - 3*ref);
        }
        idx += (ref+1)*(ref+2)/2 - 3*ref;
    }
    return idx;
}


Foam::label Foam::finiteElement::vtkTriangleDOFOffset
(
    const label ref,
    const label i,
    const label j
)
{
    return i + ref*(j - 1) - (j*(j + 1))/2;
}


Foam::label Foam::finiteElement::cartesianToVTKPrism
(
    const label i,
    const label j,
    const label k,
    const label ref
)
{
    // Cf. https://git.io/JvW0M
    label om1 = ref - 1;
    label ibdr = (i == 0);
    label jbdr = (j == 0);
    label ijbdr = (i + j == ref);
    label kbdr = (k == 0 || k == ref);

    // How many boundaries do we lie on at once?
    label nbdr = ibdr + jbdr + ijbdr + kbdr;

    // Return an invalid index given invalid coordinates
    if
    (
        i < 0 || i > ref
     || j < 0 || j > ref || i + j > ref
     || k < 0 || k > ref
    )
    {
        FatalErrorInFunction
            << "Invalid index" << endl
            << abort(FatalError);
    }

    // Vertex DOF
    if (nbdr == 3)
    {
        // ijk is a corner node. Return the proper index (somewhere in [0,5]):
        return (ibdr && jbdr ? 0 : (jbdr && ijbdr ? 1 : 2)) + (k ? 3 : 0);
    }

    label offset = 6;

    // Edge DOF
    if (nbdr == 2)
    {
        if (!kbdr)
        {
            // Must be on a vertical edge and 2 of {ibdr, jbdr, ijbdr} are true
            offset += om1*6;
            return
                offset
              + (k - 1)
              + ((ibdr && jbdr) ? 0 : (jbdr && ijbdr ? 1 : 2))*om1;
        }
        else
        {
            // Must be on a horizontal edge and kbdr plus 1 of
            // {ibdr, jbdr, ijbdr} is true
            // Skip past first 3 edges if we are on the top (k = ref) face:
            offset += (k == ref ? 3*om1 : 0);

            if (jbdr)
            {
                return offset + i - 1;
            }

            // Skip the i-axis edge
            offset += om1;
            if (ijbdr)
            {
                return offset + j - 1;
            }

            // Skip the ij-axis edge
            offset += om1;

            // if (ibdr)
            return offset + (ref - j - 1);
        }
    }

    // Skip all the edges
    offset += 9*om1;

    // Number of points on a triangular face (but not on edge/corner):
    label ntfdof = (om1 - 1)*om1/2;
    label nqfdof = om1*om1;

    // Face DOF
    if (nbdr == 1)
    {
        if (kbdr)
        {
            // We are on a triangular face.
            if (k > 0)
            {
                offset += ntfdof;
            }
            return offset + vtkTriangleDOFOffset(ref, i, j);
        }

        // Not a k-normal face, so skip them:
        offset += 2*ntfdof;

        // Face is quadrilateral (ref - 1) x (ref - 1)
        // First face is i-normal, then ij-normal, then j-normal

        // On i-normal face
        if (jbdr)
        {
            return offset + (i - 1) + om1*(k - 1);
        }

        // Skip i-normal face
        offset += nqfdof;

        // on ij-normal face
        if (ijbdr)
        {
            return offset + (ref - i - 1) + om1*(k - 1);
        }

        // Skip ij-normal face
        offset += nqfdof;

        return offset + j - 1 + om1*(k - 1);
    }

    // Skip all face DOF
    offset += 2*ntfdof + 3*nqfdof;

    // nbdr == 0: Body DOF
    return offset + vtkTriangleDOFOffset(ref, i, j) + ntfdof*(k - 1);
    // (i - 1) + (ref-1)*((j - 1) + (ref - 1)*(k - 1)));
}


Foam::label Foam::finiteElement::cartesianToVTKTensor
(
    const label idx_in,
    const label ref,
    const ElementType::Type geom
)
{
    label n = ref + 1;
    switch (geom)
    {
        case ElementType::PT:
            return idx_in;

        case ElementType::SEG:
            if (idx_in == 0 || idx_in == ref)
            {
                return idx_in ? 1 : 0;
            }
            return idx_in + 1;

        case ElementType::QUAD:
        {
            // Cf: https://git.io/JvZLT
            label i = idx_in % n;
            label j = idx_in / n;

            // Do we lie on any of the edges
            bool ibdr = (i == 0 || i == ref);
            bool jbdr = (j == 0 || j == ref);

            // Vertex DOF
            if (ibdr && jbdr)
            {
                return (i ? (j ? 2 : 1) : (j ? 3 : 0));
            }
            label offset = 4;

            // Edge DOF on j==0 or j==ref
            if (jbdr)
            {
                return (i - 1) + (j ? ref - 1 + ref - 1 : 0) + offset;
            }

            // Edge DOF on i==0 or i==ref
            else if (ibdr)
            {
                return
                    (j - 1)
                  + (i ? ref - 1 : 2*(ref - 1) + ref - 1)
                  + offset;
            }

            // Interior DOF
            else
            {
                offset += 2*(ref - 1 + ref - 1);
                return offset + (i - 1) + (ref - 1)*((j - 1));
            }
        }
        case ElementType::HEX:
        {
            // Cf: https://git.io/JvZLe
            label i = idx_in % n;
            label j = (idx_in/n) % n;
            label k = idx_in/(n*n);

            bool ibdr = (i == 0 || i == ref);
            bool jbdr = (j == 0 || j == ref);
            bool kbdr = (k == 0 || k == ref);

            // How many boundaries do we lie on at once?
            label nbdr = (ibdr ? 1 : 0) + (jbdr ? 1 : 0) + (kbdr ? 1 : 0);

            // Vertex DOF
            if (nbdr == 3)
            {
                // ijk is a corner node. Return the proper index (in [0,7])
                return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
            }

            label offset = 8;

            // Edge DOF
            if (nbdr == 2)
            {
                if (!ibdr)
                {
                // On i axis
                return
                    (i - 1)
                  + (j ? ref - 1 + ref - 1 : 0)
                  + (k ? 2*(ref - 1 + ref - 1) : 0)
                  + offset;
                }
                if (!jbdr)
                {
                // On j axis
                return
                    (j - 1)
                  + (i ? ref - 1 : 2*(ref - 1) + ref - 1)
                  + (k ? 2*(ref - 1 + ref - 1) : 0)
                  + offset;
                }

                // !kbdr, On k axis
                offset += 4*(ref - 1) + 4*(ref - 1);
                return
                    (k - 1)
                  + (ref - 1)*(i ? (j ? 3 : 1) : (j ? 2 : 0))
                  + offset;
            }

            offset += 4*(ref - 1 + ref - 1 + ref - 1);

            // Face DOF
            if (nbdr == 1)
            {
                // On i-normal face
                if (ibdr)
                {
                return
                    (j - 1)
                  + ((ref - 1)*(k - 1))
                  + (i ? (ref - 1)*(ref - 1) : 0)
                  + offset;
                }

                offset += 2*(ref - 1)*(ref - 1);

                // On j-normal face
                if (jbdr)
                {
                return
                    (i - 1)
                  + ((ref - 1)*(k - 1))
                  + (j ? (ref - 1)*(ref - 1) : 0)
                  + offset;
                }

                offset += 2*(ref - 1)*(ref - 1);

                // kbdr, On k-normal face
                return
                    (i - 1)
                  + ((ref - 1)*(j - 1))
                  + (k ? (ref - 1)*(ref - 1) : 0)
                  + offset;
            }

            // nbdr == 0: Interior DOF
            offset +=
                2*((ref - 1)*(ref - 1)
              + (ref - 1)*(ref - 1)
              + (ref - 1)*(ref - 1));

            return
                offset
              + (i - 1)
              + (ref - 1)*((j - 1) + (ref - 1)*(k - 1));
        }
        default:
            FatalErrorInFunction
                << "cartesianToVTKOrderingTensor only supports tensor "
                << "geometries." << endl
                << abort(FatalError);
        return -1;
    }
}

Foam::labelList Foam::finiteElement::vtkElementConnectivity
(
    const ElementType::Type geom,
    const label ref
)
{
    if (geom == ElementType::TRI)
    {
        labelList con((ref + 1)*(ref + 2)/2);
        label b[3];
        label idx = 0;
        for (b[1]=0; b[1] <= ref; b[1]++)
        {
            for (b[0]=0; b[0] <= ref - b[1]; b[0]++)
            {
                b[2] = ref - b[0] - b[1];
                con[barycentricToVTKTriangle(b, ref)] = idx++;
            }
        }
        return con;
    }
    else if (geom == ElementType::TET)
    {
        labelList con((ref + 1)*(ref + 2)*(ref + 3)/6);
        label idx = 0;
        label b[4];
        for (b[2] = 0; b[2] <= ref; b[2]++)
        {
            for (b[1] = 0; b[1] <= ref - b[2]; b[1]++)
            {
                for (b[0] = 0; b[0] <= ref - b[1] - b[2]; b[0]++)
                {
                    b[3] = ref-b[0]-b[1]-b[2];
                    con[barycentricToVTKTetra(b, ref)] = idx++;
                }
            }
        }
        return con;
    }
    else if (geom == ElementType::PRISM)
    {
        labelList con((ref + 1)*(ref + 1)*(ref + 2)/6);
        label idx = 0;
        for (label k = 0; k <= ref; k++)
        {
            for (label j = 0; j <= ref; j++)
            {
                for (label i = 0; i <= ref - j; i++)
                {
                    con[cartesianToVTKPrism(i, j, k, ref)] = idx++;
                }
            }
        }
        return con;
    }
    else if (geom == ElementType::PYR)
    {
        FatalErrorInFunction
            << "Lagrange pyramid elements not currently supported in VTK."
            << endl
            << abort(FatalError);
        return labelList();
    }
    else
    {
        labelList con;
        if (geom == ElementType::PT)
        {
            con.setSize(1);
        }
        else if (geom == ElementType::SEG)
        {
            con.setSize(ref + 1);
        }
        else if (geom == ElementType::PT)
        {
            con.setSize((ref + 1)*(ref + 1));
        }
        else if (geom == ElementType::HEX)
        {
            con.setSize((ref + 1)*(ref + 1)*(ref + 1));
        }

        forAll(con, idx)
        {
            con[cartesianToVTKTensor(idx, ref, geom)] = idx;
        }
        return con;
    }
}

// ************************************************************************* //

