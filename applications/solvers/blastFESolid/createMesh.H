    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::READ_IF_PRESENT
        )
    );

    const label order = args.optionLookupOrDefault("o", 1);

    List<element> elements(mesh.nCells());
    List<List<element>> boundaryElements(mesh.boundaryMesh().size());
    forAll(mesh.cells(), celli)
    {
        element& elem = elements[celli];
        elem.initialize(mesh, celli, order, GeoType::CELL);
    }

    forAll(boundaryElements, patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (!isA<emptyPolyPatch>(pp))
        {
            boundaryElements[patchi].setSize(pp.size());
        }

        forAll(boundaryElements[patchi], facei)
        {
            element& elem = boundaryElements[patchi][facei];
            elem.initialize(mesh, pp.start() + facei, order, GeoType::FACE);
        }
    }

    scalar minDx = great;
    Field<vector> nodes(mesh.points());
    if (order > 1)
    {
        DynamicList<vector> dnodes(mesh.points());
        List<labelList> edgeNodes(mesh.nEdges());
        List<labelList> faceNodes(mesh.nFaces());
        forAll(elements, celli)
        {
            element& elem = elements[celli];
            DynamicList<label> newElem(elem);

            pointField elemNodes(mesh.points(), elem);
            elemNodes = elem.getNodes(elemNodes);

            label ni = newElem.size();

            // Add edge nodes
            const label nEdgeNodes = elem.nEdgeNodes();
            const edgeList edges(elem.edges());
            forAll(edges, ei)
            {
                const edge& e = edges[ei];
                const label edgei = meshTools::findEdge(mesh, e[0], e[1]);
                if (!edgeNodes[edgei].size())
                {
                    edgeNodes[edgei].setSize(nEdgeNodes);
                    for (label eni = 0; eni < nEdgeNodes; eni++)
                    {
                        newElem.append(dnodes.size());
                        edgeNodes[edgei][eni] = dnodes.size();
                        dnodes.append(elemNodes[ni + eni]);
                    }
                }
                else
                {
                    newElem.append
                    (
                        finiteElement::matchNodes
                        (
                            SubList<vector>(elemNodes, nEdgeNodes, ni),
                            dnodes,
                            edgeNodes[edgei]
                        )
                    );
                }
                ni = newElem.size();
            }

            // Add face nodes
            const faceList faces(elem.faces());
            forAll(faces, fi)
            {
                const face& f = faces[fi];
                const cell& c = mesh.cells()[celli];
                label facei = -1;
                forAll(c, fj)
                {
                    if (f == mesh.faces()[c[fj]])
                    {
                        facei = c[fj];
                        break;
                    }
                }

                const label nFaceNodes = elem.nFaceNodes(fi);

                if (!faceNodes[facei].size())
                {
                    faceNodes[facei].setSize(nFaceNodes);
                    for (label fni = 0; fni < nFaceNodes; fni++)
                    {
                        newElem.append(dnodes.size());
                        faceNodes[facei][fni] = dnodes.size();
                        dnodes.append(elemNodes[ni + fni]);
                    }
                }
                else
                {
                    newElem.append
                    (
                        finiteElement::matchNodes
                        (
                            SubList<vector>(elemNodes, nFaceNodes, ni),
                            dnodes,
                            faceNodes[facei]
                        )
                    );
                }
                ni = newElem.size();
            }


            // Add internal nodes
            for (ni = newElem.size(); ni < elemNodes.size(); ni++)
            {
                newElem.append(dnodes.size());
                dnodes.append(elemNodes[ni]);
            }

            elem.transfer(newElem);
        }
        nodes.transfer(dnodes);

        forAll(boundaryElements, patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            forAll(boundaryElements[patchi], facei)
            {
                element& elem = boundaryElements[patchi][facei];
                DynamicList<label> newElem(elem);

                // All nodes have already been added

                // Add edge nodes
                const edgeList edges(elem.edges());
                forAll(edges, ei)
                {
                    const edge& e = edges[ei];
                    const label edgei = meshTools::findEdge(mesh, e[0], e[1]);
                    if (edge::compare(e, mesh.edges()[edgei]) > 0)
                    {
                        newElem.append(edgeNodes[edgei]);
                    }
                    else
                    {
                        newElem.append(reverseList(edgeNodes[edgei]));
                    }
                }

                // Add face nodes
                newElem.append(faceNodes[pp.start() + facei]);

                // Transfer indices
                elem.transfer(newElem);
            }
        }
    }

//     forAll(elements, celli)
//     {
//         elements[celli].setIntegrationRule
//         (
//             elements[celli].elem().ir(),
//             nodes
//         );
//         const UIndirectList<vector> pts(nodes, elements[celli]);
//         forAll(pts, i)
//             for (label j = i+1; j < pts.size(); j++)
//                 minDx = min(mag(pts[i] - pts[j]), minDx);
//     }
//     forAll(boundaryElements, patchi)
//     {
//         forAll(boundaryElements[patchi], facei)
//         {
//             boundaryElements[patchi][facei].setIntegrationRule
//             (
//                 boundaryElements[patchi][facei].elem().ir(),
//                 nodes
//             );
//         }
//     }
