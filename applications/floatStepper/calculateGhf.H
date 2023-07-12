/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
advector.surf().reconstruct();
const DynamicField<label>& interfaceLabels = advector.surf().interfaceLabels();
boolList isSurfCell(mesh.nCells(), false);
forAll(interfaceLabels, li)
{
    isSurfCell[interfaceLabels[li]] = true;
}

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();
ghf *= 0.0;
scalarField& ghfIn = ghf.primitiveFieldRef();
const scalarField& alpha1In = alpha1.primitiveField();

// Setting internal ghf values
forAll(ghfIn, facei)
{
    const label ownerCell(own[facei]);
    const label neiCell(nei[facei]);
    if (isSurfCell[ownerCell])
    {
        const point Cio = advector.surf().centre()[ownerCell];
        ghfIn[facei] = (g.value() & Cio) - ghRef.value();
        if (isSurfCell[neiCell])
        {
            const point Cin = advector.surf().centre()[neiCell];
            ghfIn[facei] += (g.value() & Cin) - ghRef.value();
            ghfIn[facei] *= 0.5;
        }
    }
    else if (isSurfCell[neiCell])
    {
        const point Cin = advector.surf().centre()[neiCell];
        ghfIn[facei] = (g.value() & Cin) -  ghRef.value();
    }
    else if (min(alpha1In[ownerCell], alpha1In[neiCell]) < .5 && max(alpha1In[ownerCell], alpha1In[neiCell]) > .5)
    {
        ghfIn[facei] = (g.value() & mesh.Cf()[facei]) - ghRef.value();
    }

}

// Setting boundary ghf values
const polyBoundaryMesh& patches = mesh.boundaryMesh();

forAll(patches, patchi)
{
    // Pout << "patches[patchi].name() = " << patches[patchi].name() << endl;
    scalarField& ghfbi = ghf.boundaryFieldRef()[patchi];
    forAll(ghfbi, facei)
    {
        const label globalFacei = facei + patches[patchi].start();
        const label ownerCell(own[globalFacei]);
        if (isSurfCell[ownerCell])
        {
            const point Cio = advector.surf().centre()[ownerCell];
            ghfbi[facei] = (g.value() & Cio) - ghRef.value();
        }
    }
}

// Setting ghf on processor patches
if (Pstream::parRun())
{
    DynamicList<label> neighProcs;
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send
    for (const polyPatch& pp : patches)
    {
        const auto* ppp = isA<processorPolyPatch>(pp);

        if (ppp && pp.nPoints())
        {
            const auto& procPatch = *ppp;
            const label nbrProci = procPatch.neighbProcNo();

            neighProcs.append(nbrProci);
            UOPstream toNbr(nbrProci, pBufs);


            const scalarField& ghfbi = ghf.boundaryField()[procPatch.index()];
            List<label> interfaceFaces;
            List<scalar> interfaceFaceValues;
            forAll(ghfbi, facei)
            {
                const label globalFacei = facei + procPatch.start();
                const label ownerCell(own[globalFacei]);
                if (isSurfCell[ownerCell])
                {
                    interfaceFaces.append(facei);
                    interfaceFaceValues.append(ghfbi[facei]);
                }
            }

            /*
            Pout << "Processor " << Pstream::myProcNo() <<
                " sending interfaceFaces = " << interfaceFaces <<
                " and interfaceFaceValues = " << interfaceFaceValues <<
                " to neighbour processor = " << nbrProci << endl;
            */

            toNbr << interfaceFaces << interfaceFaceValues;
        }
    }

    // Limited to involved neighbour procs
    pBufs.finishedNeighbourSends(neighProcs);


    // Receive and combine
    for (const polyPatch& pp : patches)
    {
        const auto* ppp = isA<processorPolyPatch>(pp);

        if (ppp && pp.nPoints())
        {
            const auto& procPatch = *ppp;

            const label nbrProci = procPatch.neighbProcNo();

            UIPstream fromNeighb(nbrProci, pBufs);
            List<label> neighbourInterfaceFaces;
            List<scalar> neighbourInterfaceFaceValues;

            fromNeighb >> neighbourInterfaceFaces >> neighbourInterfaceFaceValues;

            /*
            Pout << "Processor " << Pstream::myProcNo() <<
                " received  neighbourInterfaceFaces = " <<  neighbourInterfaceFaces <<
                " and neighbourInterfaceFaceValues = " << neighbourInterfaceFaceValues <<
                " from neighbour processor: " << nbrProci << endl;
            */

            scalarField& ghfbi = ghf.boundaryFieldRef()[procPatch.index()];

            forAll(neighbourInterfaceFaces, fi)
            {
                const label facei = neighbourInterfaceFaces[fi];
                ghfbi[facei] += neighbourInterfaceFaceValues[fi];
                const label globalFacei = facei + procPatch.start();
                const label ownerCell(own[globalFacei]);
                if (isSurfCell[ownerCell])
                {
                    ghfbi[facei] *= 0.5;
                }
            }
        }
    }
}
