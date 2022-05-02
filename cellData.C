const faceList & ff = mesh.faces(); 
const pointField & pp = mesh.points(); 

forAll ( mesh.C(), celli) 
{ 
    const cell & cc = mesh.cells()[celli]; 
    labelList pLabels(cc.labels(ff)); 
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi) 
        pLocal[pointi] = pp[pLabels[pointi]]; 
        
    scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1)); 
    
} 