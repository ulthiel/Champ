/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
    Intrinsics around Calogero-Moser cells.

    History:
        * Monday, April 14, 2014 at 12:15:05: Initial.
*/


//==============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map, v::ModTupFldElt, vstar::ModTupFldElt, y::ModTupFldElt : Complex:=false) -> AlgMatElt
{}
    if not IsRegular(G,v) then
        error "v has to be a regular vector.";
    end if;

    NumberingMap(~G);
    ReflectionLibrary(~G);
    D:=ZeroMatrix(Codomain(c), #G, #G);
    E:=VectorSpace(Codomain(c), #G);
    for i:=1 to #G do
        w := G`InverseNumberingMap(i);
        //image of e_w under D^{c,v,vstar}_y
        img := CanonicalPairing(y, vstar*Transpose(Matrix(w)))*E.i; //<y,w^-1*vstar>
        for s in G`ReflectionLibraryFlat do
            img +:= s`Eigenvalue*c(s`ReflectionClass)*CanonicalPairing(y,s`Coroot)/CanonicalPairing(v,s`Coroot)*E.G`NumberingMap(s`Element*w);
        end for;
        D[i]:=img;
    end for;

    if not Complex then
        return D;
    else
        return ChangeRing(D, ComplexField());
    end if;

end intrinsic;

//==============================================================================
intrinsic GaudinGroupElement(G::GrpMat, c::Map, v::ModTupFldElt, y::ModTupFldElt) -> AlgGrpElt
{}

    if not IsRegular(G,v) then
        error "v has to be a regular vector.";
    end if;

    B:=GroupAlgebra(Codomain(c), G);
    ReflectionLibrary(~G);
    op:=Zero(B);
    for s in G`ReflectionLibraryFlat do
        op +:= s`Eigenvalue*c(s`ReflectionClass)*CanonicalPairing(y,s`Coroot)/CanonicalPairing(v,s`Coroot)*B!s`Element;
    end for;

    return op;

end intrinsic;

//==============================================================================
intrinsic GaudinAlgebra(G::GrpMat, c::Map, v::ModTupFldElt) -> .
{}

    B:=GroupAlgebra(Codomain(c), G);
    V:=VectorSpace(G);
    return sub<B|[GaudinGroupElement(G,c,v,V.i) : i in [1..Dimension(G)]]>;

end intrinsic;

//==============================================================================
intrinsic GaudinModule(G::GrpMat, c::Map, v::ModTupFldElt, M::ModGrp) -> ModRng
{}

    V:=VectorSpace(G);

    opmats := [];
    for i:=1 to Dimension(G) do
        opmat := ZeroMatrix(Codomain(c), Dimension(M), Dimension(M));
        D:=G`GaudinGroupElement(G,c,v,V.i);
        for j:=1 to Dimension(M) do
            opmat[j] := VectorSpace(M)!(M.j*D);
        end for;
        Append(~opmats, opmat);
    end for;

    return RModule(opmats);

end intrinsic;


//==============================================================================
intrinsic GaudinSpectrum(G::GrpMat, c::Map, v::ModTupFldElt, vstar::ModTupFldElt) -> ModTupFld
{}

    formspace:=KSpace(Codomain(c), Dimension(G));
    D := [ GaudinOperator(G,c,v,vstar, VectorSpace(G).i) : i in [1..Dimension(G)] ];

    eigenvalues := [ SetToSequence(Eigenvalues(D[i])) : i in [1..#D] ];
    eigenspaces := [ [Eigenspace(D[i],eigenvalues[i][j][1]) : j in [1..#eigenvalues[i]]] : i in [1..#D]];
    commoneigenspaces := [**];
    commoneigenforms := [];
    for C in CartesianProduct([{1..#eigenvalues[i]} : i in [1..Dimension(G)]]) do
        B := eigenspaces[1][C[1]];
        for i:=2 to Dimension(G) do
            B meet:= eigenspaces[i][C[i]];
        end for;
        if Dimension(B) eq 0 then
            continue;
        end if;
        Append(~commoneigenspaces, B);
        rho := formspace![eigenvalues[i][C[i]][1] : i in [1..Dimension(G)] ];
        Append(~commoneigenforms, rho);
    end for;

    return commoneigenspaces, commoneigenforms;

end intrinsic;

//==============================================================================
intrinsic GaudinSplittingField(G::GrpMat, c::Map, v::ModTupFldElt, vstar::ModTupFldElt) -> Fld
{}

    V:=VectorSpace(G);
    for i:=1 to Dimension(G) do
        GaudinOperator(~G, c, v, vstar, V.i);
    end for;

    fields := [* SplittingField(MinimalPolynomial(GaudinOperator(G,c,v,vstar,V.i))) : i in [1..Dimension(V)] *];

    return CommonOverfield(fields);

end  intrinsic;

//==============================================================================
intrinsic CMCellModule(G::GrpMat, c::Map, rho::ModTupFldElt, v::ModTupFldElt) -> ModGrp
{}

    RegularModule(~G);

    V:=VectorSpace(G);
    for i:=1 to Dimension(G) do
        GaudinOperator(~G, c, v, Zero(V), V.i);
    end for;

    kernel := KSpace(Codomain(c), #G);
    for i:=1 to Dimension(G) do
        kernel meet:= Kernel( (GaudinOperator(G,c,v,Zero(V), V.i) - rho[i]*IdentityMatrix(Codomain(c), #G))^#G );
    end for;

    return sub<ChangeRing(G`RegularModule, Codomain(c))|kernel>;

end intrinsic;

//==============================================================================
intrinsic LeftCMCellFamilies(G::GrpMat, c::Map, v::ModTupFldElt) -> SetEnum
{}

    Representations(~G);
    LiftRepresentationsToCommonBaseField(~G);
    Modules(~G);
    print "Computing Gaudin modules";
    A:=[];
    for i:=1 to #G`Modules[0] do
        Append(~A, GaudinModule(G,c,v, G`Modules[0][i]));
    end for;

    D,S := DecompositionMatrix(A);

    L := LinkageMatrix(D)^#G`Modules[0];

    families:={ Support(L[i]) : i in [1..Nrows(L)] };

    return families;

end intrinsic;


declare attributes GrpMat:
    GaudinOperators,
    GaudinGroupElement;
