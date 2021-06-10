/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

//============================================================================
intrinsic RepTheory(H::AlgCheRes : GeneratorSets:=[{}], Rounds:=2, pExclude:={}, pRange:="Automatic", ParameterRange:="Automatic") -> List
{}

    G := H`Group;
    Representations(~G);

    print "Computing Verma modules";
    V := [* GradedModuleOld(VermaModule(H, rho)) : rho in G`Representations[0] *];

    print "Decomposing modules";
    res := RepTheory(H, V : GeneratorSets:=GeneratorSets, Rounds:=Rounds, pExclude:=pExclude, pRange:=pRange, ParameterRange:=ParameterRange);

    return res;

end intrinsic;

//============================================================================
intrinsic RepTheory(H::AlgCheRes, V::List : GeneratorSets:=[{}], Rounds:=2, pExclude:={}, pRange:="Automatic", ParameterRange:="Automatic") -> List
{}

    b,L,D,P := HeadsOfLocalModules(V:GeneratorSets:=GeneratorSets, Rounds:=Rounds, pRange:=pRange, pExclude:=pExclude, ParameterRange:=ParameterRange);

    if not b then
        error "Not successful";
    end if;

    dims := [ Dimension(L[i]) : i in [1..#L] ];
    Pseries := [ PoincareSeries(L[i]) : i in [1..#L]];

    print "Computing graded G-module structures via specialization";
    Lspec := [* Specialize(L[i], <P[1]>) : i in [1..#L] *];

    G := H`Group;
    d := Dimension(G);

    Gstruct := [* DecompositionInGradedGrothendieckGroup(Lspec[i], G, Reverse([1..Ngens(G)])) : i in [1..#L] *];

    return [* P, dims, Pseries, D, Families(D), Gstruct, L *];

end intrinsic;

//============================================================================
intrinsic Gordon(G::GrpMat) -> Rec
{}

    return CHAMP_GetFromDB(G`DBDir, "Cherednik/Gordon");

end intrinsic;
