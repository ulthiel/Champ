/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


//==============================================================================
intrinsic EssentialHyperplanes(G::GrpMat) -> SeqEnum
{}

    CherednikParameterSpace(~G);
    return ChangeUniverse(Sort(SetToSequence(Keys(CHAMP_GetFromDB(G`DBDir, "Hecke/RouquierFamilies")) diff {1})), G`CherednikParameterSpace);

end intrinsic;

//==============================================================================
intrinsic RouquierFamilies(G::GrpMat) -> Assoc
{}
    CherednikParameterSpace(~G);
    return CHAMP_GetFromDB(G`DBDir, "Hecke/RouquierFamilies");

end intrinsic;
