/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


//==============================================================================
intrinsic EssentialHyperplanes(G::GrpMat) -> SeqEnum
{}

    CherednikParameterSpace(~G);
    return ChangeUniverse(Sort(SetToSequence(Keys(CHAMP_GetFromDB("RouquierFamilies", G`DBName)) diff {1})), G`CherednikParameterSpace);

end intrinsic;

//==============================================================================
intrinsic RouquierFamilies(G::GrpMat) -> Assoc
{}
    CherednikParameterSpace(~G);
    return CHAMP_GetFromDB("RouquierFamilies", G`DBName);

end intrinsic;
