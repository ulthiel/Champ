/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for rational function fields.
*/



//==============================================================================
/*
    Intrinsic: ChangeRing

    Declaration:
        :intrinsic ChangeRing(F::FldFunRat, S::Rng) -> FldFunRat

    Parameters:
        F - A rational function field.
        S - A ring.

    Description:
        Change the coefficient ring of a rational function field +K+ to +S+. This is consistent with the
        corresponding function for polynomial rings.
*/
intrinsic ChangeRing(F::FldFunRat, S::Rng) -> FldFunRat
{}

    FS := RationalFunctionField(S, Ngens(F));
    AssignNames(~FS, Names(F));
    return FS;

end intrinsic;
