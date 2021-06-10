/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for Laurent series rings.
*/



//==============================================================================
/*
    Intrinsic: IsPalindromic

    Declaration:
        :intrinsic IsPalindromic(f::RngSerLaurElt) -> BoolElt

    Parameters:
       f - A Laurent series

    Description:
        Returns +true+ iff +f+ is polynomdromic, i.e., the sequence of coefficients can be reversed without changing +f+.
*/
intrinsic IsPalindromic(f::RngSerLaurElt) -> BoolElt
{}

    coeffs := Coefficients(f);
    return coeffs eq Reverse(coeffs);

end intrinsic;
