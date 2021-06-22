/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

intrinsic '&+'(x::Tup) -> .
{}

	y := 0;
	for i:=1 to #x do
		y +:= x[i];
	end for;
	return y;

end intrinsic;
