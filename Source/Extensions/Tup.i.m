/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

intrinsic '&+'(x::Tup) -> .
{}

	y := 0;
	for i:=1 to #x do
		y +:= x[i];
	end for;
	return y;

end intrinsic;
