/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	A extended associative array type, called dictionary. Supports sorted keys array for faster lookup.
*/

//========================================================================

declare type Dict;

declare attributes Dict:
	Keys,
	Values,
	Sorted;

//============================================================================
intrinsic Dictionary(U::. : Sorted:=false) -> Dict
{}

	D := New(Dict);

	D`Keys := [U|];
	D`Values := [];
	D`Sorted := Sorted;

	return D;

end intrinsic;

//============================================================================
intrinsic Print(D::Dict)
{}

	printf "Dictionary with universe %o", Universe(D`Keys);

end intrinsic;

//============================================================================
intrinsic Keys(D::Dict) -> SeqEnum
{}

	return D`Keys;

end intrinsic;

//============================================================================
intrinsic IsDefined(D::Dict, k::.) -> BoolElt
{}

	return k in D`Keys;

end intrinsic;

//============================================================================
intrinsic Get(D::Dict, k::.) -> .
{}

	return D`Keys[Position(D`Keys, k)];

end intrinsic;

//============================================================================
intrinsic Insert(~D::Dict, k::., v::.)
{}

	Append(~D`Keys, k);
	Append(~D`Values, v);

	if D`Sorted then
		ParallelSort(~D`Keys, ~D`Values);
	end if;

end intrinsic;
