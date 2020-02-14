/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


declare attributes AlgAss:
	VectorSpace,
	VectorSpaceMap,
	NaturalModule,
	SimpleModules,
	Center,
	CenterBasis,
	CentralCharacters,
	CentralCharacterDifferences;

//============================================================================
intrinsic VectorSpace(~A::AlgAss)
{}

	if assigned A`VectorSpace then
		return;
	end if;

	A`VectorSpace, A`VectorSpaceMap := VectorSpace(A);

end intrinsic;

//============================================================================
intrinsic NaturalModule(~A::AlgAss)
{}

	if assigned A`NaturalModule then
		return;
	end if;

	VectorSpace(~A);
	prods := BasisProducts(A : Rep:="Dense");
	opmats := [ Matrix([prods[i][j] : j in [1..Dimension(A)]]) : i in [1..Dimension(A)] ];
	A`NaturalModule := RModule(opmats);

end intrinsic;

//============================================================================
intrinsic SimpleModules(~A::AlgAss)
{}

	if assigned A`SimpleModules then
		return;
	end if;

	NaturalModule(~A);
	A`SimpleModules := Constituents(A`NaturalModule);

end intrinsic;

//============================================================================
intrinsic Center(~A::AlgAss)
{}

	if assigned A`Center then
		return;
	end if;

	A`Center := Center(A);
	A`CenterBasis := [ A!c : c in Basis(A`Center) ];

end intrinsic;

//============================================================================
intrinsic CentralCharacters(~A::AlgAss)
{}

	if assigned A`CentralCharacters then
		return;
	end if;

	Center(~A);
	SimpleModules(~A);
	A`CentralCharacters := [];
	for i:=1 to #A`SimpleModules do
		S:=A`SimpleModules[i];
		val := [ ];
		for j:=1 to #A`CenterBasis do
			c:=A`CenterBasis[j];
			cS := &+[c[k]*ActionGenerator(S,k) : k in [1..Dimension(A)] ];
			if not IsDiagonal(cS) then
				error "Non-split simple module";
			end if;
			Append(~val, cS[1][1]);
		end for;
		Append(~A`CentralCharacters, val);
	end for;

end intrinsic;

//============================================================================
intrinsic CentralCharacterDifferences(~A::AlgAss)
{}

	if assigned A`CentralCharacterDifferences then
		return;
	end if;

	A`CentralCharacterDifferences := [];
	for i:=1 to #A`SimpleModules do
		for j:=i+1 to #A`SimpleModules do
			d := [ A`CentralCharacters[i][k] - A`CentralCharacters[j][k] : k in [1..#A`CenterBasis] ];
			Append(~A`CentralCharacterDifferences, d);
		end for;
	end for;

end intrinsic;
