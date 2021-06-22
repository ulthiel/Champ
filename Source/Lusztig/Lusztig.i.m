/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/



//=============================================================================
intrinsic BetaVector(lambda::SeqEnum[RngIntElt]) -> SeqEnum
{
	The beta vector of a partition lambda as defined by Broue-Kim [BK02], section 3.
}

	h := #lambda;
	beta := [h+lambda[i]-i : i in [1..h]];
	return beta;

end intrinsic;

//=============================================================================
intrinsic BetaVectorShifted(lambda::SeqEnum[RngIntElt], m::RngIntElt) ->
SeqEnum
{
The m-shift beta[m] of the beta vector of a partition as defined by Broue-Kim [BK02].
}

	beta := BetaVector(lambda);
	for i:=1 to #beta do
		beta[i] +:= m;
	end for;
	for i:=m-1 to 0 by -1 do
		Append(~beta, i);
	end for;

	return beta;

end intrinsic;

//=============================================================================
intrinsic Content(S::SeqEnum[SeqEnum[RngIntElt]]) -> RngUPolElt
{
The content of the symbol S. This is the polynomial sum_i n_i x^i, where n_i is the multiplicity of the entry i in S.
}

	X := FlatFixed(S);
	P<x> := PolynomialRing(Integers());
	c := 0;
	for i in SequenceToSet(X) do
		c +:= NumberOfOccurrences(X,i)*x^i;
	end for;

	return c;

end intrinsic;


//=============================================================================
intrinsic BroueKimSymbol(lambda::SeqEnum, m::SeqEnum) -> SeqEnum
{
The symbol of the multipartition lambda as defined by Broue-Kim [BK02].
}

	d := #lambda;
	rk := &+[&+lambda[i] : i in [1..d]];
	hauteur := [#lambda[i] : i in [1..d]];
	hlambda := Maximum(hauteur);
	hauteurchargee := [hauteur[i] - m[i] : i in [1..d]];
	hclambda := Maximum(hauteurchargee);

	betavectors := [BetaVectorShifted(lambda[i], hclambda-hauteurchargee[i]) : i in [1..d]];

	return betavectors;

end intrinsic;
