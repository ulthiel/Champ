/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


declare attributes AlgChe:
	PoissonMatrix; //the matrix ({z_i,z_j})_{i,j}

//============================================================================
intrinsic PoissonMatrix(~H::AlgChe)
{
	The matrix with entries the Poisson brackets of the center generators of H.
}

	if assigned H`PoissonMatrix then
		return;
	end if;

	CenterGenerators(~H);
	CenterSpace(~H);
	R := H`CenterSpace;
	N := #H`CenterGenerators;
	M := ZeroMatrix(R, N, N);
	for i:=1 to N do
		for j:=i+1 to N do
			print "Computing {z"*Sprint(i)*",z"*Sprint(j)*"}";
			CenterGeneratorsPoissonBracket(~H, i,j);
			print "Computing preimage of {z"*Sprint(i)*",z"*Sprint(j)*"}";
			f := Preimage(H, H`CenterGeneratorsPoissonBrackets[<i,j>]);
			M[i,j] := f;
			M[j,i] := -f;
		end for;
	end for;

	H`PoissonMatrix := M;

end intrinsic;

//============================================================================
intrinsic RankStratification(M::AlgMatElt[RngMPol]) -> SeqEnum
{
	Returns a list with entries of the form <d,n> meaning that in the rank stratification of M there are n strata of dimension d.
}
	strata := [];
	n := Nrows(M);
	m := Ncols(M);
	o := Minimum(m,n);
	R := BaseRing(M);



	return strata;

end intrinsic;
