/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Simple extensions for sparse matrices.
*/

//============================================================================
intrinsic SparseMatrix(M::MtrxSprs) -> MtrxSprs
/*
	Just for completeness.
*/
{}

	return M;

end intrinsic;

//============================================================================
intrinsic ChangeRing(M::MtrxSprs, phi::Map) -> MtrxSprs
{
	Ring change of M using phi.
}

    N := SparseMatrix(Codomain(phi), Nrows(M), Ncols(M));
    for i:=1 to Nrows(M) do
        for j in Support(M[i]) do
            N[i][j] := phi(M[i][j]);
        end for;
    end for;

    return N;

end intrinsic;

//============================================================================
intrinsic Column(M::MtrxSprs, j::RngIntElt) -> ModTupFldElt
{
	Returns the j-th column of M.
}

    return Transpose(M)[j];

end intrinsic;
