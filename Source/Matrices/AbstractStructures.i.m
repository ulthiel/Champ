/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    In my CHAMP paper (LMS J. Comput. Math., 2015) I have defined the notion of
    abstract structures of matrices. Here, we provide some elementary intrinsics
    for this notion.
*/

//==============================================================================
intrinsic AbstractStructure(M::Mtrx : entries:=[]) -> Tup
/*
    Intrinsic: AbstractStructure

    The abstract structure of a matrix.

    Declaration:
        :intrinsic AbstractStructure(M::Mtrx : entries:=[]) -> Tup

    Parameters:
        A - A matrix.
*/
{}

    if IsEmpty(entries) then
        entries := RemoveDuplicates(Eltseq(M));
        pos := Position(entries,0);
        if pos gt 0 then
            Remove(~entries,pos);
        end if;
    end if;
    abs := ZeroMatrix(Integers(), Nrows(M), Ncols(M));
    for i:=1 to Nrows(M) do
        for j in Support(M[i]) do
            abs[i][j] := Position(entries, M[i][j]);
        end for;
    end for;
    return <abs, entries>;

end intrinsic;
