/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
    A matrix (A_{ij}) defines relations between the row indices as follows: row +i+ is linked to row +j+, written +i+ -> +j+, if there is +k+ such that A_{ik}*A_{jk} is non-zero.
*/


//==============================================================================
intrinsic LinkageMatrix(A::Mtrx) -> Mtrx
/*
    Intrinsic: LinkageMatrix

    The linkage matrix of a matrix.

    Declaration:
        :intrinsic LinkageMatrix(A::Mtrx) -> Mtrx

    Parameters:
        A - A matrix.

    Description:
        The (i,j)-entry of the linkage matrix +L+ of +A+ is non-zero if and only if +i+ is linked to +j+, i.e., +A_{ik}*A_{jk}+ is non-zero for some +k+.
*/
{}

    L := ZeroMatrix(BaseRing(A), Nrows(A), Nrows(A));

    for i:=1 to Nrows(A) do
        for j:=1 to Nrows(A) do
            if not { A[i][k] * A[j][k] : k in {1..Ncols(A)} } eq {0} then
                L[i][j] := 1;
            end if;
        end for;
    end for;

    return L;

end intrinsic;

//==============================================================================
intrinsic Families(M::Mtrx) -> SetEnum
/*
    Intrinsic: Families

    The families with respect to the linkage relation.

    Declaration:
        :intrinsic Families(M::Mtrx) -> SetEnum

    Parameters:
        A - A matrix.
*/
{}

    fams := {};
    for i:=1 to Nrows(M) do
        fam := Support(M[i]);
        glued := false;
        for X in fams do
            if not IsEmpty(X meet fam) then
                fams diff:={X};
                fams join:={X join fam};
                glued := true;
                break;
            end if;
        end for;
        if not glued then
            fams join:={fam};
        end if;
    end for;

    return fams;

end intrinsic;
