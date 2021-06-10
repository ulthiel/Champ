/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
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

intrinsic LinkageGraph(A::Mtrx) -> GrphUnd
{}

    L := LinkageMatrix(A);
    E := {};
    for i:=1 to Nrows(L) do
        for j in Support(L[i]) diff {i} do
            E join:={ {i,j} };
        end for;
    end for;
    return Graph< {1..Nrows(L)} | E>;

end intrinsic;

//==============================================================================
intrinsic Families(M::Mtrx) -> SetEnum, SeqEnum
/*
    Intrinsic: Families

    The families with respect to the linkage relation.

    Declaration:
        :intrinsic Families(M::Mtrx) -> SetEnum

    Parameters:
        A - A matrix.
*/
{}

    /*
    fams := {@@};
    for i:=1 to Nrows(M) do
        fam := SetToIndexedSet(Support(M[i]));
        glued := false;
        for X in fams do
            if not IsEmpty(X meet fam) then
                Diff(~fams, {@X@});
                fams join:={@X join fam@};
                glued := true;
                break;
            end if;
        end for;
        if not glued then
            fams join:={@fam@};
        end if;
    end for;
    */
    famspre := ConnectedComponents(LinkageGraph(M));
    fams := {@ {@ Index(f) : f in F@} : F in famspre @};
    sigma := Flat([IndexedSetToSequence(fam) : fam in fams]);

    return fams,sigma;

end intrinsic;

intrinsic Permute(A::Mtrx, sigma::SeqEnum) -> Mtrx
{}

    Asigma := [A[sigma[i],sigma[j]] : i in [1..Nrows(A)], j in [1..Ncols(A)]];
    return Matrix(BaseRing(A), Nrows(A),Ncols(A),Asigma);

end intrinsic;
