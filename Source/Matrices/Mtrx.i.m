/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Simple extensions for matrices.
*/


//============================================================================
intrinsic TensorProductMixedType(A::Mtrx, B::MtrxSprs : Rep:="Dense") -> Mtrx
{}

	if Rep eq "Dense" then
		return TensorProduct(A,Matrix(B));
	else
		return TensorProduct(SparseMatrix(A),SparseMatrix(B));
	end if;

end intrinsic;

//============================================================================
intrinsic TensorProductMixedType(A::MtrxSprs, B::Mtrx : Rep:="Dense") -> Mtrx
{}

	if Rep eq "Dense" then
		return TensorProduct(Matrix(A),B);
	else
		return TensorProduct(A,SparseMatrix(B));
	end if;

end intrinsic;

//============================================================================
intrinsic TensorProductMixedType(A::Mtrx, B::Mtrx : Rep:="Dense") -> Mtrx
{}

	if Rep eq "Dense" then
		return TensorProduct(A,B);
	else
		return TensorProduct(SparseMatrix(A),SparseMatrix(B));
	end if;

end intrinsic;

//============================================================================
intrinsic TensorProductMixedType(A::MtrxSprs, B::MtrxSprs : Rep:="Dense") -> Mtrx
{}

	if Rep eq "Dense" then
		return TensorProduct(Matrix(A),Matrix(B));
	else
		return TensorProduct(A,B);
	end if;

end intrinsic;


//==============================================================================
intrinsic ConjugateTranspose(M::AlgMatElt) -> AlgMatElt
/*
    Replaces Magma's ConjugateTranspose as this is a bit too restrictive.

    History:
        Wednesday, May 22, 2013 4:34:33 PM: Initial.
*/
{The complex conjugate transpose of the matrix M.}

    M := Transpose(M);
    for i:=1 to Nrows(M) do
        for j:=1 to Ncols(M) do
            M[i,j] := ComplexConjugate(M[i][j]);
        end for;
    end for;
    return M;

end intrinsic;

//==============================================================================
intrinsic ConjugateTranspose(M::GrpMatElt) -> AlgMatElt
/*
    History:
        Wednesday, May 22, 2013 4:44:05 PM: Initial.
*/
{The complex conjugate transpose of the matrix M.}

    return ConjugateTranspose(Matrix(M));

end intrinsic;

//==============================================================================
intrinsic IsQuadratic(M::Mtrx) -> BoolElt
/*
    History:
        Wednesday, May 22, 2013 4:31:06 PM: Initial.
*/
{True iff M is a quadratic matrix.}

    return Ncols(M) eq Nrows(M);

end intrinsic;

//==============================================================================
intrinsic IsUnitary(M::AlgMatElt) -> BoolElt
/*
    History:
        Wednesday, May 22, 2013 4:31:49 PM: Initial.
*/
{True iff M is a unitary matrix.}

    if IsQuadratic(M) eq false then
        return false;
    end if;

    R := BaseRing(M);
    n := Nrows(M);
    I := IdentityMatrix(R,n);

    return (M*ConjugateTranspose(M) eq I and ConjugateTranspose(M)*M eq I);

end intrinsic;

//==============================================================================
intrinsic IsUnitary(M::GrpMatElt) -> BoolElt
/*
    History:
        Wednesday, May 22, 2013 4:44:50 PM: Initial.
*/
{True iff M is a unitary matrix.}

    return IsUnitary(Matrix(M));

end intrinsic;

//==============================================================================
intrinsic Specialize(M::Mtrx, P::RngOrdIdl) -> Mtrx
/*
    History:
        Friday, June 21, 2013 11:26:02: Initial.
*/
{Specializes (reduces) the Matrix M modulo P.}

    k,q := ResidueClassFieldGF(P);
    Mspec := ZeroMatrix(k, Nrows(M), Ncols(M));
    for i:=1 to Nrows(M) do
        for j in Support(M[i]) do
            Mspec[i,j] := q(M[i,j]);
        end for;
    end for;

    return Mspec;

end intrinsic;

//==============================================================================
intrinsic Specialize(M::Mtrx, p::RngIntElt) -> GrpMat
/*
    History:
        Friday, June 21, 2013 11:12:01: Initial.
*/
{Specializes (reduces) the matrix M in a prime ideal over the prime number p.}

    k,q := ResidueClassFieldGF(p);
    Mspec := ZeroMatrix(k, Nrows(M), Ncols(M));
    for i:=1 to Nrows(M) do
        for j in Support(M[i]) do
            Mspec[i,j] := q(M[i,j]);
        end for;
    end for;

    return Mspec;

end intrinsic;


//==============================================================================
intrinsic Entries(M::Mtrx) -> SetEnum
/*
    History:
        Monday, September 23, 2013 16:22:31: Initial.
*/
{}

    return SequenceToSet(Eltseq(M));

end intrinsic;
