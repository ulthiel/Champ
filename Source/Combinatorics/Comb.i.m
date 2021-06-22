/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Intrinsics around combinatorics.
*/

//==============================================================================
//freeze;

//==============================================================================
intrinsic MaximalChains(S::SetEnum, lessequal::Intrinsic) -> SeqEnum
/*
    This is implemented in the really _naive_ way and has to be optimized!

    Status: Alpha

    Todo: Optimize

    History:
        Wednesday, May 8, 2013 2:15:41 PM: Initial.
*/
{The set of maximal chains of S with respect to the preorder described by lessequal.}

    maxs := [ [x] : x in MinimalElements(S,lessequal) ];
    while exists(i){i : i in [1..#maxs] | exists{y : y in S | lessequal(Last(maxs[i]),y) and Last(maxs[i]) ne y}} do
        x := Last(maxs[i]);
        covers := { y : y in S | lessequal(x,y) and x ne y };
        mincovers := MinimalElements(covers, lessequal);

        for y in mincovers do
            curchain := maxs[i];
            Append(~curchain, y);
            Append(~maxs, curchain);
        end for;

        Remove(~maxs, i);
    end while;

    return maxs;

end intrinsic;

//==============================================================================
intrinsic MaximalElements(S::SetEnum, lessequal::Intrinsic) -> SetEnum
/*
    This is implemented in the really _naive_ way and has to be optimized!

    Status: Alpha

    Todo: Optimize

    History:
        Wednesday, May 8, 2013 2:15:41 PM: Initial.
*/
{The maximal elements of S with respect to the preorder described by lessequal.}

    maxs := {};
    for x in S do
        if exists{y : y in S | lessequal(x,y) and x ne y} then
            continue;
        end if;
        maxs join:={x};
    end for;

    return maxs;

end intrinsic;


//==============================================================================
intrinsic MinimalElements(S::SetEnum, lessequal::Intrinsic) -> SetEnum
/*
    This is implemented in the really _naive_ way and has to be optimized!

    Status: Alpha

    Todo: Optimize

    History:
        Wednesday, May 8, 2013 2:15:41 PM: Initial.
*/
{The minimal elements of S with respect to the preorder described by lessequal.}

    mins := {};
    for x in S do
        if exists{y : y in S | lessequal(y,x) and x ne y} then
            continue;
        end if;
        mins join:={x};
    end for;

    return mins;

end intrinsic;



//==============================================================================
intrinsic LittlewoodRichardsonCoefficients(R::RngMPol, lambda::SeqEnum, mu::SeqEnum : Method:="LinearSystem") -> SetEnum
/*
    Computes the expansion of the product s_lambda*s_mu of Schur polynomials in R. If Method is "LinearSystem" this is done via a linear system (this is of course not an efficient way, it's just for educational purposes).


    Status: Alpha

    Todo: Implement Littlewood-Richardson rule.

    History:
        Friday, May 10, 2013 1:05:55 PM: Initial.
*/
{}
    n := Weight(lambda); //lambda is partition of n
    m := Weight(mu); //mu is partition of m
    r := Rank(R);
    K := BaseRing(R);

    if Method eq "LinearSystem" then
        cartprod := [ [mon[i] : i in [1..#mon]] : mon in CartesianProduct([ [0..n+m] : i in [1..r] ])];
        mons := [ Monomial(R,mon) : mon in cartprod | &+mon eq n+m];

        parts := [ nu : nu in Partitions(n+m) | #nu le r ];
        schurpols := [ SchurPolynomial(R,nu) : nu in parts ];

        V := VectorSpace(K, #mons);
        W := VectorSpace(K, #parts);

        product := SchurPolynomial(R,lambda)*SchurPolynomial(R,mu);

        b := V![MonomialCoefficient(product, m) : m in mons];

        A := ZeroMatrix(K, #parts, #mons);
        for i:=1 to #parts do
            for j:=1 to #mons do
                A[i][j] := MonomialCoefficient(schurpols[i],mons[j]);
            end for;
        end for;
        c := Solution(A,b);

        return { <parts[i],c[i]> : i in [1..#parts] | c[i] ne 0 };

    end if;


end intrinsic;


//=============================================================================
intrinsic NumberOfUnorderedFactorizations(n::RngIntElt, k::RngIntElt : Method:="Hughes-Shallit") -> RngIntElt
/*
	The number of unordered factorizations of a natural number n (see UnorderedFactorizations()).

	If Method is Hughes-Shallit then all factors are forced to be less or equal k, and a recursive formula by Hughes-Shallit is used.

	Reference: Knopfmacher-Mays, "Ordered and Unordered Factorizations of Integers" in The Mathematica Journal.

    Status: Beta

    Todo: Can be optimized perhaps.

    History:
        Saturday, May 4, 2013 1:34:29 PM: OK.
*/
{The number of unordered factorizations of a natural number n. If Method is Hughes-Shallit then all factors are forced to be less or equal k.}

	if n eq 1 then
		return 1;
	elif k eq 1 then
		return 0;
	else
		return &+[NumberOfUnorderedFactorizations(n div d, d) : d in Divisors(n) | d le k ];
	end if;

end intrinsic;

//=============================================================================
intrinsic NumberOfUnorderedFactorizations(n::RngIntElt : Method:="Hughes-Shallit") -> RngIntElt
{}

	return NumberOfUnorderedFactorizations(n,n : Method:=Method);

end intrinsic


//=============================================================================
intrinsic Partitions(S::SetEnum) -> SetEnum
/*
	We compute the set of partitions of a set S in a _naive_ way.

    Status: Alpha.

    Todo: Optimize.

    History:
        Saturday, May 4, 2013 1:33:32 PM: OK.
*/
{The set of partitions of the set S.}
	parts := {};
	for lambda in Partitions(#S) do
		P := CartesianProduct([ Subsets(S,lambda[i]): i in [1..#lambda]]);
		for p in P do
			if &join{ p[i] : i in [1..#lambda]} eq S then
				parts join:={ { p[i] : i in [1..#lambda]} };
			end if;
		end for;
	end for;

	return parts;

end intrinsic;

//==============================================================================
intrinsic SchurPolynomial(R::RngMPol, lambda::SeqEnum) -> RngMPolElt
/*
    History:
        Friday, May 10, 2013 12:51:56 PM: Initial.
*/
{The Schur polynomial s_lambda in R.}

    r := Rank(R);
    tabs := TableauxOfShape(lambda,r);
    s := Zero(R);
    for t in tabs do
        entries := FlatFixed(Eltseq(t));
        mon := One(R);
        for i:=1 to #entries do
            mon *:= R.entries[i];
        end for;
        s +:= mon;
    end for;

    return s;

end intrinsic;


//=============================================================================
intrinsic UnorderedFactorizations(n::RngIntElt, k::RngIntElt : Method:="Knopfmacher-Mays") -> SeqEnum
/*
	The set of all unordered factorizations (also called multiplicative partitions) of a natural number n, i.e., n=n_1*...*n_l. For n = 1 the empty partition [] is returned.

	If Method is Knopfmacher-Mays, then n_i <= k is forced for all i.

	Reference: Knopfmacher-Mays, "Ordered and Unordered Factorizations of Integers" in The Mathematica Journal.

    Status: Beta.

    Todo: Optimize.

    History:
        Saturday, May 4, 2013 1:35:20 PM: OK.
*/
{The set of all unordered factorizations of a natural number n, i.e., n=n_1*...*n_l. If Method is Knopfmacher-Mays, then n_i <= k is forced for all i.}

	if Method eq "Knopfmacher-Mays" then

		if n eq 1 then
			return [ [] ];
		elif k eq 1 then
			return [ ];
		elif IsPrime(n) then
			if k lt n then
				return [];
			else
				return [ [n] ];
			end if;
		else

			X := [];
			for d in Divisors(n) do
				if d gt k then
					continue;
				end if;
				for	f in UnorderedFactorizations(n div d, d) do
						Append(~X, f cat [d]);
				end for;
			end for;

			return X;

		end if;

	end if;


end intrinsic;

//=============================================================================
intrinsic UnorderedFactorizations(n::RngIntElt : Method:="Knopfmacher-Mays") -> SetEnum
/*
    The set of all unordered factorizations of a natural number n, i.e., n=n_1*...*n_l.

    Reference: Knopfmacher-Mays, "Ordered and Unordered Factorizations of Integers" in The Mathematica Journal.

    Status: Beta.

    Todo: Optimize.

    History:
        Saturday, May 4, 2013 1:35:20 PM: OK.
*/
{The set of all product partitions of a natural number n, i.e., of all factorizations n=n_1*...*n_k with n_i natural numbers greater than 1.}

	return UnorderedFactorizations(n,n : Method:=Method);

end intrinsic;


//==============================================================================
intrinsic FindPermutation(X::SeqEnum, Y::SeqEnum) -> SeqEnum
/*
    History:
        Saturday, September 14, 2013 23:42:34: Initial
*/
{If X and Y are subsets of the same set S try to find a permutation sigma on S such that Y[sigma[i]] = X[i] for all i.}

    if #X ne #Y then
        error "X and Y have different cardinality.";
    end if;

    if SequenceToSet(X) ne SequenceToSet(Y) then
        error "There is no bijection.";
    end if;

    sigma := [];
    remaining := [1..#Y];
    for i:=1 to #X do
        j := Position(Y, X[i]);
        if j eq 0 then
            error "There is no bijection.";
        end if;
        Append(~sigma, j);
    end for;

    if #SequenceToSet(sigma) ne #sigma then
        error "There is no bijection.";
    end if;

    return sigma;

end intrinsic;

//============================================================================
intrinsic NumberOfMonomials(n::RngIntElt, d::RngIntElt) -> RngIntElt
{The number of monomials of degree d in n variables.}

	return Binomial(n+d-1,d);

end intrinsic;

//============================================================================
intrinsic PochhammerSymbol(x::RngIntElt, n::RngIntElt) -> RngIntElt
{}

	return ArrayProduct([ x-i : i in [0..n-1] ]);

end intrinsic;
