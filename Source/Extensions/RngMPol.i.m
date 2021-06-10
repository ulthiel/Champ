/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple extensions for multivariate polynomials.
*/


declare attributes RngMPol:
    Generators,
    One,
    Zero,
    Depth,
    RegularSequence,
    IsCohenMacaulay,
    Dimension;

//==============================================================================
intrinsic Eltseq(f::RngMPolElt) -> SeqEnum
/*
    Intrinsic: Eltseq

    The sequence of variables defining a monomial.

    Declaration:
        :intrinsic Eltseq(f::RngMPolElt) -> SeqEnum

    Parameters:
    	f - A monomial of a multivariate polynomial ring

    Description:
        If +f = x_1^{r_1} ... x_n^{r_n}+ this functions returns the sequence of indices of variables giving +f+, i.e., +[1,1,...,1 (r_1 times), ..., n, n, ..., n (r_n) times]+. Applying +&*+ to this sequence gives back +f+.

*/
{}
	if not IsMonomial(f) then
		error "Eltseq only works for monomials";
	end if;

    exp := Exponents(f);
    return [ j : i in [1..exp[j]], j in [1..#exp] ];

end intrinsic;

//==============================================================================
/*
    Intrinsic: Generators

    Attaches the sequence of generators.

    Declaration:
        :intrinsic Generators(~P::RngMPol)

    Parameters:
    	P - A multivariate polynomial ring.

    Description:
        Attaches the sequence of generators to the corresponding attribute of +P+. This allows faster acces than using the dot operator.

*/
intrinsic Generators(~P::RngMPol)
{}

    if assigned P`Generators then
        return;
    else
        P`Generators := [ P.i : i in [1..Ngens(P)]];
    end if;

end intrinsic;

//==============================================================================
intrinsic GetVariableNumber(v::RngMPolElt) -> RngIntElt
/*
    Description:
        If +v+ is one of the indeterminates of a polynomial ring, return its position in the generator list.
*/
{}

    return Position(Exponents(v),1);

end intrinsic;


//==============================================================================
intrinsic SplitLinearLeft(f::RngMPolElt) -> RngMPolElt, RngMPolElt
/*
    Description:
        If +f+ is a monomial, splits off the left most linear submonomial.
*/
{}
    P := Parent(f);

    if Degree(f) eq 0 then
        return One(P), One(P);
    elif Degree(f) eq 1 then
        return f, One(P);
    else
        exp := Exponents(f);
        n:=1;
        while exp[n] eq 0 do
            n := n + 1;
        end while;
        rightpartexp := exp;
        rightpartexp[n] := rightpartexp[n] - 1;
        rightpart := Monomial(P, rightpartexp);
        leftpart := P.n;
        return leftpart, rightpart;
    end if;

end intrinsic;

//==============================================================================
intrinsic SplitLinearRight(f::RngMPolElt) -> RngMPolElt, RngMPolElt
/*
    Description:
        If +f+ is a monomial, splits off the right most linear submonomial.
*/
{}
    P := Parent(f);

    if Degree(f) eq 0 then
        return One(P), One(P);
    elif Degree(f) eq 1 then
        return f, One(P);
    else
        exp := Exponents(f);
        n:=Rank(P);
        while exp[n] eq 0 do
            n := n - 1;
        end while;
        leftpartexp := exp;
        leftpartexp[n] := leftpartexp[n] - 1;
        leftpart := Monomial(P, leftpartexp);
        rightpart := P.n;
        return leftpart, rightpart;
    end if;

end intrinsic;


//============================================================================
intrinsic RegularSequence(~I::RngMPol)
{Computes and attached a maximal regular sequence in I.}

	if assigned I`RegularSequence then
		return;
	end if;

	I`RegularSequence := RegularSequence(I);

end intrinsic;

//============================================================================
intrinsic Depth(~I::RngMPol)
{The length of one (any) maximal regular sequence in I.}

	if assigned I`Depth then
		return;
	end if;

	RegularSequence(~I);
	I`Depth := #I`RegularSequence;

end intrinsic;

//============================================================================
intrinsic Depth(I::RngMPol) -> RngIntElt
{The length of one (any) maximal regular sequence in I.}

	Depth(~I);
	return I`Depth;

end intrinsic;


//============================================================================
intrinsic Dimension(~I::RngMPol)
{The Krull dimension of I.}

	if assigned I`Dimension then
		return;
	end if;

	I`Dimension := Dimension(I);

end intrinsic;

//============================================================================
intrinsic IsCohenMacaulay(~I::RngMPol)
{True iff the depth of I is equal to the Krull dimension of I.}

	if assigned I`IsCohenMacaulay then
		return;
	end if;

	Depth(~I);
	Dimension(~I);
	I`IsCohenMacaulay := I`Depth eq I`Dimension;

end intrinsic;

//============================================================================
intrinsic IsCohenMacaulay(I::RngMPol) -> BoolElt
{True iff the depth of I is equal to the Krull dimension of I.}

	IsCohenMacaulay(~I);
	return I`IsCohenMacaulay;

end intrinsic;

//============================================================================
