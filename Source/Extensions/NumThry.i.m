/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


/*
	Simple number theoretic extensions.
*/

//==============================================================================
intrinsic pPart(n::RngIntElt, p::RngIntElt) -> RngIntElt
/*
    Intrinsic: pPart

    pPart of a number

    Declaration:
        :intrinsic pPart(n::RngIntElt, p::RngIntElt) -> RngIntElt

    Description:
        Returns +p+ ^ +a+ with +a+ maximal such that +p+ ^ +a+ divides +n+.
*/
{}

  return p^pValuation(n, p);

end intrinsic;

//==============================================================================
intrinsic pPart(G::Grp, p::RngIntElt) -> RngIntElt
/*
    Intrinsic: pPart

    pPart of a the group order

    Declaration:
        :intrinsic pPart(G::Grp, p::RngIntElt) -> RngIntElt

    Description:
        Returns the +p+-Part of the order of +G+.
*/
{}

  return p^pValuation(G, p);

end intrinsic;


//==============================================================================
intrinsic pValuation(n::RngIntElt, p::RngIntElt) -> RngIntElt
/*
    Intrinsic: pValuation

    pPart of a the group order

    Declaration:
        :intrinsic pValuation(n::RngIntElt, p::RngIntElt) -> RngIntElt

    Description:
        Returns the maximal +a+ such that +p+ ^ +a+ divides +n+.
*/
{The highest power of p dividing n.}

    fact := Factorization(n);
    res := exists(pos){ i : i in [1..#fact] | fact[i][1] eq p };
    if res eq false then
        return 0;
    else
        return fact[pos][2];
    end if;

end intrinsic;

//==============================================================================
intrinsic pValuation(G::Grp, p::RngIntElt) -> RngIntElt
/*
    Intrinsic: pValuation

    pValuation of the group order

    Declaration:
        :intrinsic pValuation(G::Grp, p::RngIntElt) -> RngIntElt

    Description:
        Returns the +p+-valuation of the order of +G+.
*/
{The highest power of p dividing the order of G.}

    return pValuation(Order(G),p);

end intrinsic;



//==============================================================================
/*
    Intrinsic: ResidueClassFieldGF

    The residue class field of the prime ideal in type GF.

    Declaration:
        :intrinsic ResidueClassFieldGF(P::RngOrdIdl) -> FldFin, Map
        and
        :intrinsic ResidueClassFieldGF(p::RngIntElt) -> FldFin, Map

    Parameters:
       P - A prime ideal

       or

       p - A prime number

    Description:
         If +P+ is a prime ideal of a ring of integers +O+, then the residue field +O/P+ together with the quotient morphism is returned. The residue field here is of type *GF*. This fixes a weird behavior in Magma, namely the residue class field of a prime ideal is not of the form GF, but abstractly a finite field. One first has to tell Magma to embed this into GF.
*/
intrinsic ResidueClassFieldGF(P::RngOrdIdl) -> FldFin, Map
{}

    L,quotmor := ResidueClassField(P);
    q := #L;
    Lfixed := GF(q);
    Embed(L,Lfixed);
    return Lfixed, quotmor;

end intrinsic;

intrinsic ResidueClassFieldGF(p::RngIntElt) -> FldFin, Map
{}

    L, quotmor := ResidueClassField(p);
    return L, quotmor;

end intrinsic;


//==============================================================================
intrinsic IsNormal(K::FldRat) -> BoolElt
/*
    Intrinsic: IsNormal

    Only return true for the rational field.

    Declaration:
        :intrinsic IsNormal(K::FldRat) -> BoolElt

    Parameters:
       K - The rational field.

    Description:
         +IsNormal+ is not implemented in Magma for rational fields (the answer is obviously yes) but this check can occur in sequential checks of number fields.
*/
{}

    return true;

end intrinsic;



//==============================================================================
/*
    Intrinsic: EmbeddingMap

    Several embedding maps for field extensions

    Declaration:
        :intrinsic EmbeddingMap(L::FldRat, K::FldAlg) -> Map
        :intrinsic EmbeddingMap(L::FldRat, K::FldFunRat) -> Map
        :intrinsic EmbeddingMap(L::FldAlg, K::FldFunRat) -> Map

    Parameters:
       L - A field.
       K - An extension field of +K+.

    Description:
         Returns the embedding of a field +L+ into an extension field +K+ for several types not supported by Magma for some reason.
*/
intrinsic EmbeddingMap(L::FldRat, K::FldAlg) -> Map
/*
    Returns the embedding of the rational numbers L into K.
*/
{}

	return hom<Rationals() -> K|>;

end intrinsic;

intrinsic EmbeddingMap(L::FldRat, K::FldFunRat) -> Map
/*
    Returns the embedding of the rational numbers L into K.
*/
{}

	return hom<Rationals() -> K|>;

end intrinsic;

intrinsic EmbeddingMap(L::FldAlg, K::FldFunRat) -> Map
/*
    Return the embedding of L into the scalars of the rational function field K over an extension field of L.
*/
{}
	Kb := BaseRing(K);
	LKbemb := EmbeddingMap(L,Kb);
	KbKemb := hom<Kb->K | [Kb.i*One(K) : i in [1..Ngens(Kb)]]>;
	LKemb := hom<L->K | [KbKemb(LKbemb(L.i)) : i in [1..Ngens(L)]]>;

	return LKemb;

end intrinsic;



//==============================================================================
/*
    Intrinsic: MinimalCyclotomicField

    A minimal cyclotomic field containing all given elements or a number field.

    Declaration:
        :intrinsic MinimalCyclotomicField(S::Setq[RngIntElt]) -> FldRat
        :intrinsic MinimalCyclotomicField(K::FldNum) -> FldCyc
        :intrinsic MinimalCyclotomicField(K::FldRat) -> FldRat

    Parameters:
        S - A sequence of integers

        or

        K - A number field or the rational field.

    Description:
        Returns a minimal cyclotomic field containing all given elements or a number field.

    History:
        * Tuesday, October 01, 2013 22:16:04: Added for FldRat.
        * Tuesday, October 01, 2013 22:08:17: Added for FldNum.
        * Friday, June 28, 2013 18:23:50: Initial.

    To-do:
        The intrinsic for type +FldNum+ is implemented in a very naive way and should be fixed. It works fine for small orders though.
*/


//==============================================================================
intrinsic MinimalCyclotomicField(S::Setq[RngIntElt]) -> FldRat
/*
    The minimal cyclotomic field containing a set of integers is simply the rational field.
*/
{}

    return Rationals();

end intrinsic;


//==============================================================================
intrinsic MinimalCyclotomicField(K::FldRat) -> FldRat
{}

    return Rationals();

end intrinsic;


//==============================================================================
intrinsic MinimalCyclotomicField(K::FldNum) -> FldCyc
{}

    //OH, I don't know how to get this theoretically. This should be fixed!
    d := Degree(K);
    n:=1;
    while EulerPhi(n) lt d do
        n +:= 1;
    end while;

    while not K subset CyclotomicField(n) do
        n +:= 1;
    end while;

    return CyclotomicField(n);

end intrinsic;




//==============================================================================
/*
    Intrinsic: CommonOverfield

    A field containing a given list of fields.

    Declaration:
        :intrinsic CommonOverfield(K::FldRat, L::FldRat) -> FldRat
        :intrinsic CommonOverfield(K::FldCyc, L::FldCyc) -> FldCyc
        :intrinsic CommonOverfield(K::FldRat, L::FldCyc) -> FldCyc
        :intrinsic CommonOverfield(K::FldCyc, L::FldRat) -> FldCyc
        :intrinsic CommonOverfield(K::FldNum, L::FldNum) -> FldCyc
        :intrinsic CommonOverfield(Q::List) -> Fld

    Parameters:
        K - A field
        L - A field
        Q - A list of fields.

    Description:
        Returns a field containing a given list of fields.

    History:
        * Wednesday, October 09, 2013 20:25:05: Added for rational fields.
        * Saturday, September 21, 2013 15:26:04: Added for list of fields.
        * Thursday, September 19, 2013 11:21:16: Added for both fields of type +FldNum+.
        * Tuesday, September 17, 2013 13:26:38: Initial.

*/

//==============================================================================
intrinsic CommonOverfield(K::FldRat, L::FldRat) -> FldRat
{}

    return Rationals();

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldCyc, L::FldCyc) -> FldCyc
/*
    Given two cyclotomic fields K and L return a common cyclotomic overfield.
*/
{}
    return MinimalCyclotomicField(Generators(K) join Generators(L));

end intrinsic;


//==============================================================================
intrinsic CommonOverfield(K::FldRat, L::FldCyc) -> FldCyc
/*
    Given two cyclotomic fields K and L return a common cyclotomic overfield.
*/
{}
    return L;

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldCyc, L::FldRat) -> FldCyc
/*
    Given two cyclotomic fields K and L return a common cyclotomic overfield.
*/
{}
    return K;

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldNum, L::FldNum) -> FldCyc
{}
    if K eq L then
        return K;
    elif K subset L then
        return L;
    elif L subset K then
        return K;
    else
        return Compositum(K,L);
    end if;

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldNum, L::FldRat) -> FldCyc
{}

    return K;

end intrinsic;


//==============================================================================
intrinsic CommonOverfield(K::FldRat, L::FldNum) -> FldCyc
{}

    return L;

end intrinsic;


//==============================================================================
intrinsic CommonOverfield(K::FldCyc, L::FldFunRat[FldCyc]) -> FldFunRat
{}

    M := CommonOverfield(K, BaseRing(L));
    return ChangeRing(L, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldFunRat[FldCyc], L::FldCyc) -> FldFunRat
{}

    M := CommonOverfield(L, BaseRing(K));
    return ChangeRing(K, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldRat, L::FldFunRat[FldCyc]) -> FldFunRat
{}

    M := CommonOverfield(K, BaseRing(L));
    return ChangeRing(L, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldFunRat[FldCyc], L::FldRat) -> FldFunRat
{}

    M := CommonOverfield(L, BaseRing(K));
    return ChangeRing(K, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldRat, L::FldFunRat[FldRat]) -> FldFunRat
{}

    M := CommonOverfield(K, BaseRing(L));
    return ChangeRing(L, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldFunRat[FldRat], L::FldRat) -> FldFunRat
{}

    M := CommonOverfield(L, BaseRing(K));
    return ChangeRing(K, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(K::FldFunRat[FldCyc], L::FldFunRat[FldCyc]) -> FldFunRat
{}

    M := CommonOverfield(BaseRing(K), BaseRing(L));
    return ChangeRing(K, M);

end intrinsic;

//==============================================================================
intrinsic CommonOverfield(Q::List) -> Fld
/*
    If Q is a list of fields, return a common overfield for all these fields if possible.
*/
{}
    K := Q[1];
    for i:=2 to #Q do
        K := CommonOverfield(K, Q[i]);
    end for;

    return K;

end intrinsic;



//==============================================================================
intrinsic HasRootOfUnity(n::RngIntElt, K::Fld) -> BoolElt, FldElt
/*
    Intrinsic: HasRootOfUnity

    Checks if a field has a root of unity.

    Declaration:
        :intrinsic HasRootOfUnity(n::RngIntElt, K::Fld) -> BoolElt, FldElt

    Parameters:
        K - A field.
        n - An integer.

    Description:
        Checks if +K+ contains a primitive +n+-th root of unity. If so, the one fixed by Magma is returned.

*/
{}

    res := true;
    try
        zeta := RootOfUnity(n, K);
    catch e
        res := false;
    end try;
    if not res then
        return false, _;
    end if;
    if zeta notin K then
        return false, _;
    else
        return true, zeta;
    end if;

end intrinsic;
