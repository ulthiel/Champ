/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple extensions for group modules.
*/

//=============================================================================
intrinsic '*'(v::ModGrpElt, f::AlgGrpElt) -> ModGrpElt
/*
    Intrinsic: '*'(v::ModGrpElt, f::AlgGrpElt)

    Action of a group ring element on a group module element.

    Declaration:
        :intrinsic '*'(v::ModGrpElt, f::AlgGrpElt) -> ModGrpElt

    Parameters:
       v - An element of a +G+-module
       f - An element of the group ring of +G+.

    Description:
        Let +M+ be an +RG+-module for a group +G+ over a commutative ring +R+. If +v+ is an element of +M+ and +f+ is an element of the group ring +RG+, then return the element of +M+ emerging from the action of +f+ on +v+.
*/
{The action of the group ring element f on the G-module element v.}

	V := Parent(v);
	w := Zero(V);
	for g in Support(f) do
		w +:= Coefficient(f,g)*(v*g);
	end for;

	return w;

end intrinsic;

//==============================================================================
intrinsic Specialize(M::ModGrp, p::RngIntElt) -> ModGrp
/*
    Intrinsic: Specialize

    Specializes a group module over a number field in a prime number.

    Declaration:
        :intrinsic Specialize(M::ModGrp, p::RngIntElt) -> ModGrp

    Parameters:
       M - A GModule
       p - A prime number

    Description:
        If +M+ is a +KG+-module over a number field +K+ a +p+ is a prime number, then randomly choose a prime ideal +P+ of the ring of integers +O+ lying above +p+ using <RandomPrimeIdeal> and reduce +M+ modulo +P+ (possibly from the localization of +O+ in +P+, the only condition is that +M+ does not have entries with denominators in +P+).
*/
{Specialize M in a prime ideal lying over the prime number p.}

    K := BaseRing(M);
    if Type(K) eq FldRat then
        L, quotmor := ResidueClassFieldGF(p);
    else
        P := RandomPrimeIdeal(K,p);
        L, quotmor := ResidueClassFieldGF(P);
    end if;
    return ChangeRing(M, quotmor);

end intrinsic;

//==============================================================================
intrinsic Entries(M::ModGrp) -> SetEnum
/*
    Intrinsic: Entries

    The set of entries of the action generators of a GModule.

    Declaration:
        :intrinsic Entries(M::ModGrp) -> SetEnum

    Parameters:
       M - A GModule

    Description:
        Returns the set of entries of the action generators of +M+.
*/
{}

    entries := {};
    for i:=1 to Ngens(M) do
        entries join:=Entries(ActionGenerator(M,i));
    end for;

end intrinsic;


//==============================================================================
intrinsic IsGModule(M::ModGrp) -> BoolElt
/*
    Intrinsic: IsGModule

    Checks if a GModule is indeed a G-module.

    Declaration:
        :intrinsic IsGModule(M::ModGrp) -> BoolElt

    Parameters:
       M - A GModule.

    Description:
        Checks if Representation(+M+) is a morphism using <IsMorphism>.
*/
{True iff M defines a module for Group(M).}

    return IsMorphism(Representation(M));

end intrinsic;


//==============================================================================
intrinsic IsFaithful(M::ModGrp) -> BoolElt
/*
    Intrinsic: IsFaithful

    Checks if a GModule is faithful.

    Declaration:
        :intrinsic IsFaithful(M::ModGrp) -> BoolElt

    Parameters:
       M - A GModule.

    Description:
        Checks if Representation(+M+) is injective using <IsInjective>.
*/
{}

    return IsFaithful(Representation(M));

end intrinsic;
