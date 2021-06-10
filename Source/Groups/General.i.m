/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for groups in general.
*/

//==============================================================================
intrinsic IspRegular(g::GrpElt, p::RngIntElt) -> BoolElt
/*
    Intrinsic: IspRegular

    Checks if an element is +p+-regular

    Declaration:
        :intrinsic IspRegular(g::GrpElt, p::RngIntElt) -> BoolElt

    Parameters:
    	g - an element of a group (of any type)
    	p - a prime number or zero

    Description:
    	Returns true iff +g+ is +p+-regular, i.e., the order of +g+ is coprime to +p+.
*/
{}

    return IsCoprime(Order(g), p);

end intrinsic;


//==============================================================================
intrinsic StrongGenerators(~G::Grp)
/*
    Intrinsic: StrongGenerators

    Strong generators of the group.

    Declaration:
        :intrinsic StrongGenerators(~G::Grp)

    Parameters:
       G - A group of type +GrpMat+ or +GrpPerm+.

    Description:
        Assigns strong generators to <Grp.StrongGenerators>.
*/
{}

    if assigned G`StrongGenerators then
        return;
    end if;

    if Type(G) eq GrpMat or Type(G) eq GrpPerm then
        G`StrongGenerators := IndexedSetToSequence(StrongGenerators(G));
        G`NumberOfStrongGenerators := #G`StrongGenerators;

        /*if [ G.i : i in [1..Ngens(G)] ] ne [ G`StrongGenerators[i] : i in [1..Ngens(G)]] then
            error "Strong generators do not refine generators!";
        end if;*/

    end if;

end intrinsic;

//==============================================================================
intrinsic StrongGeneratorWords(~G::Grp)
/*
    Intrinsic: StrongGeneratorWords

    Strong generators as words in the generators of the group.

    Declaration:
        :intrinsic StrongGeneratorWords(~G::Grp)

    Parameters:
       G - A group of type +GrpMat+ or +GrpPerm+.

    Description:
        Assigns presentations of the strong generators in the generators to <Grp.StrongGeneratorWords>.
*/
{}

    if assigned G`StrongGeneratorWords then
        return;
    end if;

    StrongGenerators(~G);

    G`StrongGeneratorWords := [ ElementToWord(G`StrongGenerators[i] : NoInverse:=true) : i in [1..G`NumberOfStrongGenerators]];

end intrinsic;


//==============================================================================
intrinsic InclusionMap(G::Grp, H::Grp) -> Map
/*
    Intrinsic: InclusionMap

    The inclusion map of a subgroup.

    Declaration:
        :intrinsic InclusionMap(G::Grp, H::Grp) -> Map

    Parameters:
       G - A group.
       H - A subgroup of +G+.

    Description:
        Returns the inclusion morphism from the subgroup +H+ of the group +G+. Internally, this intrinsic is already defined for types +GrpPC+ and +GrpGPC+.

    History:
        * Wednesday, January 15, 2014 at 13:12:33: Initial.

*/
{The map from the subgroup H of G to G.}

    return hom<H->G | [G!H.i : i in [1..Ngens(H)]]>;

end intrinsic;


//==============================================================================
/*
    Intrinsic: IsMorphism

    Checks if a map between groups is a morphism

    Declaration:
        :intrinsic IsMorphism(f::HomGrp : Method:="FPGroup") -> BoolElt
        :intrinsic IsMorphism(f::Map[Grp, AlgMat] : Method:="FPGroup") -> BoolElt

    Parameters:
       f - A morphism of groups

    Description:
        If +f+ is of type +HomGrp+ (note that when definining morphisms Magma does not check if this really defines a morphism, so +HomGrp+ does not yet mean that it's a group morphism, just that it's a map between groups constructed using +hom+), then check if +f+ is indeed a group morphism. This really checks all relations for an +FPGroup+ of the domain of +f+.

    History:
        * Thursday, February 20, 2014 at 12:41:25: Moved to GrpSystem.
        * Monday, September 23, 2013 17:48:51: Initial.

*/
intrinsic IsMorphism(f::HomGrp : Method:="FPGroup") -> BoolElt
{True iff f is a group morphism.}

    G := Domain(f);

    if f(Identity(G)) ne Identity(Codomain(f)) then
        return false;
    end if;

    for i:=1 to Ngens(G) do
        for j:=1 to Ngens(G) do
            if f(G.i*G.j) ne f(G.i)*f(G.j) then
                return false;
            end if;
        end for;
    end for;

    if Method eq "FPGroup" then
        FPGroup(~G);
        for rel in Relations(G`FPGroup) do
            lhs := Eltseq(LHS(rel));
            rhs := Eltseq(RHS(rel));

            lhseval := Identity(Codomain(f));
            for i:=1 to #lhs do
                lhseval *:= f(G.lhs[i]);
            end for;

            rhseval := Identity(Codomain(f));
            for i:=1 to #rhs do
                rhseval *:= f(G.rhs[i]);
            end for;

            if lhseval ne rhseval then
                return false;
            end if;
        end for;
    end if;

    return true;

end intrinsic;

//==============================================================================
/*
    Sometimes representations of a group have the matrix ring as codomain and are therefore not of type HomGrp. This intrinsic catches these cases.
*/
intrinsic IsMorphism(f::Map[Grp, AlgMat] : Method:="FPGroup") -> BoolElt
{}
    C := GL(Degree(Codomain(f)), BaseRing(Codomain(f)));
    G := Domain(f);
    newf := hom<G -> C | [ f(G.i) : i in [1..Ngens(G)]]>;
    return IsMorphism(newf : Method:=Method);

end intrinsic;

//==============================================================================
/*
    Intrinsic: IsInjective

    Checks if a group morphism is injective

    Declaration:
        :intrinsic IsInjective(f::HomGrp) -> BoolElt
        :intrinsic IsInjective(f::Map[Grp, AlgMat]) -> BoolElt

    Parameters:
       f - A morphism of groups

    Description:
        Checks if +f+ is injective.

    To-do:
        This is really ugly and has to be rewritten.

    History:
        * Thursday, February 20, 2014 at 12:41:25: Moved to GrpSystem.
        * Tuesday, September 24, 2013 12:23:08: Initial.

*/
intrinsic IsInjective(f::HomGrp) -> BoolElt
{}

    G := Domain(f);
    for g in Set(G) do
        if g eq Identity(G) then
            continue;
        end if;
        if f(g) eq Identity(Codomain(f)) then
            return false;
        end if;
    end for;

    return true;

end intrinsic;

//==============================================================================
intrinsic IsInjective(f::Map[Grp, AlgMat]) -> BoolElt
{}

    C := GL(Degree(Codomain(f)), BaseRing(Codomain(f)));
    G := Domain(f);
    newf := hom<G -> C | [ f(G.i) : i in [1..Ngens(G)]]>;
    return IsInjective(newf);

end intrinsic;

//==============================================================================
/*
    Intrinsic: IsFaithful

    Checks if a group morphism is injective

    Declaration:
        :intrinsic IsFaithful(f::HomGrp) -> BoolElt
        :intrinsic IsFaithful(f::Map[Grp, AlgMat]) -> BoolElt

    Parameters:
       f - A morphism of groups

    Description:
        Same as <IsInjective> for group morphisms.

*/
intrinsic IsFaithful(f::HomGrp) -> BoolElt
{}

    return IsInjective(f);

end intrinsic;

//==============================================================================
intrinsic IsFaithful(f::Map[Grp, AlgMat]) -> BoolElt
{}

    C := GL(Degree(Codomain(f)), BaseRing(Codomain(f)));
    G := Domain(f);
    newf := hom<G -> C | [ f(G.i) : i in [1..Ngens(G)]]>;
    return IsInjective(newf);

end intrinsic;

//==============================================================================
/*
    Intrinsic: Image

    The image of a group morphism.

    Declaration:
        :intrinsic Image(f::HomGrp) -> Grp
        :intrinsic Image(f::Map[Grp, AlgMat]) -> Grp

    Parameters:
       f - A morphism of groups

    Description:
        Returns the image of a group morphism +f+ as a subgroup.

    History:
        * Wednesday, June 4, 2014 at 00:54:41: Removed. Causes a conflict somewhere.
        * Tuesday, September 24, 2013 12:32:56: Initial.

*/
/*intrinsic Image(f::HomGrp) -> Grp
{}

    G := Domain(f);
    C := Codomain(f);
    return sub<C|[f(g) : g in Set(G)]>;

end intrinsic;
*/
//=============================================================================
intrinsic Image(f::Map[Grp, AlgMat]) -> Grp
{}

    C := GL(Degree(Codomain(f)), BaseRing(Codomain(f)));
    G := Domain(f);
    newf := hom<G -> C | [ f(G.i) : i in [1..Ngens(G)]]>;
    return Image(newf);

end intrinsic;
