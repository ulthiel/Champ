/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions to provide certain internal intrinsics for all types of groups.

    Some Magma basic Magma intrinsics for groups, like +Center+ for example, do not work with
    finitely presented groups, even if they are finite. When working with a quaternion group for
    example represented as a finitely presented group we would like to directly apply +Center+
    to this group and get the result as a subgroup of the same type. This transport between
    different types of groups is handled here.
*/



//==============================================================================
intrinsic PermutationGroup(~G::Grp : DoDegreeReduction:=true)
/*
    Intrinsic: PermutationGroup

    A permutation group isomorphic to the group.

    Declaration:
        :intrinsic PermutationGroup(~G::Grp : DoDegreeReduction:=true)

    Parameters:
       G - A group of an arbitrary type.

    Options:
        DoDegreeReduction - Possible choices are +true+ and +false+. If set to +true+, then an additional degree reduction is performed. Default is +true+.

    Description:
        Sets <Grp.PermutationGroup> to a permutation group isomorphic to +G+, for any type of +G+. Furthermore, <Grp.ToPermutationGroup> and <Grp.FromPermutationGroup> are set to the isomorphism to and from the permutation group. If +DoDegreeReduction+ is set to +true+ (default), then an additional degree reduction (the internal one of Magma) is perfomed and the isomorphism are adapted.

    History:
        * Wednesday, January 15, 2014 at 11:56:10: Cleaned up and fixed PC group issue. Removed non-procedural version.
*/
{}

    if assigned G`PermutationGroup then
        return;
    end if;

    if Type(G) eq GrpPerm then
        G`PermutationGroup := G;
        G`ToPermutationGroup := IdentityHomomorphism(G);
        G`FromPermutationGroup := IdentityHomomorphism(G);
    elif Type(G) eq GrpFP or Type(G) eq GrpAb then
        G`PermutationGroup, G`ToPermutationGroup := PermutationGroup(G);
        G`FromPermutationGroup := Inverse(G`ToPermutationGroup);
    elif Type(G) eq GrpMat or Type(G) eq GrpAuto then
        G`ToPermutationGroup, G`PermutationGroup := PermutationRepresentation(G);
        G`FromPermutationGroup := Inverse(G`ToPermutationGroup);
    elif Type(G) eq GrpPC then
        FPGroup(~G);
        PermutationGroup(~G`FPGroup);
        G`PermutationGroup := G`FPGroup`PermutationGroup;
        G`ToPermutationGroup := hom<G -> G`PermutationGroup | [ G`FPGroup`ToPermutationGroup(G`ToFPGroup(G.i)) : i in [1..#PCGenerators(G)]]>; //jeez!
        G`FromPermutationGroup := Inverse(G`ToPermutationGroup);
    end if;

    if DoDegreeReduction then
        S,f := DegreeReduction(G`PermutationGroup);
        if Degree(S) lt Degree(G`PermutationGroup) then
            G`PermutationGroup := S;
            if Type(G) ne GrpPC then
                N := Ngens(G);
            else
                N := #PCGenerators(G);
            end if;
            G`ToPermutationGroup := hom<G->G`PermutationGroup | [ (G`ToPermutationGroup*f)(G.i) : i in [1..N]]>;
            G`FromPermutationGroup := Inverse(G`ToPermutationGroup);
        end if;
    end if;

end intrinsic;

//==============================================================================
intrinsic FPGroup(~G::Grp)
/*
    Intrinsic: FPGroup

    An FP group isomorphic to the group.

    Declaration:
        :intrinsic FPGroup(~G::Grp)

    Parameters:
       G - A group of an arbitrary type.

    Description:
        Sets <Grp.FPGroup> to an FP group isomorphic to +G+, for any type of +G+. Furthermore, <Grp.ToFPGroup> and <Grp.FromFPGroup> are set to the isomorphism to and from the FP group.
*/
{}

    if assigned G`FPGroup then
        return;
    end if;

    if Type(G) eq GrpFP then
        G`FPGroup := G;
        G`ToFPGroup := IdentityHomomorphism(G);
        G`FromFPGroup:= IdentityHomomorphism(G);
    elif Type(G) eq GrpMat or Type(G) eq GrpPerm or Type(G) eq GrpPC or Type(G) eq GrpAuto or Type(G) eq GrpAb then
        G`FPGroup, G`FromFPGroup := FPGroup(G);
        G`ToFPGroup := Inverse(G`FromFPGroup);
    end if;

end intrinsic;


//==============================================================================
intrinsic RWSGroup(~G::Grp)
/*
    Intrinsic: RWSGroup

    An RWS group isomorphic to the group.

    Declaration:
        :intrinsic RWSGroup(~G::Grp)

    Parameters:
       G - A group of an arbitrary type.

    Description:
        Sets <Grp.RWSGroup> to an RWS group isomorphic to +G+, for any type of +G+. Furthermore, <Grp.ToRWSGroup> and <Grp.FromRWSGroup> are set to the isomorphism to and from the RWS group.

    History:
        * Wednesday, January 15, 2014 at 12:11:02: Fixed some issues. Removed non-procedural version.
*/
{}
    if assigned G`RWSGroup then
        return;
    end if;

    if Type(G) eq GrpRWS then
        G`RWSGroup := G;
        G`ToRWSGroup := IdentityHomomorphism(G);
        G`FromRWSGroup:= IdentityHomomorphism(G);
    elif Type(G) eq GrpFP then
        G`RWSGroup := RWSGroup(G);
        G`ToRWSGroup := hom<G->G`RWSGroup | [G`RWSGroup.i : i in [1..Ngens(G)]]>;
        G`FromRWSGroup := hom<G`RWSGroup->G | [G.i : i in [1..Ngens(G)]]>;
    elif Type(G) eq GrpPerm or Type(G) eq GrpMat or Type(G) eq GrpPC then
        FPGroup(~G);
        RWSGroup(~G`FPGroup);
        G`RWSGroup := G`FPGroup`RWSGroup;
        //we have to factor through a morphism from the FP group; Magma won't accept a morphism with codomain an RWS group in general!
        if Type(G) ne GrpPC then
            N := Ngens(G);
        else
            N := #PCGenerators(G);
        end if;
        G`ToRWSGroup := G`ToFPGroup*hom<G`FPGroup->G`RWSGroup | [G`RWSGroup.i : i in [1..N]]>; //we have to define it in this way!
        G`FromRWSGroup := Inverse(G`ToRWSGroup);
    end if;

end intrinsic;

//==============================================================================
intrinsic PushforwardAttribute(~f::Map[Grp,Grp], attr::MonStgElt)
/*
    Intrinsic: PushforwardAttribute

    Pushes forward an attribute along a group morphism.

    Declaration:
        :intrinsic PushforwardAttribute(~f::Map[Grp,Grp], attr::MonStgElt)

    Parameters:
       f - A group morphism.
       attr - A string describing an attribute of a group.

    Description:
        If +f+ is an isomorphism of groups, push forward the attribute +attr+ from the domain of +f+ to the codomain of +f+. So far, +attr+ can be "Set", "Center", or "Classes".
*/
{}

    G := Domain(f);
    H := Codomain(f);

    if attr eq "Set" then
        H`Set := {f(g) : g in G};
    elif attr eq "Center" then
        H`Center := sub<H | [f(G`Center.i) : i in [1..Ngens(G`Center)]]>;
    elif attr eq "Classes" then
        H`Classes := [ <G`Classes[i][1], G`Classes[i][2], f(G`Classes[i][3])> : i in [1..#G`Classes]];
    end if;

end intrinsic;

//==============================================================================
intrinsic Set(~G::Grp : Method:="PermutationGroup")
/*
    Intrinsic: Set

    The underlying set of the group.

    Declaration:
        :intrinsic Set(~G::Grp : Method:="PermutationGroup")

    Parameters:
       G - A group of an arbitrary type.

    Options:
        Method - The method to be used for computing the set for group types not supported so far. Possible choices are "PermutationGroup" and "RWSGroup". Default is "PermutationGroup".

    Description:
        Sets <Grp.Set> to the underlying set of +G+, for any type of +G+.

    History:
        * Wednesday, January 15, 2014 at 12:17:52: Removed non-procedural versions.
*/
{}

    if assigned G`Set then
        return;
    end if;

    if Type(G) eq GrpPerm or Type(G) eq GrpMat or Type(G) eq GrpPC or Type(G) eq GrpRWS then
        G`Set := Set(G);
    elif Type(G) eq GrpFP or Type(G) eq GrpAuto then
        if Method eq "PermutationGroup" then
            PermutationGroup(~G);
            Set(~G`PermutationGroup);
            PushforwardAttribute(~G`FromPermutationGroup, "Set");
        elif Method eq "RWSGroup" then
            RWSGroup(~G);
            PushforwardAttribute(~G`FromRWSGroup, "Set");
        end if;
    end if;

end intrinsic;

//==============================================================================
intrinsic NumberingMap(~G::Grp)
/*
    Intrinsic: NumberingMap

    A numbering of the group elements.

    Declaration:
        :intrinsic NumberingMap(~G::Grp)

    Parameters:
       G - A group.

    Description:
        Sets <Grp.NumberingMap> to a bijective map from +G+ to +{1,...,|G|}+. Uses the internal function for types +GrpAb+, +GrpMat+, +GrpPC+, +GrpPerm+. The BGBS should not be changed afterwards as this changes the map.

        If +G+ is of type +GrpFP+, then the numbering map is not of type +Map+ but it is a user program. The reason for this is that a general group element does not look a priori like one in G`Set, only after applying relations. This is why we first have to map an element of G to its permutation group and then apply the numbering map of the permutation group to get a useful function. Hence, in this case, its a user program.

    History:
        * Wednesday, January 15, 2014 at 18:19:05: Changed numbering map to user program for type +GrpFP+.

*/
{}
    if Type(G) eq GrpAb or Type(G) eq GrpMat or Type(G) eq GrpPC or Type(G) eq GrpPerm or Type(G) eq GrpPermCox then
        //StrongGenerators(~G);   //would be done automatically but now it's fixed and assigned
        G`NumberingMap := NumberingMap(G);
        G`InverseNumberingMap := Inverse(G`NumberingMap);
        Set(~G);
    elif Type(G) eq GrpFP then
        PermutationGroup(~G);
        NumberingMap(~G`PermutationGroup);
        Set(~G);
        //the problem is that a general group element does not look a priori like one in G`Set, only after applying relations. This is why we first have to map an element of G to its permutation group and then apply the numbering map of the permutation group to get a useful function. Hence, in this case, its a user program.
        G`InverseNumberingMap := map<{1..#G`Set}->G`Set | [<G`PermutationGroup`NumberingMap(G`ToPermutationGroup(x)), x> : x in G`Set] >;
        G`NumberingMap := func<x|G`PermutationGroup`NumberingMap(G`ToPermutationGroup(x))>;
    end if;

end intrinsic;

//==============================================================================
intrinsic Center(~G::Grp)
/*
    Intrinsic: Center

    The center of a group.

    Declaration:
        :intrinsic Center(~G::Grp : Method:="PermutationGroup")

    Parameters:
       G - A group.

    Description:
        Sets <Grp.Center> to the center of +G+ which is computed using an isomorphic permutation group  for non-supported group types. It is constructed as a subgroup of +G+. Moreover, <Grp.CenterEmbedding> is the inclusion morphism from the center of +G+ into +G+ and <Grp.CenterGenerators> are the images of the generators of the center under this embedding.

    History:
        * Wednesday, January 15, 2014 at 12:56:03: Included CenterEmbedding and CenterGenerators.
        * Wednesday, January 15, 2014 at 12:20:24: Removed non-procedural version.
*/
{}

    if assigned G`Center then
        return;
    end if;

    if Type(G) eq GrpMat or Type(G) eq GrpPerm or Type(G) eq GrpAb or Type(G) eq GrpPC then
        G`Center := Center(G);
        G`CenterEmbedding := InclusionMap(G, G`Center);
    elif Type(G) eq GrpFP then
        PermutationGroup(~G);
        Center(~G`PermutationGroup);
        G`Center := sub<G | [ G`FromPermutationGroup(G`PermutationGroup`Center.i) : i in [1..Ngens(G`PermutationGroup`Center)] ]>;
        G`CenterEmbedding := hom<G`Center -> G | [ G`FromPermutationGroup(G`PermutationGroup`Center.i) : i in [1..Ngens(G`PermutationGroup`Center)] ]>;
    end if;

    G`CenterGenerators := {@ G`CenterEmbedding(G`Center.i) : i in [1..Ngens(G`Center)] @};

end intrinsic;


//==============================================================================
intrinsic CommutatorSubgroup(~G::Grp)
/*
    Intrinsic: CommutatorSubgroup

    The commutator subgroup of a group.

    Declaration:
        :intrinsic CommutatorSubgroup(~G::Grp)

    Parameters:
       G - A group.

    Description:
        Sets <Grp.CommutatorSubgroup> to the commutator subgroup of +G+. It is constructed as a subgroup of +G+. Moreover, <Grp.CommutatorSubgroupEmbedding> is the inclusion morphism from the center of +G+ into +G+ and <Grp.CommutatorSubgroupGenerators> are the images of the generators of the center under this embedding.
*/
{}

    if assigned G`CommutatorSubgroup then
        return;
    end if;

    G`CommutatorSubgroup := CommutatorSubgroup(G);
    G`CommutatorSubgroupEmbedding := InclusionMap(G, G`CommutatorSubgroup);
    G`CommutatorSubgroupGenerators := {@ G`CommutatorSubgroupEmbedding(G`CommutatorSubgroup.i) : i in [1..Ngens(G`CommutatorSubgroup)] @};

end intrinsic;


//=============================================================================
/*
    Namespace: Grp

    Additions to the category +Grp+.
*/
declare attributes Grp:
    PermutationGroup,
    ToPermutationGroup,
    FromPermutationGroup,
    RWSGroup,
    ToRWSGroup,
    FromRWSGroup,
    FPGroup,
    ToFPGroup,
    FromFPGroup,

    Set,
    NumberingMap,
    InverseNumberingMap,

    Center,
    CenterEmbedding,
    CenterGenerators,

    CommutatorSubgroup,
    CommutatorSubgroupEmbedding,
    CommutatorSubgroupGenerators,

    StrongGenerators,
    NumberOfStrongGenerators,
    StrongGeneratorWords,
    WordEmbedding,

    NumberOfGenerators,
    Generators

    ;

/*
    Attribute: PermutationGroup

    A permutation group isomorphic to the group.
*/

/*
    Attribute: ToPermutationGroup

    The isomorphism from the group to the permutation group.
*/

/*
    Attribute: FromPermutationGroup

    The isomorphism from the permutation group to the group.
*/

/*
    Attribute: RWSGroup

    An RWS group isomorphic to the group.
*/

/*
    Attribute: ToRWSGroup

    The isomorphism from the group to the RWS group.
*/


/*
    Attribute: FromRWSGroup

    The isomorphism from the RWS group to the group.
*/

/*
    Attribute: FPGroup

    An FP group isomorphic to the group.
*/

/*
    Attribute: ToFPGroup

    The isomorphism from the group to the FP group.
*/

/*
    Attribute: FromFPGroup

    The isomorphism from the FP group to the group.
*/

/*
    Attribute: Set

    The underlying set of the group.
*/

/*
    Attribute: NumberingMap

    A numbering of the elements of the group. Usually a map but a user program for type +GrpFP+.
*/

/*
    Attribute: InverseNumberingMap

    Inverse of <Grp.NumberingMap>.

/*
    Attribute: Center

    The center of the group.
*/

/*
    Attribute: CenterEmbedding

    The embedding of the center of the group into the group.
*/

/*
    Attribute: CenterGenerators

    Generators of the center of the group.
*/

/*
    Attribute: CommutatorSubgroup

    The commutator subgroup of the group.
*/

/*
    Attribute: CommutatorSubgroupEmbedding

    The embedding of the commutator of the group into the group.
*/

/*
    Attribute: CommutatorSubgroupGenerators

    Generators of the commutator subgroup of the group.
*/

/*
    Attribute: StrongGenerators

    Strong generators of the group
*/

/*
    Attribute: NumberOfStrongGenerators

    The number of strong generators of the group.
*/

/*
    Attribute: StrongGeneratorWords

    The strong generators as words in the generators.
*/
