/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


/*
	Intrinsics for attaching the absolutely irreducible representations.
*/


//===================================================================
/*
    Namespace: Grp

    Additions to the category +Grp+.
*/
declare attributes Grp:
    Modules, //An integer indexed associative array containing the absolutely irreducible modules in each characteristic.
    ModuleNames, //An integer indexed associative array containing names of the modules in Modules.
    Representations, //Same as <Modules> but containing actual group morphisms.
    RepresentationNames, //Same as <ModuleNames> for Modules but for Representations.
    RegularModule, //The regular module
	CharactersDualityInvolution,
	RepresentationsDualityInvolution;

//============================================================================
intrinsic Representations(~G::Grp, p::RngIntElt : UseDB:=true, Check:=false)
/*
    Intrinsic: Representations

    Compute or load the absolutely irreducible representations of specified characteristic.

    Declaration:
        :intrinsic Representations(~G::Grp, p::RngIntElt : UseDB:=true, Check:=false)

    Parameters:
       G - A group.
       p - The characteristic.

    Options:
        UseDB - Possible choices are +true+ or +false+. If +true+, then try to load representations from data base.

    Description:
        Compute or load the absolutely irreducible representations of +G+ of characteristic +p+.

    History:
        * Sunday, September 15, 2013 16:25:19: Initial.
        * Thursday, September 19, 2013 11:17:20: Introduced characteristic.

*/
{}

    if assigned G`Representations and IsDefined(G`Representations, p) then
        return;
    end if;

    if not assigned G`Representations then
        G`Representations := AssociativeArray(Integers()); //indexed by characteristic
    end if;
    if not assigned G`Modules then
        G`Modules := AssociativeArray(Integers()); //indexed by characteristic
    end if;
    if not assigned G`RepresentationNames then
        G`RepresentationNames := AssociativeArray(Integers()); //indexed by characteristic
    end if;

    //first, load characters if p=0
    if p eq 0 then
        CharacterTable(~G : UseDB:=UseDB, Check:=Check);
    end if;

    //first try to find representations in the database
    foundindb := false;
    if UseDB and assigned G`DBDir and CHAMP_ExistsInDB(G`DBDir, "Representations/Representations_"*Sprint(p)) then
        reps := CHAMP_GetFromDB(G`DBDir, "Representations/Representations_"*Sprint(p));
        G`Representations[p] := [* *];
        for i:=1 to #reps do
            K := BaseRing(reps[i][1]);
            d := Ncols(reps[i][1]);
            f := hom< G->GeneralLinearGroup(d, K) | reps[i]>;
            if Check then
                if not IsMorphism(f) then
                    error "(Representations.) One representation does not define a morphism!";
                end if;
                if not IsIrreducible(f) then
                    error "(Representations.) One representation is not irreducible!";
                end if;
            end if;
            Append(~G`Representations[p], f);
        end for;

        RepresentationNames(~G,p);
        foundindb := true;

        //if character table is already assigned we have to remap the representations perhaps
        if p eq 0 and assigned G`CharacterTable then
            chars := [ Character(G`Representations[0][i]) : i in [1..#G`Representations[0]] ];
            sigma := [ Position(chars, G`CharacterTable[i]) : i in [1..#G`CharacterTable] ];
            if #sigma ne #chars then
                error "Something went wrong.";
            end if;
            if sigma ne [1..#sigma] then
                reps := G`Representations[0];
                G`Representations[0] := [ reps[sigma[i]] : i in [1..#G`Representations[0]] ];
            end if;
            //if the characters have names, we use them also for the representations
            if assigned G`CharacterNames then
                G`RepresentationNames[p] := G`CharacterNames;
            end if;
        end if;
    end if;

    if foundindb then
        return;
    end if;

    //if we could not find representations in the database we try to construct them.

    //splitting field by Brauer
    K := SplittingField(G,p);

    if IsSolvable(G) then
        G`Modules[p] := AbsolutelyIrreducibleModules(G,K);
    else
        G`Modules[p] := Constituents(RegularModule(G,K));
    end if;

    //try to reduce scalars in positive characteristic
    if p gt 0 then
        for i:=1 to #G`Modules[p] do
            L := MinimalField(G`Modules[p][i]);
            G`Modules[p][i] := ChangeRing(G`Modules[p][i], L);
        end for;
    end if;

    G`Representations[p] := [ Representation(G`Modules[p][i]) : i in [1..#G`Modules[p]] ];

    //match with character table
    if p eq 0 then
        CharacterTable(~G);
        repchars := [ Character(G`Representations[0][i]) : i in [1..#G`CharacterTable]];
        sigma := [ Position(repchars, G`CharacterTable[i]) : i in [1..#G`CharacterTable]];
        if sigma ne [1..#G`CharacterTable] then
            G`Representations[0] := [ G`Representations[0][sigma[i]] : i in [1..#G`CharacterTable] ];
            G`Modules[0] := [ G`Modules[0][sigma[i]] : i in [1..#G`CharacterTable] ];
        end if;
    end if;

    RepresentationNames(~G,p);

end intrinsic;

//============================================================================
intrinsic RepresentationNames(~G::Grp, p::RngIntElt)
{}

    if assigned G`RepresentationNames and IsDefined(G`RepresentationNames, p) then
        return;
    end if;

    Representations(~G,p);

    if not assigned G`RepresentationNames then
        G`RepresentationNames := AssociativeArray(Integers());
    end if;

    if not assigned G`ModuleNames then
        G`ModuleNames := AssociativeArray(Integers());
    end if;

    if assigned G`DBDir and CHAMP_ExistsInDB(G`DBDir, "Representations/RepresentationNames_"*Sprint(p)) then
        G`RepresentationNames[p] := CHAMP_GetFromDB(G`DBDir, "Representations/RepresentationNames_"*Sprint(p));
        G`ModuleNames[p] := G`RepresentationNames[p];
    end if;

end intrinsic;

//============================================================================
intrinsic RepresentationNames(~G::Grp)
{}

    RepresentationNames(~G,0);

end intrinsic;

//============================================================================
intrinsic ModuleNames(~G::Grp,p::RngIntElt)
{}

    RepresentationNames(~G,p);

end intrinsic;

//============================================================================
intrinsic Representations(~G::Grp : UseDB:=true, Check:=false)
{}

    Representations(~G,0: UseDB:=UseDB, Check:=Check);

end intrinsic;

//============================================================================
intrinsic Modules(~G::Grp, p::RngIntElt : UseDB:=true, Check:=false)
{}

    if assigned G`Modules and IsDefined(G`Modules, p) then
        return;
    end if;

    if not assigned G`Modules then
        G`Modules := AssociativeArray(Integers()); //indexed by characteristic
    end if;

    if not assigned G`ModuleNames then
        G`ModuleNames := AssociativeArray(Integers()); //indexed by characteristic
    end if;

    Representations(~G,p:UseDB:=UseDB, Check:=Check);

    if assigned G`Representations and IsDefined(G`Representations, p) then
        G`Modules[p] := [ GModule(G`Representations[p][i]) : i in [1..#G`Representations[p]] ];
        if IsDefined(G`RepresentationNames, p) then
            G`ModuleNames[p] := G`RepresentationNames[p];
        end if;
    end if;

end intrinsic;

//============================================================================
intrinsic Modules(~G::Grp : UseDB:=true)
{}

    Modules(~G,0:UseDB:=UseDB);

end intrinsic;

//============================================================================
intrinsic Modules(G::Grp,p::RngIntElt : UseDB:=true) -> SeqEnum
{}

    Modules(~G,p: UseDB:=UseDB);
    return G`Modules;

end intrinsic;

//============================================================================
intrinsic Modules(G::Grp : UseDB:=true) -> SeqEnum
{}

    Modules(G,0 : UseDB:=UseDB);

end intrinsic;

//============================================================================
intrinsic RegularModule(~G::Grp)
/*
    Intrinsic: RegularModule

    The regular module.

    Declaration:
        :intrinsic RegularModule(G::Grp, R::Rng) -> ModGrp

    Parameters:
       G - A group.
       R - A ring.

    Description:
        Returns the regular +RG+-module for a group +G+ and a ring +R+. To understand the action explicitly we have to enumerate the elements of +G+ and have to keep track of this. To this end, <NumberingMap> is computed and the regular module is computed with respect to the enumeration given by this set.

    History:
        * Wednesday, January 15, 2014 at 14:13:46: Using NumberingMap now.
        * Friday, September 20, 2013 00:52:13: Initial.

*/
{}

    if assigned G`RegularModule then
        return;
    end if;

    R:=Integers();
    NumberingMap(~G);
    opmats := [];
    if Type(G) eq GrpPC then
        N := #PCGenerators(G);
    else
        N := Ngens(G);
    end if;
    for k:=1 to N do
        A := ZeroMatrix(R, #G, #G);
        for i:=1 to #G do
            j := G`NumberingMap( G`InverseNumberingMap(i)*G.k );
            A[i][j] := 1;
        end for;
        Append(~opmats, A);
    end for;

    G`RegularModule := GModule(G, opmats);

end intrinsic;

//============================================================================
intrinsic RegularModule(G::Grp, R::Rng) -> ModGrp
{}

    RegularModule(~G);
    return ChangeRing(G`RegularModule, R);

end intrinsic;


//============================================================================
intrinsic DiagonalizeRepresentation(rho::HomGrp, num::RngIntElt) -> HomGrp
{Over a field of characteristic zero perfomes a base change and possibly field extension so that G.num acts diagonally on rho.}

    if Characteristic(BaseRing(Codomain(rho))) eq 0 then
        G:=Domain(rho);
        M := rho(G.num);
        K := MinimalCyclotomicField(SplittingField(CharacteristicPolynomial(M)));
        M := ChangeRing(M, K);
        J,T := JordanForm(M);
        N := Nrows(M);
        fj := hom<G -> GeneralLinearGroup(N, K) | [T*ChangeRing(rho(G.i),K)*T^-1 : i in [1..Ngens(G)]]>;
        return fj;
    else
        error "Not implemented yet";
    end if;

end intrinsic;

//============================================================================
intrinsic DiagonalizeRepresentations(~G::Grp, p::RngIntElt, num::RngIntElt)
/*
    Intrinsic: DiagonalizeRepresentations

    Declaration:
        :intrinsic DiagonalizeRepresentations(~G::Grp, p::RngIntElt, num::RngIntElt)
        :intrinsic DiagonalizeRepresentations(~G::Grp, num::RngIntElt)
        :intrinsic DiagonalizeRepresentations(~G::GrpMat)

    Parameters:
       G - A group.
       p - The characteristic.
       num - An integer.

    History:
        * Tuesday, April 1, 2014 at 17:18:02: Renamed from JordanizeRepresentations to DiagonalizeRepresentations.
        * Tuesday, October 01, 2013 22:13:29: Initial.

    To-do:
        * So far only implemented for characteristic zero (uses +MinimalCyclotomicField+).

*/
{}
    if not assigned G`Modules then
        G`Modules := AssociativeArray(Integers());
    end if;

    G`Modules[p] := [];

    for i:=1 to #G`Representations[p] do
		f := G`Representations[p][i];
        fj := DiagonalizeRepresentation(f, num);
		G`Representations[p][i] := fj;
        Append(~G`Modules[p], GModule(fj));
	end for;

end intrinsic;

//============================================================================
intrinsic DiagonalizeRepresentations(~G::Grp, num::RngIntElt)
/*
    History:
        * Tuesday, April 1, 2014 at 17:20:45: Initial.
*/
{}

    DiagonalizeRepresentations(~G, 0, num);

end intrinsic;

//============================================================================
intrinsic DiagonalizeRepresentations(~G::GrpMat)
/*

    History:
        * Tuesday, April 1, 2014 at 17:21:27: Initial.
*/
{}

    //detect if one generator of G is diagonal
    i:=1;
    while i le Ngens(G) do
        if IsDiagonal(Matrix(G.i)) then
            break;
        else
            i+:=1;
        end if;
    end while;

    if i gt Ngens(G) then
        i:=1;
    end if;

    DiagonalizeRepresentations(~G,i);

end intrinsic;

//============================================================================
intrinsic LiftRepresentationsToCommonBaseField(~G::Grp,p::RngIntElt)
/*
    Intrinsic: LiftRepresentationsToCommonBaseField

    Lifts the absolutely irreducible representations of characteristic +p+ to a common base field.

    Declaration:
        :intrinsic LiftRepresentationsToCommonBaseField(~G::Grp,p::RngIntElt)
        and
        :intrinsic LiftRepresentationsToCommonBaseField(~G::Grp)

    Parameters:
       G - A group.
       p - The characteristic.

    Description:
        Uses <CommonOverfield> to determine a field containing all entries of the matrices describing the absolutely irreducible modules of +G+ of characteristic +p+ as stored in <Grp.Representations> or <Grp.Modules>. In case the parameter +p+ is not provided, this intrinsic is called with +p=0+.

    History:
        * Monday, October 07, 2013 23:44:00: Initial.

*/
{}
    Modules(~G,p);
    Representations(~G,p);

    if not assigned G`Representations then
        G`Representations := AssociativeArray(Integers());
        G`Representations[p] := [];
    end if;

    K := CommonOverfield([* BaseRing(G`Modules[p][i]) : i in [1..#G`Modules[p]] *]);
    G`Modules[p] := [ChangeRing(G`Modules[p][i], K) : i in [1..#G`Representations[p]] ];
    G`Representations[p] := [ Representation(G`Modules[p][i]) : i in [1..#G`Representations[p]]];

end intrinsic;

//============================================================================
intrinsic LiftRepresentationsToCommonBaseField(~G::Grp)
/*
    History:
        Monday, October 07, 2013 23:46:03: Initial.
*/
{}

    LiftRepresentationsToCommonBaseField(~G,0);

end intrinsic;



//============================================================================
intrinsic SplittingField(G::Grp, p::RngIntElt) -> Fld
/*
    History:
        * Tuesday, April 1, 2014 at 17:37:41: Initial
*/
{A splitting field for G in characteristic p (this will not necessarily be a minimal one).}

    //Using Brauer theorems
    m := Exponent(G);
    if p eq 0 then
        K := CyclotomicField(m);
    elif not IsDivisibleBy(Order(G), p) then
        f := Order(p,m);
        K := GF(p^f);
    else
        C := CyclotomicField(m);
        O := RingOfIntegers(C);
        I := ideal<O|p>;
        fact := Factorization(I);
        K := ResidueClassFieldGF(fact[1][1]); //some maximal ideal over p
    end if;

    return K;

end intrinsic;







//============================================================================
procedure AddToIrreduciblesList(~Irr, M)

    isinlist:=false;
    for N in Irr do
        if IsIsomorphic(N,M) then
            isinlist:=true;
            break;
        end if;
    end for;
    if not isinlist then
        Append(~Irr, M);
    end if;

end procedure;

//============================================================================
intrinsic GaloisConjugate(M::ModGrp, sigma::Map) -> ModGrp
{}

    return GModule(Group(M), [ChangeRing(ActionGenerator(M,i), sigma) : i in [1..#ActionGenerators(M)]]);

end intrinsic;

//============================================================================
intrinsic FindIrreducibles(G::Grp, ~Irr::SeqEnum[ModGrp] : Limit:=0, Galois:={})
{}

    tensorlength := 0;

    while #Irr lt Limit or Limit eq 0 do
        tensorlength +:=1;
        newadded:=false;
        space := CartesianProduct([{1..#Irr} : i in [1..tensorlength]]);
        for tensor in space do
            Perfomed join:={tensor};
            M := Irr[tensor[1]];
            for i:=2 to tensorlength do
                M := TensorProduct(M, Irr[tensor[i]]);
            end for;
            consts := Constituents(M);
            for M in consts do
                AddToIrreduciblesList(~Irr, M);
                AddToIrreduciblesList(~Irr, Dual(M));
                for sigma in Galois do
                    AddToIrreduciblesList(~Irr, GaloisConjugate(M,sigma));
                end for;
            end for;
        end for;
        PrintAndDelete("Number of irreducibles: "*Sprint(#Irr));
    end while;


end intrinsic;

//============================================================================
intrinsic GrothendieckRing(~G::Grp : SaveToDB:=false)
{}

    if assigned G`GrothendieckRing then
        return;
    end if;

    if assigned G`DBDir then
        if CHAMP_ExistsInDB(G`DBDir, "Representations/GrothendieckRing") then
            G`GrothendieckRing := CHAMP_GetFromDB(G`DBDir, "Representations/GrothendieckRing");
            return;
        end if;
    end if;

    CharacterTable(~G);
    V := VectorSpace(Rationals(), #G`CharacterTable);
    structs := [];
    for i in [1..#G`CharacterTable] do
        for j in [i..#G`CharacterTable] do
            PrintPercentage(i, #G`CharacterTable);
            dec := V!Decomposition(G`CharacterTable, G`CharacterTable[i]*G`CharacterTable[j]);
            for k in Support(dec) do
                Append(~structs, <i,j,k,dec[k]>);
                if i ne j then
                    Append(~structs, <j,i,k,dec[k]>);
                end if;
            end for;
        end for;
    end for;
    G`GrothendieckRing := AssociativeAlgebra<Rationals(), Dimension(V)|structs : Check:=false>;

    if SaveToDB then
        str := "/*\n\tSplit characteristic zero Grothendieck ring\n\tDate: "*Sprint(Date())*"\n*/\n";
        str *:= "structs:="*Sprint(structs, "Magma")*";\n";
        str *:= "dim:="*Sprint(Dimension(V))*";\n";
        str *:= "\nreturn AssociativeAlgebra<Rationals(), dim|structs : Check:=false>;";
        Write(CHAMP_GetDBDir()*G`DBDir*"GrothendieckRing.m", str);
    end if;

end intrinsic;

//============================================================================
intrinsic CharacterField(G::Grp) -> Fld
{The characterf field of G.}

    CharacterTable(~G);
    Ks := [* CharacterField(G`CharacterTable[i]) : i in [1..#G`CharacterTable] *];
    K := Ks[1];
    for i:=2 to #Ks do
        K := Compositum(K, Ks[i]);
    end for;
    if Type(K) eq FldCyc and CyclotomicOrder(K) eq 1 then
        K := Rationals();
    end if;
    return K;

end intrinsic;

//============================================================================
intrinsic NumberOfpRegularClasses(G::Grp, p::RngIntElt) -> RngIntElt
{}

    Classes(~G);
    return #[ i : i in [1..#G`Classes] | IsDivisibleBy(Order(G`Classes[i][3]),p) eq false ]; //number of p-regular classes

end intrinsic;

//============================================================================
intrinsic SplitsOverRationals(G::Grp) -> BoolElt
{True iff G splits over the rational numbers.}

    K := CharacterField(G);
    if Type(K) ne FldRat then
        return false;
    end if;
    CharacterTable(~G);
    for i:=1 to #G`CharacterTable do
        if SchurIndex(G`CharacterTable[i]) ne 1 then
            return false;
        end if;
    end for;

    return true;

end intrinsic;

//============================================================================
intrinsic CharactersDualityInvolution(~G::Grp)
/*
    Intrinsic: CharactersDualityInvolution

    The involution defined by taking duals on the absolutely irreducible characters.

    Declaration:
        :intrinsic CharactersDualityInvolution(~G::Grp)

    Parameters:
       G - A group.

    Description:
        This assigns to +G`CharactersDualityInvolution+ the map on the indices of
        characters defined by taking duals.
*/
{}

	if assigned G`CharactersDualityInvolution then
		return;
	end if;

	CharacterTable(~G);
	G`CharactersDualityInvolution := [];
	for i:=1 to #G`CharacterTable do
		chi := G`CharacterTable[i];
		chidual := ComplexConjugate(chi);
		pos := Position(G`CharacterTable, chidual);
		if pos eq 0 then
			error "Something went wrong!";
		end if;
		Append(~G`CharactersDualityInvolution, pos);
	end for;

end intrinsic;

//============================================================================
intrinsic RepresentationsDualityInvolution(~G::Grp, p::RngIntElt)
/*
    Intrinsic: RepresentationsDualityInvolution

    The involution defined by taking duals on the absolutely irreducible characteristic p representations.

    Declaration:
        :intrinsic RepresentationsDualityInvolution(~G::Grp, p::RngIntElt)

    Parameters:
       G - A group.
       p - The characteristic.

    Description:
        This assigns to +G`RepresentationsDualityInvolution[p]+ the map on the indices of
        absolutely irreducible characteristic p representations stored in +G`Modules[p]+ defined
        by taking duals.

    History:
        * Monday, October 07, 2013 23:44:00: Initial.

*/
{}

	if not assigned G`RepresentationsDualityInvolution then
		G`RepresentationsDualityInvolution := AssociativeArray(Integers());
	end if;

	if IsDefined(G`RepresentationsDualityInvolution, p) then
		return;
	end if;

	if p eq 0 then
		CharactersDualityInvolution(~G);
		G`RepresentationsDualityInvolution[0] := G`CharactersDualityInvolution;
		return;
	else
    	Modules(~G,p);
    	G`RepresentationsDualityInvolution[p] := [];
    	cands := {1..#G`Modules[p]};
    	for i:=1 to #G`Modules[p] do
        	//have to find position of dual of G`Modules[p][i] first
        	Mdual := Dual(G`Modules[p][i]);
        	identified := false;
        	for j in cands do
            	K := CommonOverfield([* BaseRing(Mdual), BaseRing(G`Modules[p][j]) *]);
            	if IsIsomorphic(ChangeRing(Mdual,K), ChangeRing(G`Modules[p][j],K)) then
                	identifed := true;
                	cands diff:={j};
                  	Mdualpos := j;
                	break;
            	end if;
        	end for;
        	if not identifed then
    	    	error "Something went wrong";
        	end if;
        	Append(~G`RepresentationsDualityInvolution[p], Mdualpos);
        end for;
    end if;

end intrinsic;
