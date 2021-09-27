/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Handlers for explicit types of reflection groups.
*/

intrinsic ComplexReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt : Realization:="CHEVIE") -> GrpMat
{}

  label := "G"*Sprint(m)*"_"*Sprint(p)*"_"*Sprint(n)*"_"*Realization;

  if not CHAMP_ExistsInDB("ReflectionGroups/"*label, "GrpMat") then
    error "Group is not in database. Use export/import scripts in Champ-DB/Tools to add group from CHEVIE.";
  end if;

  G := CHAMP_GetFromDB("ReflectionGroups/"*label, "GrpMat");
  G`IsReflectionGroup := true;
  G`IsPolynomial := true;

  return G;

end intrinsic;

intrinsic ComplexReflectionGroup(n::RngIntElt : Realization:="CHEVIE") -> GrpMat
{}

  label := "G"*Sprint(n)*"_"*Realization;

  G := CHAMP_GetFromDB("ReflectionGroups/"*label, "GrpMat");
  G`IsReflectionGroup := true;
  G`IsPolynomial := true;

  return G;

end intrinsic;

//============================================================================
intrinsic ExceptionalComplexReflectionGroup(n::RngIntElt : Realization:="CHEVIE") -> GrpMat
{
	Returns an explicit realization of the Shephard–Todd exceptional complex reflection group G_n. Currently, only realization CHEVIE is supported which returns the transposed of ComplexReflectionGroup(n) as in *gap3-jm5*. In the database more information about these groups are stored.
}

    if Realization eq "CHEVIE" then
    	G := CHAMP_GetFromDB("ReflectionGroups/G"*Sprint(n)*"_CHEVIE", "GrpMat"); //automatically sets DBDir
        DualGroup(~G);
        Degrees(~G);
        G`DualGroup`Degrees := G`Degrees;
    else
    	error "Realization not known";
    end if;

    G`IsReflectionGroup := true;
    G`DualGroup`IsReflectionGroup := true;
    G`IsPolynomial := true;
    G`DualGroup`IsPolynomial := true;

    return G;

end intrinsic;

//============================================================================
intrinsic SymmetricReflectionGroup(n::RngIntElt : UseDB:=true) -> GrpMat, HomGrp
{
	The natural irreducible reflection representation of S_(n+1) over K, together with an isomorphism from S_(n+1). K can be any field of characteristic coprime to the degree n. If K is not provided, then K is taken to be the rational numbers. In this case and for n between 2 and 22 and if UseDB is true, we import the data from *gap3-jm5*.
}

	G, phi := SymmetricReflectionGroup(n, Rationals() : UseDB:=UseDB);

    return G, phi;

end intrinsic;

//============================================================================
intrinsic SymmetricReflectionGroup(n::RngIntElt, K::Fld : UseDB:=true) -> GrpMat, HomGrp
{}

    if n lt 2 then
        error "n >= 2 required.";
    end if;

    if Characteristic(K) ne 0 and IsDivisibleBy(n,Characteristic(K)) then
        error "The characteristic has to be coprime to the degree n.";
    end if;

    S := SymmetricGroup(n);

	if n le 22 and Type(K) eq FldRat and UseDB then
		G := CHAMP_GetFromDB("ReflectionGroups/S"*Sprint(n)*"_CHEVIE", "GrpMat");
		DualGroup(~G);
    else
    	mats := [ ];
   	 	for i:=1 to n-1 do
       		//construct the reflection corresponding to the transposition s_{i,i+1}.
       	 	M := DiagonalMatrix(K, [1:i in [1..n-1]]);
        	M[i][i] := -1;
        	if i gt 1 then
            	M[i][i-1] := 1;
        	end if;
        	if i lt n-1 then
            	M[i][i+1] := 1;
        	end if;
        	//Append(~mats, Transpose(M)); //have to transpose matrices!
        	Append(~mats, M); //have to transpose matrices!
    	end for;
    	G := MatrixGroup<n-1, K | mats>;
    end if;


    //to define the morphism we need the symmetric group with generators the transpositions
    T := sub<S|[ S!(i,i+1) : i in [1..n-1]]>;
    phi := hom<T -> G | [G.i : i in [1..n-1]]>;
    res, psi := IsIsomorphic(S,T);

    G`IsReflectionGroup := true;
    DualGroup(~G);
    G`IsPolynomial := true; //Kemper-Malle!
    G`DualGroup`IsPolynomial := true;

    return G,hom<S->G | [phi(S.i) : i in [1..Ngens(S)]]>;

end intrinsic;

//============================================================================
intrinsic CyclicReflectionGroup(m::RngIntElt, K::Fld) -> GrpMat, Map
{
	The natural reflection representation of the cyclic group of order m over the field K of characteristic coprime to m. Here, natural means that it is the one acting by the primitive m-th root of unity fixed by Magma. If K is not provided, it is taken to be the cyclotomic field of order m. We automatically attach the irreducible representations.
}

    if m lt 2 then
        error "m >= 2 required.";
    end if;

    if Characteristic(K) ne 0 and IsDivisibleBy(m,Characteristic(K)) then
        error "The characteristic has to be coprime to m.";
    end if;

    zeta := RootOfUnity(m, K);
    if Parent(zeta) ne K then
        error "The base field K has to have a primitive m-th root of unity.";
    end if;

    G := MatrixGroup<1,K | [zeta]>;
    C := AbelianGroup([m]);

    //representations
    p := Characteristic(K);
    G`Representations := AssociativeArray(Integers());
    G`Representations[p] := [];
    for r:=1 to m do
        g := Matrix(K,1,1,[zeta^r]);
        phi := hom<G->GL(1,K) | [g] >;
        Append(~G`Representations[p], phi);
    end for;

    return G, hom<C->G | [zeta]>;

end intrinsic;


//============================================================================
intrinsic CyclicReflectionGroup(m::RngIntElt) -> GrpMat, Map
{}

    return CyclicReflectionGroup(m, CyclotomicField(m));

end intrinsic;


//============================================================================
intrinsic TypeBReflectionGroup(n::RngIntElt, K::Fld : UseDB:=true) -> GrpMat
{
	A reflection representation for the Weyl group of type B_n. If n is between 2 and 22, K is the rational field and UseDB is true, then G(2,1,n) is loaded as in gap3-jm5.
}

	if n le 12 and Type(K) eq FldRat and UseDB then
		G := CHAMP_GetFromDB("ReflectionGroups/B"*Sprint(n)*"_CHEVIE", "GrpMat");
		DualGroup(~G);
        G`IsReflectionGroup := true;
    	DualGroup(~G);
    	G`IsPolynomial := true;
   	 	G`DualGroup`IsPolynomial := true;
		return G;
	else
		return ImprimitiveReflectionGroup(2,1,n,K);
	end if;

end intrinsic;

//============================================================================
intrinsic TypeBReflectionGroup(n::RngIntElt : UseDB:=true) -> GrpMat
{}

	return TypeBReflectionGroup(n, Rationals() : UseDB:=UseDB);

end intrinsic;



//============================================================================
intrinsic TypeDReflectionGroup(n::RngIntElt, K::Fld : UseDB:=true) -> GrpMat
{
	A reflection representation for the Weyl group of type +Dn+. If n is between 2 and 22, +K+ is the rational field and +UseDB+ is true, then +G(2,2,n)+ is loaded as in gap3-jm5.
}

	if n le 12 and Type(K) eq FldRat and UseDB then
		G := CHAMP_GetFromDB("ReflectionGroups/D"*Sprint(n)*"_CHEVIE", "GrpMat");
		DualGroup(~G);
        G`IsReflectionGroup := true;
    	DualGroup(~G);
    	G`IsPolynomial := true;
   	 	G`DualGroup`IsPolynomial := true;
		return G;
	else
		return ImprimitiveReflectionGroup(2,2,n,K);
	end if;

end intrinsic;

//============================================================================
intrinsic TypeDReflectionGroup(n::RngIntElt : UseDB:=true) -> GrpMat
{}

	return TypeDReflectionGroup(n, Rationals() : UseDB:=UseDB);

end intrinsic;




//============================================================================
intrinsic ImprimitiveReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt, K::Fld) -> GrpMat
{
	The imprimitive Shephard-Todd reflection group G(m,p,n)
}

    permmats := [];
    for i:=1 to n-1 do
        s := IdentityMatrix(K,n);
        s[i][i] := 0;
        s[i+1][i+1] := 0;
        s[i][i+1] := 1;
        s[i+1][i] := 1;
        Append(~permmats,s);
    end for;

    t := IdentityMatrix(K,n);
    res, zeta := HasRootOfUnity(m,K);
    if not res then
        error "Field does not have a primitive "*Sprint(m)*"-th root of unity.";
    end if;
    t[1][1] := zeta;

    G := MatrixGroup<n, K|[ t^p, t^-1*permmats[1]*t ] cat permmats >;

    return G;

    /*C := CyclicGroup(m);
    T := DirectProduct([C : i in [1..n]]);
    gamma := hom<T->C | [ (C.1)^(m div p) : i in [1..n]]>;
    A := Kernel(gamma);

    return A;*/

end intrinsic;

//============================================================================
intrinsic ImprimitiveReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt : Sparse:=true, Realization:="Thesis") -> GrpMat
{}
    if Realization eq "Thesis" then
    	K := CyclotomicField(m:Sparse:=true);
        return ImprimitiveReflectionGroup(m,p,n,K);
    else
    	error "Chosen realization not known.";
    end if;

end intrinsic;

//============================================================================
intrinsic ImprimitiveReflectionGroupCenterOrder(m::RngIntElt, p::RngIntElt, n::RngIntElt) -> RngIntElt
{
	The order of the center of the imprimitive reflection group G(m,p,n).
}

    return GreatestCommonDivisor([ k*m : k in [1..n-1]] cat [(n*m) div p]);

end intrinsic;

//============================================================================
intrinsic DihedralReflectionGroup(m::RngIntElt, K::Fld) -> GrpMat
{
	The dihedral reflection group of order 2*m over K. We automaticlly attach the irreducible representations.
}

	if not m ge 3 then
		error "m >= 3 needed.";
	end if;
	if Characteristic(K) eq 2 then
		error "Characteristic 2 not allowed";
	end if;
	if not Characteristic(K) eq 0 and IsCoprime(Characteristic(K), m) then
		error "Characteristic has to be coprime to m.";
	end if;
	res, zeta := HasRootOfUnity(m,K);
	if not res then
		error "Field has to have primitive m-th root of unity";
	end if;
	x := Matrix(K, 2, 2, [0,1,1,0]);
	y := Matrix(K, 2, 2, [0, zeta^-1, zeta, 0]);
	G := MatrixGroup<2,K | [x,y]>;

	p := Characteristic(K);
	Classes(~G);
	G`Representations := AssociativeArray(Integers());
	G`Representations[p] := [];
	G`RepresentationNames := AssociativeArray(Integers());
	G`RepresentationNames[p] := [];
	if p eq 0 then
		G`CharacterTable := [];
		G`CharacterNames := [];
	end if;
	if IsOdd(m) then
		for i:=1 to Integers()!((m-1)/2) do
			rhox := x;
			rhoy := Matrix(K, 2, 2, [0, zeta^-i, zeta^i, 0]);
			rho := hom<G->GL(2,K) | [rhox,rhoy]>;
			Append(~G`Representations[p], rho);
			if p eq 0 then
				Append(~G`CharacterTable, Character(rho));
				Append(~G`CharacterNames, "rho"*Sprint(i));
			end if;
			Append(~G`RepresentationNames[p], "rho"*Sprint(i));
		end for;
		triv := hom<G->GL(1,K) | [ [1], [1] ]>;
		Append(~G`Representations[0], triv);
		if p eq 0 then
			Append(~G`CharacterTable, Character(triv));
			Append(~G`CharacterNames, "1");
		end if;
		Append(~G`RepresentationNames[p], "1");
		eps := hom<G->GL(1,K) | [ [-1], [-1] ]>;
		Append(~G`Representations[p], eps);
		if p eq 0 then
			Append(~G`CharacterTable, Character(eps));
			Append(~G`CharacterNames, "eps");
		end if;
		Append(~G`RepresentationNames[p], "eps");
	else
		for i:=1 to Integers()!((m-2)/2) do
			rhox := x;
			rhoy := Matrix(K, 2, 2, [0, zeta^-i, zeta^i, 0]);
			rho := hom<G->GL(2,K) | [rhox,rhoy]>;
			Append(~G`Representations[p], rho);
			if Characteristic(K) eq 0 then
				Append(~G`CharacterTable, Character(rho));
				Append(~G`CharacterNames, "rho"*Sprint(i));
			end if;
			Append(~G`RepresentationNames[p], "rho"*Sprint(i));
		end for;
		triv := hom<G->GL(1,K) | [ [1], [1] ]>;
		Append(~G`Representations[p], triv);
		if p eq 0 then
			Append(~G`CharacterTable, Character(triv));
			Append(~G`CharacterNames, "1");
		end if;
		Append(~G`RepresentationNames[p], "1");

		eps := hom<G->GL(1,K) | [ [-1], [-1] ]>;
		Append(~G`Representations[p], eps);
		if p eq 0 then
			Append(~G`CharacterTable, Character(eps));
			Append(~G`CharacterNames, "eps");
		end if;
		Append(~G`RepresentationNames[p], "eps");

		epsbar := hom<G->GL(1,K) | [ [1], [-1] ]>;
		Append(~G`Representations[p], epsbar);
		if p eq 0 then
			Append(~G`CharacterTable, Character(epsbar));
			Append(~G`CharacterNames, "eps1");
		end if;
		Append(~G`RepresentationNames[p], "eps1");

		gamma := hom<G->GL(1,K) | [ [-1], [1] ]>;
		Append(~G`Representations[p], gamma);
		if p eq 0 then
			Append(~G`CharacterTable, Character(gamma));
			Append(~G`CharacterNames, "eps2");
		end if;
		Append(~G`RepresentationNames[p], "eps2");
	end if;

	return G;

end intrinsic;


//============================================================================
intrinsic DihedralReflectionGroup(m::RngIntElt) -> GrpMat, GrpFPCox
{}

	G := DihedralReflectionGroup(m, CyclotomicField(m));

	//G`aFunction := [ 1 : i in [1..#G`CharacterTable] ];
	pos := Position(G`CharacterNames, "1");
	//G`aFunction[pos] := 0;
	pos := Position(G`CharacterNames, "eps");
	//G`aFunction[pos] := m;

	N:=Matrix(Integers(),2,2,[1,m,m,1]);
	W := CoxeterGroup(GrpFPCox,N);

	P1 := sub<G|G.1>;
	CharacterTable(~P1);
	if P1`CharacterTable[1](P1.1) ne -1 then
		X := P1`CharacterTable;
		P1`CharacterTable[1] := X[2];
		P1`CharacterTable[2] := X[1];
	end if;
	P1`CharacterNames := ["1", "psi1"];
	P2 := sub<G|G.2>;
	CharacterTable(~P2);
	if P2`CharacterTable[1](P2.1) ne -1 then
		X := P2`CharacterTable;
		P2`CharacterTable[1] := X[2];
		P2`CharacterTable[2] := X[1];
	end if;
	P2`CharacterNames := ["1", "psi2"];

	G`ParabolicSubgroups := [sub<G|>, P1, P2, G];

	return G,W;

end intrinsic;

//============================================================================
intrinsic RL(n::RngIntElt, q::RngIntElt) -> GrpMat
{
	The subgroup of GL(n,q) generated by all the reflections.
}

    G := GL(n,q);
    Reflections(~G);
    H := sub<G|G`Reflections>;
    H`IsReflectionGroup := true;
    return H;

end intrinsic;

//============================================================================
intrinsic RU(n::RngIntElt, q::RngIntElt) -> GrpMat
{The subgroup of GU(n,q) generated by all the reflections.}

    G := GU(n,q);
    Reflections(~G);
    H := sub<G|G`Reflections>;
    H`IsReflectionGroup := true;
    return H;

end intrinsic;


//============================================================================
intrinsic FindModularDiagonalizableReflectionGroups(Orders::SetEnum : AvoidPrimes:={2}) -> SeqEnum
{
	Finds all modular diagonalizable reflection groups of given order.
}

    groups := [];
    for o in Orders do
        print o;
        for i:=1 to NumberOfSmallGroups(o) do
            PrintPercentage(i, NumberOfSmallGroups(o));
            G := PermutationGroup(SmallGroup(o, i));
            for p in SequenceToSet(PrimeFactors(o)) diff AvoidPrimes do
                Modules(~G, p);
                for j:=1 to #G`Modules[p] do
                    if not IsFaithful(G`Modules[p][j]) then
                        continue;
                    end if;
                    H := Image(Representation(G`Modules[p][j]));
                    if IsDiagonalizableReflectionGroup(H) then
                        print "Found one: "*Sprint(<o,i>);
                        Append(~groups, <o,i,p,H>);
                    end if;
                end for;
            end for;
        end for;
    end for;

    return groups;

end intrinsic;

//============================================================================
intrinsic ReflectionRepresentations(G::Grp, p::RngIntElt : DiagonalizableOnly:=false) -> SeqEnum
{}

    Representations(~G,p);

    reps := [];
    for j:=1 to #G`Representations[p] do
        if not IsFaithful(G`Representations[p][j]) then
            continue;
        end if;
        H := Image(G`Representations[p][j]);
        if (not DiagonalizableOnly and IsReflectionGroup(H)) or (DiagonalizableOnly and IsDiagonalizableReflectionGroup(H)) then
            Append(~reps, H);
        end if;
    end for;

    return reps;

end intrinsic;
