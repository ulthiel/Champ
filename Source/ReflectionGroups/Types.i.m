//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Constructors for (complex) reflection groups.
//
//==============================================================================

//==============================================================================
// G1
//==============================================================================
intrinsic SymmetricReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat, HomGrp
{
	The standard irreducible reflection representation of the symmetric group, together with an isomorphism to Magma's symmetric group (mapping the generating reflections to the transpositions). The model is the same as in CHEVIE (which is also the same as the Lehrer-Taylor model in Magma).
}
	name := "S"*Sprint(n)*"_"*Model;

	G := CHAMP_GetFromDB("ReflectionGroups", name);

	G`IsReflectionGroup := true;

	//Generators of the above will correspond to the transpositions.
	//We also want to return an isomorphism to SymmetricGroup.
	//To define the morphism we need the symmetric group with generators the transpositions
	S := SymmetricGroup(n);
	T := sub<S|[ S!(i,i+1) : i in [1..n-1]]>;
	phi := hom<T -> G | [G.i : i in [1..n-1]]>;
	res, psi := IsIsomorphic(S,T);

	return G, hom<S->G | [phi(S.i) : i in [1..Ngens(S)]]>;;

end intrinsic;

intrinsic TypeAReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{
	Alias for SymmetricReflectionGroup.
}
	return SymmetricReflectionGroup(n : Model:=Model);

end intrinsic;



//==============================================================================
// G2
//==============================================================================
intrinsic CyclicReflectionGroup(m::RngIntElt) -> GrpMat, Map
{
	The standard reflection representation of the cyclic group of order m. Here, standard means that it is the one acting by the primitive m-th root of unity fixed by Magma. We automatically attach the irreducible representations with a standard labeling.
}

	if m lt 2 then
		error "m >= 2 required.";
	end if;

	K := CyclotomicField(m);
	zeta := K.1;

	G := MatrixGroup<1,K | [zeta]>;
	C := AbelianGroup([m]);

	//representations
	G`Representations := [];
	for r:=1 to m do
		g := Matrix(K,1,1,[zeta^r]);
		phi := hom<G->GL(1,K) | [g] >;
		Append(~G`Representations, phi);
	end for;

	G`IsReflectionGroup := true;

	return G, hom<C->G | [zeta]>;

end intrinsic;


//==============================================================================
// G4..G37
//==============================================================================
intrinsic ExceptionalComplexReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{
	Returns an explicit model of the Shephard–Todd exceptional complex reflection group G_n. The default model is "CHEVIE". The Model "LT" (Lehrer-Taylor) is the one used in Magma.
}

	name := "G"*Sprint(n)*"_"*Model;

	if Model eq "CHEVIE" then
		G := CHAMP_GetFromDB("ReflectionGroups", name);
	elif Model eq "LT" then
		G := ShephardTodd(n);
		G`DBName := name;
	else
		error "Unknown model";
	end if;

	G`IsReflectionGroup := true;

	return G;

end intrinsic;


//==============================================================================
// Weyl groups of type B and D
//==============================================================================

//============================================================================
intrinsic TypeBReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{
	Returns an explicit model of the Weyl group of type Bn as a reflection group. The default model is "CHEVIE". The Model "LT" (Lehrer-Taylor) is the one used in Magma.
}

	name := "B"*Sprint(n)*"_"*Model;

	if Model eq "CHEVIE" then
		G := CHAMP_GetFromDB("ReflectionGroups", name);
	elif Model eq "LT" then
		G := ShephardTodd(2,1,n);
		G`DBName := name;
	else
		error "Unknown model";
	end if;

	G`IsReflectionGroup := true;

	return G;

end intrinsic;

//============================================================================
intrinsic TypeDReflectionGroup(n::RngIntElt : Model := "CHEVIE") -> GrpMat
{
	Returns an explicit model of the Weyl group of type Dn as a reflection group. The default model is "CHEVIE". The Model "LT" (Lehrer-Taylor) is the one used in Magma.
}

	name := "D"*Sprint(n)*"_"*Model;

	if Model eq "CHEVIE" then
		G := CHAMP_GetFromDB("ReflectionGroups", name);
	elif Model eq "LT" then
		G := ShephardTodd(2,2,n);
		G`DBName := name;
	else
		error "Unknown model";
	end if;

	G`IsReflectionGroup := true;

	return G;

end intrinsic;


//==============================================================================
// Dihedral group
//==============================================================================

//============================================================================
intrinsic DihedralReflectionGroup(m::RngIntElt) -> GrpMat
{
	The standard reflection representation of the dihedral reflection group of order 2*m. We automaticlly attach the irreducible representations.
}

	if not m ge 3 then
		error "m >= 3 needed.";
	end if;

	K := CyclotomicField(m);
	zeta := K.1;

	x := Matrix(K, 2, 2, [0,1,1,0]);
	y := Matrix(K, 2, 2, [0, zeta^-1, zeta, 0]);
	G := MatrixGroup<2,K | [x,y]>;

	Classes(~G);
	G`Representations := [];
	G`RepresentationNames := [];
	G`CharacterTable := [];
	G`CharacterNames := [];

	if IsOdd(m) then
		for i:=1 to Integers()!((m-1)/2) do
			rhox := x;
			rhoy := Matrix(K, 2, 2, [0, zeta^-i, zeta^i, 0]);
			rho := hom<G->GL(2,K) | [rhox,rhoy]>;
			Append(~G`Representations, rho);
			Append(~G`CharacterTable, Character(rho));
			Append(~G`CharacterNames, "rho"*Sprint(i));
			Append(~G`RepresentationNames, "rho"*Sprint(i));
		end for;
		triv := hom<G->GL(1,K) | [ [1], [1] ]>;
		Append(~G`Representations, triv);
		Append(~G`CharacterTable, Character(triv));
		Append(~G`CharacterNames, "1");
		Append(~G`RepresentationNames, "1");
		eps := hom<G->GL(1,K) | [ [-1], [-1] ]>;
		Append(~G`Representations, eps);
		Append(~G`CharacterTable, Character(eps));
		Append(~G`CharacterNames, "eps");
		Append(~G`RepresentationNames, "eps");
	else
		for i:=1 to Integers()!((m-2)/2) do
			rhox := x;
			rhoy := Matrix(K, 2, 2, [0, zeta^-i, zeta^i, 0]);
			rho := hom<G->GL(2,K) | [rhox,rhoy]>;
			Append(~G`Representations, rho);
			Append(~G`CharacterTable, Character(rho));
			Append(~G`CharacterNames, "rho"*Sprint(i));
			Append(~G`RepresentationNames, "rho"*Sprint(i));
		end for;
		triv := hom<G->GL(1,K) | [ [1], [1] ]>;
		Append(~G`Representations, triv);
		Append(~G`CharacterTable, Character(triv));
		Append(~G`CharacterNames, "1");
		Append(~G`RepresentationNames, "1");

		eps := hom<G->GL(1,K) | [ [-1], [-1] ]>;
		Append(~G`Representations, eps);
		Append(~G`CharacterTable, Character(eps));
		Append(~G`CharacterNames, "eps");
		Append(~G`RepresentationNames, "eps");

		epsbar := hom<G->GL(1,K) | [ [1], [-1] ]>;
		Append(~G`Representations, epsbar);
		Append(~G`CharacterTable, Character(epsbar));
		Append(~G`CharacterNames, "eps1");
		Append(~G`RepresentationNames, "eps1");

		gamma := hom<G->GL(1,K) | [ [-1], [1] ]>;
		Append(~G`Representations, gamma);
		Append(~G`CharacterTable, Character(gamma));
		Append(~G`CharacterNames, "eps2");
		Append(~G`RepresentationNames, "eps2");
	end if;

	//parabolic subgroups
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

	return G;

end intrinsic;
