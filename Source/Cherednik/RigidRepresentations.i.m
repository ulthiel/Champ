/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Rigid representations
*/

//========================================================================
intrinsic IsRigid(chi::AlgChtrElt, c::Map) -> BoolElt
{Checks if chi is a (weakly!) c-rigid character.}

	G := Group(Parent(chi));
	ReflectionLibrary(~G);
	R := Codomain(c);
	for i,j in [1..Dimension(G)] do
		sum := 0;
		for s in G`ReflectionLibraryFlat do
			sum +:= c(s`ReflectionClass)*R!(CherednikCoefficient(i,j,s)*chi(s`Element));
		end for;
		if not sum eq 0 then
			return false;
		end if;
	end for;

	return true;

end intrinsic;

//========================================================================
intrinsic IsRigid(rho::Map, c::Map) -> BoolElt
{Checks if rho is a c-rigid representation.}

	G := Domain(rho);
	L := BaseRing(G);
	K := BaseRing(Codomain(rho));
	M := CommonOverfield(K,L);
	R := ChangeRing(Codomain(c), M);
	if Characteristic(K) eq 0 then
		if not IsRigid(Character(rho),c) then
			return false;
		end if;
	end if;

	for i,j in [1..Dimension(G)] do
		sum := 0;
		for s in G`ReflectionLibraryFlat do
			sum +:= R!c(s`ReflectionClass)*(R!CherednikCoefficient(i,j,s))*ChangeRing(rho(s`Element), R);
		end for;
		if not sum eq 0 then
			return false;
		end if;
	end for;

	return true;

end intrinsic;

//========================================================================
intrinsic RigidRepresentations(G::GrpMat, c::Map) -> SetEnum
{The c-rigid representations of G.}

	K := BaseRing(G);
	p := Characteristic(K);
	Representations(~G, p);

	return { i : i in [1..#G`Representations[p]] | IsRigid(G`Representations[p][i], c) };

end intrinsic;

//========================================================================
intrinsic RigidCharacters(G::GrpMat, c::Map) -> SetEnum
{The (weakly!) c-rigid characters of G.}

	CharacterTable(~G);

	return { i : i in [1..#G`CharacterTable] | IsRigid(G`CharacterTable[i], c) };

end intrinsic;

//========================================================================
intrinsic RigidityVariety(rho::HomGrp  : Type:="EG") -> RngMPol
{The variety on which rho becomes c-rigid.}

	G := Domain(rho);
	c := CherednikParameter(G : Type:="EG", Rational:=false);
	R := Codomain(c);
	F := {};
	for i,j in [1..Dimension(G)] do
		sum := 0;
		for s in G`ReflectionLibraryFlat do
			sum +:= c(s`ReflectionClass)*(R!CherednikCoefficient(i,j,s))*ChangeRing(rho(s`Element), R);
		end for;
		F join:=Entries(sum);
	end for;

	return GroebnerBasis(ideal<R|F>);

end intrinsic;
