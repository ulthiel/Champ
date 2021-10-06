/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Compute the center of the rational Cherednik algebra at t=0.

	Joint work with Cédric Bonnafé (Montpellier).
*/


declare attributes AlgChe:
	CenterGenerators, 	//generators of the center
	CenterGeneratorsOfDegreeZero,	//generators of the center of degree zero
	CalogeroMoserFamilies, 	//Calogero-Moser families. This is an associative array indexed by the exceptional hyperplanes.
	CenterGeneratorsPoissonBrackets,	//{z_i,z_j} for the center generators z_i,z_j such that deg(z_i) + deg(z_j) = 0.
	CuspidalCalogeroMoserFamilies;	//the cuspidal Calogero-Moser families


//============================================================================
intrinsic IsRegular(W::GrpMat, v::ModTupFldElt) -> BoolElt
{v is regular for the action of W if v is not contained in any of the reflectin hyperplanes of W.}

	ReflectionLibrary(~W);
	for s in W`ReflectionLibraryFlat do
		if v in s`Hyperplane then
			return false;
		end if;
	end for;
	return true;

end intrinsic;


//============================================================================
intrinsic RegularVector(W::GrpMat) -> ModTupFldElt
{
Returns a regular vector, i.e., one not fixed by any reflection.
}

	//We try to find a regular vector randomly.

	ReflectionLibrary(~W);
	V := VectorSpace(W);
	n := Dimension(V);
	d := 0;
	while true do
		d +:= 1;
		coeffs := CartesianProduct([ [0..d] cat [-j : j in [0..d]] : i in [1..n]]);
		for c in coeffs do
			v := V![c[i] : i in [1..n]];
			if IsRegular(W,v) then
				return v;
			end if;
		end for;
	end while;

end intrinsic;
//============================================================================
intrinsic Truncation(h::AlgCheElt) -> RngMPolElt
{
The truncation is defined in [BR15]. This is simply the coefficient of the trivial group element of h.
}

	return Coefficient(h`Element, Identity(h`Parent`Group));

end intrinsic;

//============================================================================
intrinsic TruncationInverse(H::AlgChe, f::RngMPolElt : yreg:=0, Threads:=1) -> AlgCheElt
{
	f is an element of k[V^* + V] (note ordering V^*+V). Compute the inverse of f in the center of H under the truncation operator.
}

	P := Parent(f);
	R := BaseRing(H`GroupAlgebra);
	W := H`Group;
	bideg := Bidegree(f);
	delta := Minimum({bideg[1], bideg[2]});

	phi := hom<P->R | [R.i : i in [1..Ngens(R)]]>;

	z := phi(f)*H`GroupAlgebra!Identity(W);

	if yreg eq 0 then
		yreg := RegularVector(W); //a "nice" random choice
	end if;

	yregR := &+[yreg[i]*H`yxAlgebra.i : i in [1..Dimension(W)]];

	for r:=1 to delta+1 do
		print "Step "*Sprint(r)*" of "*Sprint(delta+1);
		newz := phi(f)*H`GroupAlgebra!Identity(W); //this will be the new z

		countw := 0;
		for w in W do	//construct the coefficient of w of newz
			countw +:= 1;
			PrintPercentage(countw, #W);
			if w eq Identity(W) then
				continue;
			end if;
			rhs := Zero(R);

			refls := {k : k in [1..#W`ReflectionLibraryFlat] | w*W`ReflectionLibraryFlat[k]`Element^-1 in Support(z)};

			count := 0;
			for k in refls do
				s := W`ReflectionLibraryFlat[k];
				a := Coefficient(z, w*s`Element^-1);
				b := &+[ yreg[i]*Commutator(H, i, a, k) : i in Support(yreg) ];
				rhs +:= b;
			end for;

			if rhs eq 0 then
				continue;
			else
				c := rhs div (SymplecticDoublingAction(yregR, w) - yregR);
			end if;
			newz +:= c*H`GroupAlgebra!w;
		end for;
		z := newz;
		//return z;
		//PrintPercentage(r, delta+1);
		print "";
	end for;

	zH := Zero(H);
	zH`Element := z;

	if Truncation(zH) ne phi(f) then
		error "Something wrong with truncation inverse.";
	end if;

	return zH;

end intrinsic;


//============================================================================
intrinsic CenterGenerators(~H::AlgChe : UseDB:=true, SaveToDB:=false)
{Add alls center generators.}

	W := H`Group;
	SymplecticDoublingFundamentalInvariants(~W);

	for i:=1 to #W`SymplecticDoublingFundamentalInvariants do
		if assigned H`CenterGenerators and H`CenterGenerators[i] ne Zero(H) then
			continue;
		else
			print "Deforming center generator "*Sprint(i)*" of "*Sprint(#W`SymplecticDoublingFundamentalInvariants);
			CenterGenerator(~H,i : UseDB:=UseDB, SaveToDB:=SaveToDB);
		end if;
	end for;

end intrinsic;

intrinsic CenterGenerators(H::AlgChe : UseDB:=true, SaveToDB:=false) -> SeqEnum
{}

	CenterGenerators(~H);
	return H`CenterGenerators;

end intrinsic;

//============================================================================
intrinsic CenterGenerator(~H::AlgChe, i::RngIntElt : UseDB:=true, SaveToDB:=false)
{Compute and add i-th generator of the center}

	W := H`Group;

	if not assigned H`CenterGenerators then
		SymplecticDoublingFundamentalInvariants(~W);
		H`CenterGenerators := [* Zero(H) : i in [1..#W`SymplecticDoublingFundamentalInvariants] *];
	end if;

	if H`CenterGenerators[i] eq Zero(H) then
		gotfromdb := false;
		if UseDB and not SaveToDB and assigned H`DBDir then
			if CHAMP_ExistsInDB(H`DBDir, "CenterGenerators/"*Sprint(i)) then
				z := CHAMP_GetFromDB(H`DBDir, "CenterGenerators/"*Sprint(i));
				H`CenterGenerators[i] := H!z;
				delete z;	//perhaps not necessary
				print "Found center generator in DB.";
				gotfromdb := true;
			end if;
		end if;
		if not gotfromdb then
			H`CenterGenerators[i] := TruncationInverse(H, W`SymplecticDoublingFundamentalInvariants[i]);
		end if;
	end if;

	if SaveToDB then
		if not assigned H`DBDir then
			error "Algebra has no database directory assigned (needed for saving). Cannot save.";
		end if;
		gen, type := IsGeneric(H);
		CHAMP_SaveToDB(Rprint(H`CenterGenerators[i]), H`DBDir, "CenterGenerators/"*Sprint(i));
	end if;

end intrinsic;

//============================================================================
intrinsic CenterGeneratorsOfDegreeZero(~H::AlgChe : UseDB:=true, SaveToDB:=false)
{Adds all center generators of degree zero.}

	W := H`Group;
	SymplecticDoublingFundamentalInvariants(~W);
	deg0gens := [i : i in [1..#W`SymplecticDoublingFundamentalInvariants] | Bidegree(W`SymplecticDoublingFundamentalInvariants[i])[1] eq Bidegree(W`SymplecticDoublingFundamentalInvariants[i])[2]];

	count := 0;
	for i in deg0gens do
		count +:= 1;
		print "Deforming center generator "*Sprint(count)*" of "*Sprint(#deg0gens);
		CenterGenerator(~H, i : SaveToDB:=SaveToDB, UseDB:=UseDB);
	end for;

	if not assigned H`CenterGeneratorsOfDegreeZero then
		H`CenterGeneratorsOfDegreeZero := [* Zero(H) : i in [1..#deg0gens] *];
	end if;

	for i:=1 to #deg0gens do
		j := deg0gens[i];
		H`CenterGeneratorsOfDegreeZero[i] := H`CenterGenerators[j];
	end for;

end intrinsic;

intrinsic CenterGeneratorsOfDegreeZero(H::AlgChe : UseDB:=true, SaveToDB:=false)
-> List
{}

	CenterGeneratorsOfDegreeZero(~H);
	return H`CenterGeneratorsOfDegreeZero;

end intrinsic;

//============================================================================
intrinsic CenterGeneratorsPoissonBracket(~H::AlgChe, i::RngIntElt, j::RngIntElt : UseDB:=true, SaveToDB:=false)
{The Poisson bracket of z_i and z_j.}

	W := H`Group;

	if not assigned H`CenterGeneratorsPoissonBrackets then
		SymplecticDoublingFundamentalInvariants(~W);
		H`CenterGeneratorsPoissonBrackets := AssociativeArray({<k,l> : k,l in [1..#W`SymplecticDoublingFundamentalInvariants]});
	end if;

	if not IsDefined(H`CenterGeneratorsPoissonBrackets, <i,j>) then
		if IsDefined(H`CenterGeneratorsPoissonBrackets, <j,i>) then
			H`CenterGeneratorsPoissonBrackets[<i,j>] := -H`CenterGeneratorsPoissonBrackets[<j,i>];
		elif i eq j then
			H`CenterGeneratorsPoissonBrackets[<i,j>] := Zero(H);	//because Poisson brackets are anti-symmetric
		else
			gotfromdb := false;
			if UseDB and not SaveToDB and assigned W`DBDir then
				gen, type := IsGeneric(H);
				if gen and assigned type and CHAMP_ExistsInDB(W`DBDir, "Cherednik/Generic/"*type*"CenterGeneratorsPoissonBrackets/"*Sprint(i)*"_"*Sprint(j)) then
					z := CHAMP_GetFromDB(W`DBDir, "Cherednik/Generic/"*type*"CenterGeneratorsPoissonBrackets/"*Sprint(i)*"_"*Sprint(j));
					H`CenterGeneratorsPoissonBrackets[<i,j>] := H!z;
					delete z;	//perhaps not necessary
					print "Found Poisson bracket in DB.";
					gotfromdb := true;
				end if;
			end if;
			if not gotfromdb then
				//print "HERE";
				CenterGenerator(~H, i);
				CenterGenerator(~H, j);
				H`CenterGeneratorsPoissonBrackets[<i,j>] := PoissonBracket(H`CenterGenerators[i], H`CenterGenerators[j]);
			end if;
		end if;
	end if;

	if SaveToDB then
		if not assigned W`DBDir then
			error "Group has no database directory assigned (needed for saving).";
		end if;
		gen, type := IsGeneric(H);
		if gen and assigned type then
			CHAMP_SaveToDB(Rprint(H`CenterGeneratorsPoissonBrackets[<i,j>]), W`DBDir, "CherednikGeneric"*type*"CenterGeneratorsPoissonBracket"*Sprint(i)*"_"*Sprint(j));
		else
			error "Saving not yet possible for this type of base ring.";
		end if;
	end if;

end intrinsic;


//============================================================================
intrinsic ActionOnVermaModule(chi::AlgChtrElt, z::AlgCheElt) -> RngElt
{The scalar by which the center element z acts on the baby Verma module attached to chi. We assume that z is bihomogeneous here!}

	//We assume that z is bihomogeneous here!
	H := z`Parent;
	act := Zero(H`BaseRing);
	zW := GroupAlgebraPart(z);
	R := H`BaseRing;
	for w in Support(zW) do
		act +:= Coefficient(zW,w)*(R!chi(w));
	end for;
	act *:= 1/(R!chi(1));

	return act;

end intrinsic;

intrinsic CenterGeneratorsAvailable(H::AlgChe) -> BoolElt
{}

	W := H`Group;
	SymplecticDoublingFundamentalInvariants(~W);
	inv := W`SymplecticDoublingFundamentalInvariants;
	for i:=1 to #inv do
		if not CHAMP_ExistsInDB(H`DBDir, "CenterGenerators/"*Sprint(i)) then
			return false;
		end if;
	end for;
	return true;

end intrinsic;
