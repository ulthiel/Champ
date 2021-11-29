/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
  Some tools for Calogero-Moser families
*/

//============================================================================
intrinsic GenericCalogeroMoserFamilies(W::GrpMat) -> SetEnum
{}

	c := CherednikParameter(W);
	return CalogeroMoserFamilies(W,c);

end intrinsic;

//============================================================================
/*
	Just a helper function for nice printing of families.
*/
function FamiliesString(W, F)

	str := "";
	for i:=1 to #F do
		str *:= Sprint(F[i]);
		if i lt #F then
			str *:= ", ";
		end if;
	end for;
	return str;

end function;

//============================================================================
intrinsic CalogeroMoserFamiliesTry(W::GrpMat, c::Map) -> SetEnum
{Tries to determine the CM families by simple tricks (using Euler element and supersingularity).}

	eu := {@ f[1] : f in EulerFamilies(W,c) @};

	//sort euler families
	eusorted := {@ @};
	worked := {};
	for i:=1 to #W`CharacterTable do
		if i in worked then
			continue;
		end if;
		for f in eu do
			if i in f then
				eusorted join:={@f@};
				worked join:={j : j in f};
			end if;
		end for;
	end for;
	eu := eusorted;
	print "The Euler families are:";
	IndentPush();
	print FamiliesString(W,eu);
	IndentPop();

	//singleton euler families are cm families
	cm := {@ f : f in eu | #f eq 1 @};
	print "\nSingleton Euler families are CM families, so the following are already CM families:";
	IndentPush();
	print FamiliesString(W,cm);
	IndentPop();
	Diff(~eu, cm);

	if IsEmpty(eu) then
		print "\nSucessfully determined the CM families. They are:";
		IndentPush();
		print FamiliesString(W,cm);
		IndentPop();
		return cm;
	end if;

	//supersingularity
	SupersingularRepresentations(~W);
	supsing := SequenceToIndexedSet(W`SupersingularRepresentations);
	print "\nThe supersingular characters are:";
	IndentPush();
	print FamiliesString(W,supsing);
	IndentPop();

	//good euler families
	goodeuler := {@ @};
	for f in eu do
		if #f eq 2 then
			if f[1] in supsing or f[2] in supsing then
				goodeuler join:={@f@};
			end if;
		elif #f eq 3 then
			if f[1] in supsing and f[2] in supsing and f[3] in supsing then
				goodeuler join:={@f@};
			end if;
		end if;
	end for;
	cm join:=goodeuler;
	Diff(~eu, goodeuler);
	print "\nThe following Euler families are CM families due to supersingularity:";
	IndentPush();
	print FamiliesString(W,goodeuler);
	IndentPop();

	if IsEmpty(eu) then
		print "\nSucessfully determined the CM families. They are:";
		IndentPush();
		print FamiliesString(W,cm);
		IndentPop();
		return cm;
	end if;

	//rigids
	rigids := {};

	return cm;

end intrinsic;

//============================================================================
intrinsic Residue(lambda::SeqEnum) -> SeqEnum
/*
    k-Residue of the partition X.
*/
{The residue of the partition lambda}

    R<x> := LaurentSeriesRing(Integers());
    res := Zero(R);

    for i:=1 to #lambda do
        for j:=1 to lambda[i] do
            res +:= x^(j-i);
        end for;
    end for;

    return res;

end intrinsic;


//============================================================================
intrinsic CalogeroMoserFamilies(m::RngIntElt, n::RngIntElt, H::SeqEnum[RngIntElt], kappa:RngIntElt) -> SeqEnum
{
	Calogero-Moser families for G(m,1,n) at H,kappa as defined by Martino.
}

    R<x> := LaurentSeriesRing(Integers());
    a := [0 : i in [1..m]];
    for i:=1 to m do
        a[i] := &+[H[j] : j in [1..i]];
    end for;
    //print a;
    residues := [];
    parts := Multipartitions(n,m);
    for lambda in parts do
        res := Zero(R);
        for i:=1 to 2 do
            res +:= x^a[i]*Evaluate(Residue(lambda[i]), x^(-kappa));
        end for;
        Append(~residues,res);
    end for;

    families := {};
    resset := SequenceToSet(residues);
    for res in resset do
        resset diff:={res};
        families join:= {{ parts[i] : i in [1..#parts] | residues[i] eq res}};
    end for;

    return families;

end intrinsic;

//============================================================================
intrinsic CalogeroMoserFamilies(~H::AlgChe : SaveToDB:=false, UseDB:=true)
{
	The Calogero-Moser families of H, determined via central characters evaluated at degree-0 center generators. The result is an assocaitive array indexed by the equations for the varieties of the exceptional locus, and by 1 for the generic situation.
}

	if assigned H`CalogeroMoserFamilies then
		return;
	end if;

	W := H`Group;
	R := H`BaseRing;

	if UseDB and assigned W`DBDir and CHAMP_ExistsInDB(W`DBDir, "Cherednik/Gordon") then
		gordon := CHAMP_GetFromDB(W`DBDir, "Cherednik/Gordon");
		try
			H`CalogeroMoserFamilies := AssociativeArray(SequenceToSet(gordon`BlGen));
			for h in Keys(H`CalogeroMoserFamilies) do
				H`CalogeroMoserFamilies[h] := gordon`Data[{h}]`CMFamilies;
			end for;
			return;
		catch e
			;
		end try;
	end if;

	CenterGeneratorsOfDegreeZero(~H);
	CharacterTable(~H`Group);

	P := PolynomialRing(Rationals(), Rank(R));
	AssignNames(~P, Names(R));
	H`CalogeroMoserFamilies := AssociativeArray(P);

	print "Evaluating central characters";

	centralchars := [ ];
	for i:=1 to #W`CharacterTable do
		chi := W`CharacterTable[i];
		act := [ActionOnVermaModule(chi,z) : z in H`CenterGeneratorsOfDegreeZero ];
		Append(~centralchars, act);
		PrintPercentage(i, #W`CharacterTable);
	end for;

	print "";

	//generic CM families
	H`CalogeroMoserFamilies[1] := { };
	centralcharsset := SequenceToSet(centralchars);
	for v in centralcharsset do
		fam := { i : i in [1..#W`CharacterTable] | centralchars[i] eq v };
		H`CalogeroMoserFamilies[1] join:= { fam };
	end for;



	//determine exceptional locus
	//by theory this is a subscheme of codimension one, stable under C^* action. we conjecture it is in fact a union of rational hyperplanes but we don't know yet. we check if this is indeed the case and quit with an error if not.

	print "Determining exceptional locus";
	N := #H`CenterGeneratorsOfDegreeZero;

	/*print "Components:";

	for twocentralchars in Subsets(centralcharsset, 2) do
		twocentralcharsseq := SetToSequence(twocentralchars);
		variety := [ twocentralcharsseq[1][i] - twocentralcharsseq[2][i] : i in [1..N] ];
		I := ideal<R | variety>;
		IndentPush();
		print I, Dimension(I);
		print RadicalDecomposition(I);
		IndentPop();
		print "";
	end for;*/

	hyperplanes := {};
	if Type(R) ne RngMPol then
		error "Need polynomial ring as base ring to determine exceptional locus.";
	end if;

	for twocentralchars in Subsets(centralcharsset, 2) do
		twocentralcharsseq := SetToSequence(twocentralchars);
		variety := [ twocentralcharsseq[1][i] - twocentralcharsseq[2][i] : i in [1..N] ];
		I := ideal<R | variety>;

		components := RadicalDecomposition(I);

		for J in components do
			//only consider codimension one components
			if Dimension(J) ne Rank(R) - 1 then
				continue;
			end if;
			if #Basis(J) ne 1 then
				error "Something wrong.";
			end if;
			f := Basis(J)[1];
			if not (IsHomogeneous(f) and Degree(f) eq 1) then
				error "Exceptional locus not union of hyperplanes!";
			end if;
			//f := P!NormalizeRationalHyperplaneEquation(f);
			hyperplanes join:={f};
		end for;

	end for;

	//now, determine families on the generic points of the hyperplanes
	print "Determining families";
	for f in hyperplanes do
		I := ideal<R|f>;
		Q,q := R/I;
		specializedcentralchars := [ [q(z[i]) : i in [1..N]] : z in centralchars ];
		specializedcentralcharsset := { z : z in specializedcentralchars};
		fams := { { i : i in [1..#specializedcentralchars] | specializedcentralchars[i] eq z } : z in specializedcentralcharsset };
		fP := P!NormalizeRationalHyperplaneEquation(P!f);
		H`CalogeroMoserFamilies[fP] := fams;
	end for;

	if SaveToDB then
		gen, type := IsGeneric(H);
		if not gen or type ne "GGOR" then
			error "Saving only for generic GGOR type enabled";
		end if;
		//check if Gordon record exists, compare and/or fill up
		res := CHAMP_ExistsInDB(W`DBDir, "Gordon");
		if res then
			print "Found Gordon record in DB";
			gordon := CHAMP_GetFromDB(W`DBDir, "Cherednik/Gordon");
			gordonkeys := {P!f : f in Keys(gordon)};
			if gordonkeys eq {1} then
				print "In existing Gordon record only the generic situation is covered";
				if not gordon[1]`CMFamilies eq H`CalogeroMoserFamilies[1] then
					error "Families do not coincide";
				else
					print "Generic families coincide";
				end if;
			else
				if gordonkeys ne Keys(H`CalogeroMoserFamilies) then
					error "Computed exceptional locus does not coincide with Gordon record!";
				end if;
				print "Exceptional loci coincide";
				for f in Keys(gordon) do
					if not gordon[f]`CMFamilies eq H`CalogeroMoserFamilies[P!f] then
						error "Families do not coincide on "*Sprint(f);
					end if;
				end for;
				print "Families coincide, not saving anything";
				return;
			end if;
		else
			print "No Gordon record found in DB, creating new one";
			gordon := AssociativeArray(P);
		end if;

		hyperplanes join:={1};

		//add data
		print "Saving data";
		GordonRec := recformat<SimpleDims, SimplePSeries, SimpleGModStruct, SimpleGradedGModStruct, VermaDecomposition, CMFamilies, CuspidalCMFamilies>;
		for f in hyperplanes do
			fP := P!NormalizeRationalHyperplaneEquation(P!f);
			if fP notin Keys(gordon) then
				gordon[fP] := rec<GordonRec|>;
			end if;
			gordon[fP]`CMFamilies := H`CalogeroMoserFamilies[fP];
		end for;

		WriteGordonRecord(W, gordon);

	end if;

end intrinsic;

intrinsic CalogeroMoserFamilies(H::AlgChe) -> SetEnum
{}

	CalogeroMoserFamilies(~H);
	return H`CalogeroMoserFamilies;

end intrinsic;

intrinsic CalogeroMoserHyperplanes(H::AlgChe) -> SetEnum
{}

	CalogeroMoserFamilies(~H);
	return SetToSequence(Keys(H`CalogeroMoserFamilies) diff {1});

end intrinsic;

intrinsic CalogeroMoserHyperplanes(W::GrpMat : Type:="GGOR") -> SetEnum
{}

	ReflectionClasses(~W);
	if #W`ReflectionClasses eq 1 then
		return {1};
	end if;

	if Type eq "GGOR" and CHAMP_ExistsInDB(W`DBDir*"/Cherednik", "Gordon") then
		gord := CHAMP_GetFromDB(W`DBDir*"/Cherednik", "Gordon");
		return gord`BlGen;
	else
		H := RationalCherednikAlgebra(W,0 : Type:=Type);
		return CalogeroMoserHyperplanes(H);
	end if;

end intrinsic;

intrinsic CalogeroMoserFamilies(W::GrpMat : Type:="GGOR") -> SeqEnum
{}

	if Type eq "GGOR" and CHAMP_ExistsInDB(W`DBDir*"/Cherednik", "Gordon") then
		gordon := CHAMP_GetFromDB(W`DBDir*"/Cherednik", "Gordon");
		hyp := CalogeroMoserHyperplanes(W);
		fams := AssociativeArray(Universe(hyp));
		for h in hyp do
			fams[h] := gordon`Data[{h}]`CMFamilies;
		end for;
		fams[1] := gordon`Data[{1}]`CMFamilies;
		return fams;
	else
		H := RationalCherednikAlgebra(W,0 : Type:=Type);
		return CalogeroMoserFamilies(H);
	end if;

end intrinsic;
// 
// intrinsic CalogeroMoserHyperplanesAvailable(W::GrpMat) -> BoolElt
// {}
//
// 	try
// 		hyp := CalogeroMoserHyperplanes(W);
// 		return true;
// 	catch e
// 		return false;
// 	end try;
//
// end intrinsic;
//
// intrinsic CalogeroMoserFamiliesAvailable(W::GrpMat) -> BoolElt
// {}
//
// 	try
// 		fams := CalogeroMoserFamilies(W);
// 		return true;
// 	catch e
// 		return false;
// 	end try;
//
// end intrinsic;
