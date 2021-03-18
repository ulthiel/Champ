//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Run GAP3 programs and get results.
//
//==============================================================================

//==============================================================================
// Generic function to run GAP3 code
//==============================================================================
intrinsic GAP3(code::MonStgElt) -> MonStgElt
{This runs the given code in GAP3 and returns the final computation result as a string. The GAP3 command has to be set via the system environment variable CHAMP_GAP3.}

	in_file_name := "/tmp/CHAMP_GAP.g";
	out_file_name := "/tmp/CHAMP_GAP.txt";
	code *:= "\nPrintTo(\""*out_file_name*"\", last);";
	Write(in_file_name, code : Overwrite:=true);
	gap_cmd := GetEnv("CHAMP_GAP3");
	ret := System(gap_cmd*" -q < "*in_file_name*" > /dev/null");
	res := Read(out_file_name);
	return res;

end intrinsic;

//==============================================================================
// Wrapper for ComplexReflectionGroup
//==============================================================================
intrinsic GAP3_ComplexReflectionGroup(n::RngIntElt) -> GrpMat
{Returns ComplexReflectionGroup(n) from GAP3 (this is the exceptional complex reflection group).}

	GAP3_Code := Sprintf("ComplexReflectionGroup(%o);", Sprint(n));
	mats_seq := eval "E := func<n | RootOfUnity(n)>; return "*GAP3(GAP3_Code*"\nlast.matgens;");

	//We transpose matrices since GAP acts from the left
	mats := [ Transpose(Matrix(M)) : M in mats_seq ];

	//Change integral matrix groups to rationals
	A := Universe(mats);
	K := BaseRing(A);
	if Type(K) eq RngInt then
		n := Degree(A);
		B := MatrixAlgebra(Rationals(),n);
		ChangeUniverse(~mats, B);
	end if;

	W := MatrixGroup(mats);
	W`GAP3_Code := GAP3_Code;

	return W;

end intrinsic;

intrinsic GAP3_ComplexReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt) -> GrpMat
{Returns ComplexReflectionGroup(m,p,n) from GAP3.}

	GAP3_Code := Sprintf("ComplexReflectionGroup(%o,%o,%o);", Sprint(m),Sprint(p),Sprint(n));
	mats_seq := eval "E := func<n | RootOfUnity(n)>; return "*GAP3(GAP3_Code*"\nlast.matgens;");

	//We transpose matrices since GAP acts from the left
	mats := [ Transpose(Matrix(M)) : M in mats_seq ];

	//Change integral matrix groups to rationals
	A := Universe(mats);
	K := BaseRing(A);
	if Type(K) eq RngInt then
		n := Degree(A);
		B := MatrixAlgebra(Rationals(),n);
		ChangeUniverse(~mats, B);
	end if;

	W := MatrixGroup(mats);
	W`GAP3_Code := GAP3_Code;

	return W;

end intrinsic;

//==============================================================================
// Import of ComplexReflectionGroup to the CHAMP database
//==============================================================================
intrinsic GAP3_ImportComplexReflectionGroup(n::RngIntElt)
{}

	W_GAP := GAP3_ComplexReflectionGroup(n);

	// Test whether group is already in the database.
	// If so, check if this is the same model. Error if not.
	try
		W_CHAMP := ComplexReflectionGroup(n);
		if [ W_GAP.i : i in [1..Ngens(W_GAP)] ] eq [ W_CHAMP.i : i in [1..Ngens(W_CHAMP)] ] then
			print "Group is already in the database.";
		else
			error "A different model is in the database.";
		end if;
	catch e
		name := ComplexReflectionGroupName(n)*"_CHEVIE";
		CHAMP_SaveToDB(Sprint(W_GAP, "Magma"), "ReflectionGroups", name);
		print "Group imported.";
	end try;

end intrinsic;

intrinsic GAP3_ImportComplexReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt)
{}

	W_GAP := GAP3_ComplexReflectionGroup(m,p,n);

	// Test whether group is already in the database.
	// If so, check if this is the same model. Error if not.
	try
		W_CHAMP := ComplexReflectionGroup(m,p,n);
		if [ W_GAP.i : i in [1..Ngens(W_GAP)] ] eq [ W_CHAMP.i : i in [1..Ngens(W_CHAMP)] ] then
			print "Group is already in the database.";
		else
			error "A different model is in the database.";
		end if;
	catch e
		name := ComplexReflectionGroupName(m,p,n)*"_CHEVIE";
		CHAMP_SaveToDB(Sprint(W_GAP, "Magma"), "ReflectionGroups", name);
		print "Group imported.";
	end try;

end intrinsic;

intrinsic GAP3_ChevieClassData(~W::GrpMat)
{This function assumes that W is a ComplexReflectionGroup and assigns the class data that is in ChevieClassInfo to the corresponding attributes of W.}

	if not assigned W`GAP3_Code then
		error "Group needs to have attribute GAP3_Code assigned to construct it";
	end if;

	GAP3_Code := "W:="*W`GAP3_Code*"\n";
	class_words := eval GAP3(GAP3_Code*"\nChevieClassInfo(W).classtext;");
	class_names := eval GAP3(GAP3_Code*"\nChevieClassInfo(W).classnames;");

	class_reps := [ WordToElement(W,x) : x in class_words ];

	if not assigned W`Classes then
		class_reps := [ WordToElement(W,x) : x in class_words ];
		W`Classes := class_reps;
	end if;

	sigma := ClassPermutation(W, class_words);

	W`ClassWords := [ class_words[sigma[i]] : i in [1..#sigma] ];
	W`ClassNames := [ class_names[sigma[i]] : i in [1..#sigma] ];

end intrinsic;
