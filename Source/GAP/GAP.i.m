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

intrinsic GAP3(code::MonStgElt) -> MonStgElt
{This runs the given code in GAP3. The final computation result in the code should be a string assigned to the variable MAGMA_RESULT. This is then returned by the function. The GAP3 command has to be set via the environment variable CHAMP_GAP3.}

	in_file_name := "/tmp/CHAMP_GAP.g";
	out_file_name := "/tmp/CHAMP_GAP.txt";
	code *:= "\nPrintTo(\""*out_file_name*"\", MAGMA_RESULT);";
	Write(in_file_name, code : Overwrite:=true);
	gap_cmd := GetEnv("CHAMP_GAP3");
	ret := System(gap_cmd*" -q < "*in_file_name*" > /dev/null");
	res := Read(out_file_name);
	return res;

end intrinsic;

intrinsic GAP3_ComplexReflectionGroup(n::RngIntElt) -> GrpMat
{Returns ComplexReflectionGroup(n) from GAP3 (this is the exceptional complex reflection group).}

	mats_str := GAP3(Sprintf("MAGMA_RESULT:=ComplexReflectionGroup(%o).matgens;", Sprint(n)));
	mats_seq := eval "E := func<n | RootOfUnity(n)>; return "*mats_str;

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
	return W;

end intrinsic;

intrinsic GAP3_ComplexReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt) -> GrpMat
{Returns ComplexReflectionGroup(m,p,n) from GAP3.}

	mats_str := GAP3(Sprintf("MAGMA_RESULT:=ComplexReflectionGroup(%o,%o,%o).matgens;", Sprint(m),Sprint(p),Sprint(n)));
	mats_seq := eval "E := func<n | RootOfUnity(n)>; return "*mats_str;

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
	return W;

end intrinsic;
