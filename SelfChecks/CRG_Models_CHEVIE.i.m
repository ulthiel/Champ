//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Test whether models of complex reflection groups are equal to the CHEVIE
// ones (uses GAP3 interface, so make sure GAP3 works and the environment
// variable CHAMP_GAP3 for the GAP3 command is set).
//
//==============================================================================

print("Testing type A");
for n:= 1 to 25 do
	W := TypeAReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(1,1,n+1);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

print("Testing type B");
for n:= 2 to 25 do
	W := TypeBReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(2,1,n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

print("Testing type D");
for n:= 4 to 25 do
	W := TypeDReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(2,2,n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

print("Testing dihedrals");
for m:= 5 to 25 do
	W := DihedralReflectionGroup(m);
	W_GAP := GAP3_ComplexReflectionGroup(m,m,2);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

print("Testing cyclic");
for m:= 3 to 25 do
	W := CyclicReflectionGroup(m);
	W_GAP := GAP3_ComplexReflectionGroup(m,1,1);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

print("Testing exceptionals");
for n:= 4 to 37 do
	W := ComplexReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
end for;

quit;
