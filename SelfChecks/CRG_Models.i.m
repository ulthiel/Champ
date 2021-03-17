//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Test models of complex reflection groups.
//
//==============================================================================

print("Testing type A");
for n:= 2 to 25 do
	W := TypeAReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(1,1,n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
	if n le 10 then
		assert IsIsomorphic(W, ShephardTodd(1,1,n));
	end if;
end for;

print("Testing type B");
for n:= 2 to 25 do
	W := TypeBReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(2,1,n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
	if n le 10 then
		assert IsIsomorphic(W, ShephardTodd(2,1,n));
	end if;
end for;

print("Testing type D");
for n:= 2 to 25 do
	W := TypeDReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(2,2,n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
	if n le 10 then
		assert IsIsomorphic(W, ShephardTodd(2,2,n));
	end if;
end for;

print("Testing dihedrals");
for n:= 3 to 12 do
	W := DihedralReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(n,n,2));
end for;

print("Testing cyclic");
for n:= 3 to 12 do
	W := CyclicReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(n,1,1));
end for;

print("Testing exceptionals");
for n:= 4 to 37 do
	W := ExceptionalReflectionGroup(n);
	W_GAP := GAP3_ComplexReflectionGroup(n);
	assert [ W.i : i in [1..Ngens(W)] ] eq [ W_GAP.i : i in [1..Ngens(W_GAP)] ];
	assert IsIsomorphic(W, ShephardTodd(n));
end for;

quit;
