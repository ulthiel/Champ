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
for n:= 2 to 10 do
	W := TypeAReflectionGroup(n);
	assert IsIsomorphic(W, SymmetricGroup(n));
	assert Dimension(W) eq n-1;
end for;

print("Testing type B");
for n:= 2 to 10 do
	W := TypeBReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(2,1,n));
end for;

print("Testing type D");
for n:= 2 to 10 do
	W := TypeDReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(2,2,n));
end for;

print("Testing dihedrals");
for n:= 3 to 10 do
	W := DihedralReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(n,n,2));
end for;

print("Testing cyclic");
for n:= 3 to 10 do
	W := CyclicReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(n,1,1));
end for;

print("Testing exceptionals");
for n:= 4 to 37 do
	W := ExceptionalComplexReflectionGroup(n);
	assert IsIsomorphic(W, ShephardTodd(n));
end for;

quit;
