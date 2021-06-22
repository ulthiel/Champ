/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Intrinsics for affine algebras with action of a group (typically coordinate or coinvariant algebras of a matrix group).
*/


//============================================================================
intrinsic Action(f::RngMPolElt, g::GrpMatElt :
	Dual:=false  //if false, then action on K[V], if true then action on K[V^*]
	) -> RngMPolElt
{
	The action of an element g of a matrix group G on a polynomial in as many variables as the dimension of G.
}

	n := Ncols(g);
	P := Parent(f);

	//direct method
	if Dual then
		return f^Transpose(g^-1);
	else
		return f^g;
	end if;

	//alternative method
/*
	res := Zero(P);
	mons := Monomials(f);
	for mon in mons do
		exp := Exponents(mon);
		mong := ArrayProduct([(P`Generators[i]^exp[i])^g : i in [1..n] | exp[i] ne 0] : OneElement:=P`One);
		res +:= MonomialCoefficient(f,mon)*mong;
	end for;

	return res;
*/

end intrinsic;

//============================================================================
intrinsic Action(f::RngMPolResElt, g::GrpMatElt : Dual:=false) -> RngMPolResElt
{
	If f is an element of the affine algebra A=P/I and g acts on P, then compute the action of g on f. Make sure that I is g-invariant, this is not checked.
}

	//there are two options here: we can lift f back to P, compute the action and project, or we can compute this directly in A. The latter is more efficient.

	G := Parent(g);

	//quickest method
	n := Dimension(G);
	A := Parent(f);
	P := OriginalRing(A);
	mons := Monomials(f);
	res := Zero(A);

	if Dual then
		g := Transpose(g^-1);
	end if;

	for m in mons do
		exp := Exponents(m);
		res +:= MonomialCoefficient(f,m)*ArrayProduct([ (A!(P`Generators[i]^g))^exp[i] : i in [1..n]]);
	end for;

	return res;

	//lift method
/*
	A := Parent(f);
	P := OriginalRing(A);
	return A!(Action(P!f,g : Dual:=Dual));
*/

end intrinsic;

//============================================================================
intrinsic GradedGModule(A::RngMPolRes, G::GrpMat : Rep:="Sparse", Verbose:=false, Dual:=false) -> ModGr
{
	If A=P/I is a finite-dimensional graded quotient of the coordinate algebra of G, return the corresponding module of type ModGrGrp.
}

	Generators(~G);
	Basis(~A);
	K := BaseRing(A);

	//group algebra sits in degree zero
	algdegrees := [ 0 : i in [1..Ngens(G)] ];

	//initialize operation matrices for each component
	if Rep eq "Dense" then
		opmats := [* <i,d,ZeroMatrix(K,A`ComponentDimension[d],A`ComponentDimension[d])> : i in [1..Ngens(G)], d in A`Support *];
	else
		opmats := [* <i,d,SparseMatrix(K,A`ComponentDimension[d],A`ComponentDimension[d])> : i in [1..Ngens(G)], d in A`Support *];
	end if;

	indices := [ <o[1],o[2]> : o in opmats ];

	total := Ngens(G)*#A`Support;	//total number of operations
	count := 0;	// operation counter
	for i in [1..Ngens(G)] do
		for d in A`Support do
			if Verbose then
				count +:= 1;
				PrintPercentage(count, total);
			end if;
			//operation of g_i on homogeneous component of degree A`Degrees[j]
			p := Position(indices, <i,d>);
			for k in [1..A`ComponentDimension[d]] do
				//operation of g_i on k-th basis vector of homogeneous component of degree d
				b := A`Basis[A`ComponentBasis[d][k]];
				v := A`VectorSpaceMap(Action(b, G`Generators[i] : Dual:=Dual));
				for l in Support(v) do
					//remap l-th basis vector to basis vector of the j-th component via A`BasisToComponentBasis[l][2]
					opmats[p][3][k][A`BasisToComponentBasis[l][2]] := v[l];
				end for;
			end for;
		end for;
	end for;

	if Verbose then
		print "";
	end if;

	return GradedGModule(BaseRing(A),G,[ <d,A`ComponentDimension[d]> : d in A`Support], opmats : Rep:=Rep);

end intrinsic;
