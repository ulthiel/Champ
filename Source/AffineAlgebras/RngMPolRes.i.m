/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Simple extensions for affine algebras (algebras of finite type over a field, type RngMPolRes).
*/


declare attributes RngMPolRes:
	Zero,
	One,
    Basis, // A basis of the algebra sorted by degree. This will usually be MonomialBasis but may be differently assigned by user.
    MonomialBasis, //basis selected by MonomialBasis (this is the basis used by Magma internally)
    Generators, //Contains the sequence of generators.
    PoincareSeries, //The Poincare series of the algebra.
    StructureConstantAlgebra, //The corresponding structure constant algebra.
    VectorSpace, //The corresponding vector space.
    VectorSpaceMap, //The map from the algebra into its vector space.
    Support, //The sorted list of the support of A
	ComponentBasis, //associative array such that d-th entry is the indexed set of numbers of basis vectors of the homogeneous component of degree d
	ComponentDimension, //associative array with dimensions of the homogeneous components
	BasisToComponentBasis, //sequence such that i-th entry is <d,k> where d is the degree of the i-th basis element and k is the index of the basis element within the homogeneous component of degree d
	GradedVectorSpace, //associative array indexed by support so that d-th entry is the vector space corresponding to the homogeneous component of degree d
	GradedVectorSpaceMap;

//============================================================================
intrinsic Generators(~A::RngMPolRes)
{
	Attaches the sequence of generators to the corresponding attribute of P. This allows faster acces than using the dot operator.
}

    if assigned A`Generators then
        return;
    else
        A`Generators := [ A.i : i in [1..Ngens(A)]];
    end if;

end intrinsic;


//============================================================================
intrinsic MonomialBasis(~A::RngMPolRes : SortByDegree:=true)
{
	Assigns a monomial basis computed with MonomialBasis to the corresponding attribute and sorts this basis by degree. Furthermore, VectorSpace and VectorSpaceMap is assigned.
}

	//we need this anyways to get VectorSpaceMap (which will be modified according to A`Basis)
	if not assigned A`MonomialBasis then
		A`MonomialBasis := MonomialBasis(A);

		//all the following is independent of the chosen basis
		A`Support := RemoveDuplicates(Sort([ Degree(m) : m in A`MonomialBasis ]));
		A`ComponentDimension := AssociativeArray(A`Support);
		A`GradedVectorSpace := AssociativeArray(A`Support);
		for d in A`Support do
			currentcomponent := {@ i : i in [1..#A`MonomialBasis] | Degree(A`MonomialBasis[i]) eq d@};
			A`ComponentDimension[d] := #currentcomponent;
			A`GradedVectorSpace[d] := KSpace(BaseRing(A), #currentcomponent);
		end for;

		//initialize
		A`ComponentBasis := AssociativeArray(A`Support);
		A`BasisToComponentBasis := AssociativeArray({1..#A`MonomialBasis});
	end if;

end intrinsic;


//============================================================================
intrinsic Basis(~A::RngMPolRes : SortByDegree:=true)
{}

	if not assigned A`Basis then
		MonomialBasis(~A);
		if SortByDegree then
			//MonomialBasis is sorted, but not by degree so let's take care of this
			A`Basis := {@@};
			for d in A`Support do
				 A`Basis join:=Sort({@ b : b in A`MonomialBasis | Degree(b) eq d@});
			end for;
		else
			A`Basis := A`MonomialBasis;
		end if;
		SetBasis(~A, A`Basis);
	end if;

end intrinsic;

//============================================================================
intrinsic Basis(A::RngMPolRes : SortByDegree:=true) -> SetIndx
{Finds a vector space basis of A. If SortByDegree is true (default), the basis is sorted by degree.}

    Basis(~A : SortByDegree:=SortByDegree);
    return A`Basis;

end intrinsic;

//============================================================================
intrinsic SetBasis(~A::RngMPolRes, basis::SetIndx)
{Sets VectorSpaceMap and GradedVectorSpace according to currently chosen basis of A.}

	MonomialBasis(~A);
	A`Basis := basis;

	//now we have to take the possible base change into account.
	sigma := [ Position(A`Basis,m) : m in A`MonomialBasis ];
	if sigma eq [1..#A`Basis] then //in this case no change
		A`VectorSpace, A`VectorSpaceMap := VectorSpace(A);
	elif SequenceToSet(sigma) eq {1..#A`Basis} then //permutation
		A`VectorSpace, vmap := VectorSpace(A);
		A`VectorSpaceMap := function(f);
			v := vmap(f);
			w := Zero(A`VectorSpace);
			for i in Support(v) do;
				w[sigma[i]] := v[i];
			end for;
			return w;
		end function;
	else //in this case basis looks completely different; the computation is not that efficient
		A`VectorSpace, vmap := VectorSpace(A);
		phi := hom<A`VectorSpace->A`VectorSpace | [ vmap(A`Basis[i]) : i in [1..#A`Basis] ]>;
		psi := Inverse(phi);
		A`VectorSpaceMap := function(f);
			return psi(vmap(f));
		end function;
	end if;

	for d in A`Support do
		currentcomponent := {@ i : i in [1..#A`Basis] | Degree(A`Basis[i]) eq d@};
		A`ComponentBasis[d] := currentcomponent;
		for i:=1 to #currentcomponent do
			A`BasisToComponentBasis[currentcomponent[i]] := <d,i>;
		end for;
	end for;

	A`GradedVectorSpaceMap := function(f);
		v := A`VectorSpaceMap(f);
		supp := Support(v);
		degs := { A`BasisToComponentBasis[i][1] : i in supp };
		vd := < <d, &+[ v[i]*A`GradedVectorSpace[d].(A`BasisToComponentBasis[i][2]) : i in supp | A`BasisToComponentBasis[i][1] eq d]> : d in degs >;
		return vd;
	end function;

	//now, do a self check
	for i:=1 to #A`Basis do
		assert A`VectorSpaceMap(A`Basis[i]) eq A`VectorSpace.i;
		d := A`BasisToComponentBasis[i][1];
		k := A`BasisToComponentBasis[i][2];
		assert Degree(A`Basis[i]) eq d;
	end for;

end intrinsic;

//============================================================================
intrinsic StructureConstantAlgebra(~A::RngMPolRes : Rep:="Sparse")
{The structure constant algebra of A with respect to the assigned basis.}

    if assigned A`StructureConstantAlgebra then
		return;
	end if;

	VectorSpace(~A);
	Q := [];
	for i:=1 to #A`Basis do
		for j:=1 to #A`Basis do
			v := A`VectorSpaceMap(A`Basis[i]*A`Basis[j]);
			for k in Support(v) do
				Append(~Q, <i,j,k,v[k]>);
			end for;
		end for;
	end for;
	A`StructureConstantAlgebra := AssociativeAlgebra<BaseRing(A), #A`Basis | Q : Rep:="Sparse", Check:=false>;

end intrinsic;

//============================================================================
intrinsic StructureConstantAlgebra(A::RngMPolRes : Rep:="Sparse") -> AlgAss
{The structure constant algebra of A.}

    StructureConstantAlgebra(~A : Rep:=Rep);
    return A`StructureConstantAlgebra;

end intrinsic;


//============================================================================
intrinsic PoincareSeries(~A::RngMPolRes)
{The Poincare series of A.}

    A`PoincareSeries := HilbertSeries(A);

end intrinsic;

//============================================================================
intrinsic PoincareSeries(A::RngMPolRes) -> FldFunRatUElt
{}

    PoincareSeries(~A);
    return A`PoincareSeries;

end intrinsic;


//============================================================================
intrinsic IsMonomial(f::RngMPolResElt) -> BoolElt
{True iff f is represented as a monomial.}

    return #Monomials(f) eq 1;

end intrinsic;


//============================================================================
intrinsic GradedModule(A::RngMPolRes : Rep:="Sparse", Verbose:=false) -> ModGr
{
	If A=P/I is a finite-dimensional graded affine algebra, return the corresponding module of type ModGr.
}

	Basis(~A);
	K := BaseRing(A);
	Generators(~A);

	//group algebra sits in degree zero
	algdegrees := [ Degree(A`Generators[i]) : i in [1..Ngens(A)] ];

	//initialize operation matrices for each component
	if Rep eq "Dense" then
		opmats := [* <i,d,ZeroMatrix(K,A`ComponentDimension[d],A`ComponentDimension[d+algdegrees[i]])> : i in [1..Ngens(A)], d in A`Support | d+algdegrees[i] in A`Support *];
	else
		opmats := [* <i,d,SparseMatrix(K,A`ComponentDimension[d],A`ComponentDimension[d+algdegrees[i]])> : i in [1..Ngens(A)], d in A`Support | d+algdegrees[i] in A`Support *];
	end if;

	indices := [ <o[1],o[2]> : o in opmats ];

	total := Ngens(A)*#A`Support;
	count := 0;
	for i in [1..Ngens(A)] do
		for d in A`Support do
			if Verbose then
				count +:= 1;
				PrintPercentage(count, total);
			end if;
			//operation of x_i on homogeneous component of degree A`Degrees[j]
			p := Position(indices, <i,d>);
			for k in [1..A`ComponentDimension[d]] do
				//operation of g_i on k-th basis vector of homogeneous component of degree d
				b := A`Basis[A`ComponentBasis[d][k]];
				v := A`VectorSpaceMap(b*A`Generators[i]);
				//w := Zero(A`GradedVectorSpace[j]);
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

	return GradedModule(BaseRing(A),algdegrees,[ <d,A`ComponentDimension[d]> : d in A`Support], opmats : Rep:=Rep);

end intrinsic;
