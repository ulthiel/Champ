/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	A new graded module structure.
*/


declare type ModGr[ModGrElt];

declare attributes ModGr:
	Rep,	//"Dense" or "Sparse"
    BaseRing,			// Basering R
    ModuleComponents,		// Associative Array - List of the components M_i in M,
				//			Key elements : Degree of the components
    Support,			// List of the nonzero grades of the module components
    AlgebraDegrees,		// List of the degrees of the algebra generators
    Dimension,			// Dimension of the module M
    ComponentDimension,		// List of the dimensions of the components
    Matrices,		// List of AssociativeArrays with the support of the graded module as key-elements. The i-th array element contains matrices how the ith algebra generator operates on the module component with degree d.
    Zero,	// the Zero element (avoids creation every time)
    FullSpace, //this is the underlying full vector space of the module (i.e., the direct sum of all components). Coercion from an element here in graded module creates the decomposition into components.
    FullSpaceBasisToComponentBasis; //i-th basis vector of full space corresponds to j-th basis vector in degree d component. The i-th entry of this array is the tuple <d,k>. This is used for coercion of FullSpace elements into graded module.

declare attributes ModGrElt:
	Parent,			// Parent of the the vector
	Components,		// AssociativeArray with the degree of the components as key elements
				//	- Contains the homogenous parts of the vector
	Support;		// List with the nonzero degrees of the vector



//============================================================================
intrinsic GradedModule(R::Rng, AlgebraDegrees::SeqEnum[RngIntElt], Components::SeqEnum[Tup], Q::List : Rep:="Sparse") -> ModGr

{
Graded module over basering R.

Parameters:
R: Basering or basefield;
AlgebraDegrees: List of degrees of the algebra generators;
Components: Lists of 2-tuples: First entry is the degree, second entry is the dimension of the component;
Q: The action of the algebra generators on the the components.
}

	M := New(ModGr);

	//M`Basering
	M`BaseRing := R;

	//M`AlgebraDegrees
	M`AlgebraDegrees := AlgebraDegrees;

	// M`Support
	M`Support := [ Components[i][1] : i in [1..#Components] ];

	// M`ComponentDimension
	M`ComponentDimension := [ Components[i][2] : i in [1..#Components] ];

	// Set matrix representation type
	M`Rep := Rep;

	// Full dimension
	M`Dimension := ArraySum(M`ComponentDimension : ZeroElement:=Zero(Integers()));

	// Full space (direct sum of all spaces)
	if ISA(Type(M`BaseRing), Fld) then
		M`FullSpace := KSpace(M`BaseRing, M`Dimension);
	else
		M`FullSpace := RSpace(M`BaseRing, M`Dimension);
	end if;

	// M`ModuleComponents
	M`ModuleComponents := AssociativeArray(M`Support);

	if ISA(Type(M`BaseRing), Fld) eq true then	//Checks whether BaseRing R is a field or not.
		for i:=1 to #Components do
			d := M`Support[i];
			dim := M`ComponentDimension[i];
			M`ModuleComponents[d] := KSpace(M`BaseRing,dim); // KSpace
			end for;

	 elif ISA(Type(M`BaseRing), Rng) eq true then	// Checks whether BaseRing R is a ring or not.
		for i:=1 to #Components do
			d := M`Support[i];
			dim := M`ComponentDimension[i];
			M`ModuleComponents[d] := RSpace(M`BaseRing,dim); // RSpace
		end for;
	else
		error "Input "*Sprint(R)*" is not a Ring or Field.";
	end if;

	// M`FullSpaceBasisToComponentBasis
	M`FullSpaceBasisToComponentBasis := [ <M`Support[i],M`ComponentDimension[i]> : i in [1..#M`Support] ];

	// M`Matrices
	M`Matrices := < AssociativeArray(M`Support) : i in [1..#M`AlgebraDegrees] >;

	for q in Q do
		i:= q[1];
		d:= q[2];
		m:= q[3];
		if ISA( Type(m), Mtrx ) eq true and Rep eq "Dense" then
					M`Matrices[i][d] := m;
				elif ISA( Type(m), Mtrx ) eq true and Rep eq "Sparse" then
					M`Matrices[i][d] := SparseMatrix(m);
				elif ISA( Type(m), MtrxSprs ) eq true and Rep eq "Sparse" then
					M`Matrices[i][d] := m;
				elif  ISA( Type(m), MtrxSprs ) eq true and Rep eq "Dense" then
					M`Matrices[i][d] := Matrix(m);
				else
					if ExtendedType(m) eq SeqEnum[Tup] then	// If the input for Matrices are sequences
						M`Matrices[i][d] := SparseMatrix(M`BaseRing, ComponentDimension(M,d+M`AlgebraDegrees[i]), M`ComponentDimension[d], m );
						if Rep eq "Dense" then
							M`Matrices[i][d] := Matrix(M`Matrices[i][d]);
						end if;

					else
						M`Matrices[i][d] := Matrix(M`BaseRing, ComponentDimension(M,d+M`AlgebraDegrees[i]), M`ComponentDimension[d], m );

						if Rep eq "Sparse" then
							M`Matrices[i][d] := SparseMatrix(M`Matrices[i][d]);
						end if;
					end if;
				end if;
	end for;

	// Set Zero element
	M`Zero := New(ModGrElt);
	M`Zero`Parent := M;
	M`Zero`Support := [];
	M`Zero`Components := AssociativeArray(M`Support);

	return M;

end intrinsic;

//============================================================================
intrinsic ComponentDimension(M::ModGr, d::RngIntElt) -> RngIntElt
{Returns the dimension of the module component of grade d of a graded module M.}

	return Dimension(M`ModuleComponents[d]);

end intrinsic;

intrinsic Support(M::ModGr) -> SeqEnum
{}

	return M`Support;

end intrinsic;

intrinsic BaseRing(M::ModGr) -> Rng
{}

	return M`BaseRing;

end intrinsic;

//============================================================================
intrinsic VectorComponent(v::ModGrElt, d::RngIntElt) -> ModTupFldElt
{Returns the component of a vector of a graded module with degree d.}

	if not d in v`Support then
		return Zero(v`Parent);
	else
		return v`Components[d];
	end if;

end intrinsic;

//============================================================================
intrinsic Degree(v::ModGrElt) -> RngIntElt
{}

	d := v`Support;
	if #d eq 1 then
		return d[1];
	else
		error "Not homogeneous";
	end if;

end intrinsic;

//============================================================================
intrinsic SupportSort(M::ModGr, support::SeqEnum) -> SeqEnum
{Sorts a (support-)sequence in the order corresponding to the order given in the graded module M.}

	if not #SequenceToSet(support) eq #support then
		error "Sequence "*Sprint(support)*" contains integers of the same value.";
	end if;
	for d in support do
		if not d in M`Support then
			error "Degree "*Sprint(d)*" is not contained in the support of the graded module.";
		end if;
	end for;
	s:= Sort(support, func<x,y| Position(M`Support,x)-Position(M`Support,y)> );
	return s;

end intrinsic;

//============================================================================
intrinsic IsCoercible(M::ModGr, x::.) -> BoolElt, .
{Checks whether x is coercible into M and returns the result if so.}

	//coercion of element of full space
	if Type(x) eq ModTupFldElt or Type(x) eq ModTupRngElt then
		if not Parent(x) eq M`FullSpace then
			return false, "Illegal coercion.";
		end if;
		v := New(ModGrElt);
		supp := {<i,M`FullSpaceBasisToComponentBasis[i]> : i in Support(x)}; //entries of the form <i,<d,j>> where i in support of full space and <d,j> is degree of e_i j is the corresponding basis element in M_d
		v`Support := SetToSequence({ s[2][1] : s in supp }); //degrees in support
		v`Components := AssociativeArray(M`Support);
		for d in v`Support do
			vd := Zero(M`ModuleComponents[d]);
			for s in { s : s in supp | s[2][1] eq d } do
				i := s[1];
				j := s[2][2];
				vd[j] := x[i];
			end for;
			v`Components[d] := vd;
		end for;
		return true, v;
	end if;

	if ISA(Type(x),Tup) eq true then		// x is a tuple

		for i:=1 to #x do
			if not #x[i] eq 2 then
				return false, " Illegal Coersion.";
			end if;
		end for;

		Degrees := [x[i][1] : i in [1..#x]];
		Components := [* x[i][2] : i in [1..#x]*];

		if not SequenceToSet(Degrees) subset SequenceToSet(M`Support) then
			return false, "Degrees do not match.";
		end if;

		for i:=1 to #x do		// Checks whether the components of x is coersible in M
			d := Degrees[i];
			v := Components[i];
			if not IsCoercible(M`ModuleComponents[d],v) then
				return false, "Basering or dimension of the components do not match.";
			end if;
		end for;

		v := New(ModGrElt);	// Creates the new graded module element
		v`Parent := M;
		v`Support := Degrees;
		v`Components := AssociativeArray(v`Support);
		for i:= 1 to #x do
			mdeg := x[i][1];
			v`Components[mdeg] := x[i][2];
		end for;
		return true, v;

	elif ISA(Type(x),ModGrElt) eq true then		// x is a graded module element
		if not SequenceToSet(x`Support) subset SequenceToSet(M`Support) then
			return false, "Degrees do not match.";
		else
			for d in x`Support do
				if not IsCoercible(M`ModuleComponents[d], x`Components[d]) then
					return false, "Basering or dimension of the components do not match.";
				end if;
			end for;
			return true, x;
		end if;
	else
		return false, "Illegal coercion.";
	end if;
end intrinsic;

//============================================================================
intrinsic Parent(m::ModGrElt) -> ModGr
{}

	return m`Parent;

end intrinsic;

//============================================================================
intrinsic Zero(M::ModGr) -> ModGrElt
{Returns the zero element of a graded module.}

	return M`Zero;

end intrinsic;

//============================================================================
intrinsic Basis(M::ModGr) -> List
{Returns a list with the basis vectors of M.}

	if #M`Support eq 0 then
		return [];
	end if;
	B:= [M.i : i in [1..M`Dimension]];
	return B;

end intrinsic;

//============================================================================
intrinsic '*' (r::RngElt, m::ModGrElt) -> ModGrElt
{Scalar multiplication of an element m of a graded module with an element r of its base ring.}
	if r eq 0 then
		return Zero(m`Parent);
	end if;

	x := New(ModGrElt);
	x`Parent := m`Parent;
	x`Support := m`Support;
	x`Components := m`Components;
	for d in m`Support do
		x`Components[d] *:= r;
	end for;
	return x;
end intrinsic;

intrinsic '*' (m::ModGrElt, r::RngElt) -> ModGrElt
{Scalar multiplication of an element m of a graded module with an element r of its base ring.}
	if r eq 0 then
		return Zero(m`Parent);
	end if;

	x := New(ModGrElt);
	x`Parent := m`Parent;
	x`Support := m`Support;
	x`Components := m`Components;
	for d in m`Support do
		x`Components[d] *:= r;
	end for;
	return x;
end intrinsic;

//============================================================================
intrinsic '+' (n1::ModGrElt, n2::ModGrElt) -> ModGrElt
{Addition of two elements of a graded module.}
	if #n1`Support ge #n2`Support then	// Compares the length of n1 and n2
		m1 := n1;
		m2 := n2;
	else
		m1 := n2;
		m2 := n1;
	end if;

		x := New(ModGrElt);
		x`Parent := m1`Parent;
		x`Support := m1`Parent`Support;
		N:= #x`Support;

		for i in [1..N] do		// Removes the entry from the support, if not in m1 or m2
			d := x`Support[i];
			if  not d in m1`Support and not d in m2`Support then
			Exclude(~x`Support, d);
			end if;
		end for;

		x`Components := AssociativeArray(x`Support);

		for d in x`Support do				// Adding the Components
			if d in m1`Support and d in m2`Support then
				x`Components[d] := m1`Components[d] + m2`Components[d];
			elif d in m1`Support and not d in m2`Support then
				x`Components[d] := m1`Components[d];
			else
				x`Components[d] := m2`Components[d];
			end if;

			if x`Components[d] eq 0 then		// Removes the component if zero
				Remove(~x`Components,d);
				Remove(~x`Support,d);
			end if;
		end for;

	return x;

end intrinsic;

//============================================================================
intrinsic '-' (m1::ModGrElt, m2::ModGrElt) -> ModGrElt
{Subtraction of two elements of a graded module.}

	x := New(ModGrElt);
	x`Parent := m1`Parent;
	x`Support := m1`Parent`Support;
	N := #x`Support;

	for i in [1..N] do				// Support of x
		d:= x`Support[i];
		if not d in m1`Support and not d in m2`Support then
		Exclude(~x`Support, d);
		end if;
	end for;

	x`Components := AssociativeArray(x`Support);

	for d in x`Support do				// Adding the components
		if d in m1`Support and d in m2`Support then
			x`Components[d] := m1`Components[d] - m2`Components[d];
		elif d in m1`Support and not d in m2`Support then
			x`Components[d] := m1`Components[d];
		else
			x`Components[d] := -m2`Components[d];
		end if;
		if x`Components[d] eq 0 then	// Removes the component if zero
			Remove(~x`Components,d);
			Remove(~x`Support,d);
		end if;

	end for;

	return x;

end intrinsic;


//============================================================================
intrinsic '.' (M::ModGr, q::Tup) -> ModGrElt
{For a pair q := <d,j> this intrinsic returns the j-th basis vector of the component of degree d.}
	d := q[1];
	j := q[2];
	x := New(ModGrElt);
	x`Parent := M;
	x`Support:= [d];
	x`Components := AssociativeArray(x`Support);
	x`Components[d] := Basis(M`ModuleComponents[d])[j];

	return x;

end intrinsic;


//============================================================================
intrinsic '.' (M::ModGr, j::RngIntElt) -> ModGrElt
{Returns the j-th basis vector of M via counting the dimensions in the same order
as in ModuleComponents.}

	if not j in [1 .. M`Dimension] then
		error "Argument 2 ("*Sprint(j)*") should be in the range "*Sprint([1 .. M`Dimension])*".";
	end if;

	x := New(ModGrElt);
	x`Parent := M;

	pos := 1;			//current position
	a := 0;				//lower bond
	b := M`ComponentDimension[pos];	//upper bond
	flag := false;

	while flag eq false do
		if j gt a and j le b then		// Looking for the position of j in ComponentDimension
			d := M`Support[pos];
			x`Support := [d];
			x`Components := AssociativeArray(x`Support);
			x`Components[d] := Basis(M`ModuleComponents[d])[j-a];
			flag := true;
		else					// Shifts position +1
			a +:= M`ComponentDimension[pos];
			b +:= M`ComponentDimension[pos+1];
			pos +:= 1;
		end if;
	end while;
	return x;
end intrinsic;


//============================================================================
intrinsic Print(M::ModGr)
{Prints a graded module}

	if M`Dimension eq 0 then
		   	printf "Graded module of dimension 0 and components "*Sprint( [<M`Support[i],M`ComponentDimension[i] > : i in [1..#M`Support]] )*" over algebra with generator degrees "*Sprint(M`AlgebraDegrees)*" over "*Sprint(M`BaseRing);
	else

   	printf "Graded module of dimension "*Sprint(&+[M`ComponentDimension[i] : i in [1..#M`Support]])*" and components "*Sprint( [<M`Support[i],M`ComponentDimension[i] > : i in [1..#M`Support]] )*" over algebra with generator degrees "*Sprint(M`AlgebraDegrees)*" over "*Sprint(M`BaseRing);
	end if;

end intrinsic;

//============================================================================
intrinsic Print(m::ModGrElt)
{Prints an element of a graded module M}
   printf " "*Sprint(< <d,m`Components[d]> : d in SupportSort(m`Parent, m`Support) >)*" ";

end intrinsic;

//============================================================================
intrinsic '^' (m::ModGrElt, i::RngIntElt) -> ModGrElt
{Operation of the i-th Algebra generator on the graded module element m.}


	M := m`Parent;
	adeg := M`AlgebraDegrees[i];				// Degree of the operating algebra generator
	S := [ d+adeg : d in m`Support | d+adeg in M`Support ];	// Support of x=m*A.i

	if #S eq 0 then			// Checks whether support of x empty
		return Zero(M);
	end if;

	x := New(ModGrElt);
	x`Parent := M;
	x`Support := S;
	x`Components := AssociativeArray(x`Support);

	for d in x`Support do
		x`Components[d] := m`Components[d-adeg]*M`Matrices[i][d-adeg];
		if x`Components[d] eq Zero(M`ModuleComponents[d]) then	// Removes the component if zero
			Remove(~x`Components,d);
			Exclude(~x`Support,d);
		end if;
	end for;

	return x;

end intrinsic;


//============================================================================
intrinsic 'eq' (m1::ModGrElt, m2:ModGrElt) -> BoolElt
{Compares two elements m1 and m2 of a graded module M. Returns true if they are equal and false if they are not.}

	if not m1`Support eq m2`Support then
		return false;
	end if;
	for d in m1`Support do
		if not m1`Components[d] eq m2`Components[d] then
			return false;
		end if;
	end for;
	return true;
end intrinsic;





//============================================================================
//	Spinning
//============================================================================


//============================================================================
intrinsic Spin(M::ModGr, U::List : Verbose:=true) -> List
{Returns a basis of the graded submodule of a graded module M generated by homogenous elements listet in U.}

	for m in U do
		if not #m`Support eq 1 then
			error "Vector "*Sprint(m)*" is not homogenous.";
		end if;
	end for;

	//this will be the span
	span := AssociativeArray(M`Support);

	for d in M`Support do
		Ud := {u`Components[d] : u in U | u`Support[1] eq d};
		span[d] := sub<M`ModuleComponents[d]|Ud>;
	end for;

	SignAdeg := {Sign(d): d in M`AlgebraDegrees};

	//special treatment when all algebra generators are of non-negative or non-positive degree
	if #SignAdeg le 2 and not SignAdeg eq {-1,1} then
		if SignAdeg subset {0,1} then
			Support := Sort(M`Support);
		else
			Support := Reverse(Sort(M`Support));
		end if;

		for d in Sort(M`Support) do
			//we compute U_d
			gens := [ j : j in [1..#M`AlgebraDegrees] | d-M`AlgebraDegrees[j] in M`Support and IsDefined(span,d-M`AlgebraDegrees[j]) and not Dimension(span[d-M`AlgebraDegrees[j]]) eq 0 ]; //all generators which map to U_d
			images := {};
			for j in gens do
				e := M`AlgebraDegrees[j];
				for b in Basis(span[d-e]) do
					v := New(ModGrElt);
					v`Parent := M;
					v`Support := [d-e];
					v`Components := AssociativeArray(M`Support);
					v`Components[d-e] := b;
					w := v^j;
					if #w`Support ne 0 then
						images join:={w`Components[d]};
					end if;
				end for;
			end for;
			span[d] +:= sub<M`ModuleComponents[d]|images>;
			//span[d] := sub<M`ModuleComponents[d]|SequenceToSet(Basis(span[d])) join images>; //no time difference to above
		end for;
		B := [* *];
		for d in Keys(span) do
			for b in Basis(span[d]) do
				v := New(ModGrElt);
				v`Parent := M;
				v`Support := [d];
				v`Components := AssociativeArray(M`Support);
				v`Components[d] := b;
				Append(~B,v);
			end for;
		end for;
		return B;
	end if;

	//general version (here we can positive, negative, zero degree generators)

	//this will be our basis
	B := [**];
	for d in Keys(span) do
		for b in Basis(span[d]) do
			m := New(ModGrElt);
			m`Parent := M;
			m`Support := [d];
			m`Components := AssociativeArray(M`Support);
			m`Components[d] := b;
			Append(~B, m);
		end for;
	end for;
	tospin := B;

	//here, we record the components which are of full dimension
	fulldegrees := {d : d in Keys(span) | Dimension(span[d]) eq ComponentDimension(M,d) };

	currentdim := #B;
	if currentdim eq M`Dimension then
		return B;
	end if;

	while not #tospin eq 0 do
		i := 1;//Random({1..#tospin});
		v := tospin[i];
		Remove(~tospin,i);	// Removes the vector that has been spinned.
		for i:=1 to #M`AlgebraDegrees do
			d := v`Support[1];	// note that v homogeneous
			adeg := M`AlgebraDegrees[i];
			s := d+adeg;
			if s notin M`Support then
				continue;
			end if;
			if s in fulldegrees then
				continue;
			end if;
			w := v^i;
			if #w`Support eq 0 then
				continue;
			end if;
			if w`Components[s] in span[s] then
				continue;
			end if;
			Append(~B,w);
			currentdim +:= 1;
			if currentdim eq M`Dimension then	//generate the full module, so can quit
				if Verbose then
					PrintAndDelete("Dimension: "*Sprint(currentdim)*"                                          ");
					print "";
				end if;
				return B;
			end if;
			span[s] +:= sub<M`ModuleComponents[s]|w`Components[s]>;
			Append(~tospin, w);
			if Dimension(span[s]) eq ComponentDimension(M,s) then
				fulldegrees join:={s};
			end if;
		end for;
		if Verbose then
			PrintAndDelete("Current dimension: "*Sprint(currentdim)*" ("*Sprint(#tospin)*" to spin)        ");
		end if;
	end while;

	if Verbose then
		print "";
	end if;

	return B;

end intrinsic;

//============================================================================
intrinsic Spin(M::ModGr, m::ModGrElt : Verbose:=true) -> List
{Returns the basis of the graded submodule of a graded module M generated by the homogenous element m.}

	return Spin(M,[*m*] : Verbose:=Verbose);

end intrinsic;




//============================================================================
//	Submodules
//============================================================================


//============================================================================
intrinsic RestrictMatrix(A::., U1::ModTupFld, U2::ModTupFld : Rep:= "Parent") -> Mtrx
{
	Restricts the matrix A: V -> W to a matrix B: U1 -> U2, where U1 is a subspace of V and U2 is a subspace of W such that A(U1) subset U2.
}

	image:= [ Vector(Coordinates(U2,b*A)) : b in Basis(U1) ];

	if Rep eq "Sparse" or (Rep eq "Parent" and Type(A) eq MtrxSprs) then
		B:= SparseMatrix(Matrix(image));
	else
		B:= Matrix(image);
	end if;

	return B;

end intrinsic;



//============================================================================
intrinsic Submodule(M::ModGr, U::List : DoSpinning:=true, Rep:="Parent", Verbose:=true)-> ModGr
{Creates the graded submodule of M generated by U.}

	assert ISA(Type(M`BaseRing), Fld);

	// Spin the vectors if U is not a basis.
	if DoSpinning then
		U := Spin(M,U : Verbose:=Verbose);
	end if;

	//determine components of U
	Usupp:= SetToSequence({u`Support[1] : u in U}); // support of U
	Ucomponents:= AssociativeArray(Usupp);		// components of U (we also include the zero-dimensionals as this makes the code below easier)
	for d in Usupp do
		Md := M`ModuleComponents[d];
		Ucomponents[d] := sub<Md | [u`Components[d] : u in U | u`Support[1] eq d]>;
	end for;

	matrices := [**];
	for i:=1 to #M`AlgebraDegrees do
		adeg := M`AlgebraDegrees[i];
		for d in Usupp do
			s := adeg+d;
			if s notin Usupp then
				continue;
			end if;
			m := RestrictMatrix(M`Matrices[i][d], Ucomponents[d], Ucomponents[s] : Rep:=Rep);
			Append(~matrices, <i,d,m>);
		end for;
	end for;

	Ucomponentdims := [<d,Dimension(Ucomponents[d])> : d in Usupp];

	if Rep eq "Parent" then
		Rep := M`Rep;
	end if;

	UM := GradedModule(M`BaseRing, M`AlgebraDegrees, Ucomponentdims, matrices : Rep:=Rep);

	return UM;

end intrinsic;


//============================================================================
intrinsic Submodule(M::ModGr, m::ModGrElt : Rep:="Parent", Verbose:=true)-> ModGr
{}

	return Submodule(M, [*m*] : Rep:=Rep, Verbose:=true, DoSpinning:=true);

end intrinsic;

//============================================================================
//	Quotients
//============================================================================

intrinsic QuotientMatrix(A::., U1::ModTupFld, U2::ModTupFld : Rep:= "Parent") -> Mtrx
{If A defines a morphism f:V1->V2 with f(U1) subset U2, compute the matrix of the map V1/U1->V2/U2}

	V1 := Generic(U1);
	V2 := Generic(U2);
	//Q1,q1 := quo<V1|U1>;
	Q2,q2 := quo<V2|U2>;
	C1 := Complement(V1,U1);
	image := [ q2(v*A) : v in Basis(C1) ];

	if Rep eq "Sparse" or (Rep eq "Parent" and Type(A) eq MtrxSprs) then
		B:= SparseMatrix(Matrix(image));
	else
		B:= Matrix(image);
	end if;

	return B;

end intrinsic;

//============================================================================
intrinsic Quotient(M::ModGr, U::List: DoSpinning:= true, Rep:="Parent") -> ModGr
{Creates the quotient module of M by the submodule generated by the elements in U.}

	assert ISA(Type(M`BaseRing), Fld);

	// Spin the vectors of U to a submodule
	if DoSpinning eq true then
		U:=Spin(M,U);
	end if;

	//determine components of U
	Usupp:= SetToSequence({u`Support[1] : u in U}); // support of U
	Ucomponents:= AssociativeArray(Usupp);		// components of U (we also include the zero-dimensionals as this makes the code below easier)
	for d in M`Support do
		Md := M`ModuleComponents[d];
		Ucomponents[d] := sub<Md | [u`Components[d] : u in U | u`Support[1] eq d]>;
	end for;

	//determine support of Q
	Qsupp := {};
	for d in M`Support do
		Md := M`ModuleComponents[d];
		Ud := Ucomponents[d];
		if Dimension(Ud) ne Dimension(Md) then
			Qsupp join:={d};
		end if;
	end for;

	//set components of Q
	Qcomponents:=AssociativeArray(Qsupp);

	for d in Qsupp do
		Md := M`ModuleComponents[d];
		Ud := Ucomponents[d];
		Qd := quo<Md|Ud>;
		Qcomponents[d] := Qd;
	end for;

	matrices := [**];
	for i:=1 to #M`AlgebraDegrees do
		adeg := M`AlgebraDegrees[i];
		for d in Qsupp do
			s := adeg+d;
			if s notin Qsupp then
				continue;
			end if;
			m := QuotientMatrix(M`Matrices[i][d], Ucomponents[d], Ucomponents[s] : Rep:=Rep);
			Append(~matrices, <i,d,m>);
		end for;
	end for;

	Qcomponentdims := [<d,Dimension(Qcomponents[d])> : d in Qsupp];

	if Rep eq "Parent" then
		Rep := M`Rep;
	end if;

	Q := GradedModule(M`BaseRing, M`AlgebraDegrees, Qcomponentdims, matrices : Rep:=Rep);

	return Q;

end intrinsic;

//============================================================================
intrinsic ActionGenerators(M::ModGr : Rep:="Parent") -> SeqEnum
{
	The RModule defined by M.
}

	//We have to be a bit careful since in ModGr we do not save the action matrices which are zero. For RModule we need all of them.

	dim := M`Dimension;

	if Rep eq "Parent" then
		Rep := M`Rep;
	end if;
	if Rep eq "Dense" then
		Q := [ ZeroMatrix(M`BaseRing, dim, dim) : i in [1..#M`AlgebraDegrees] ];
	else
		Q := [ SparseMatrix(M`BaseRing, dim, dim) : i in [1..#M`AlgebraDegrees] ];
	end if;

	V := VectorSpace(M`BaseRing, dim);

	for i:=1 to #M`AlgebraDegrees do
		for j:=1 to #M`Support do
			sourcedegree := M`Support[j];
			sourcedim := M`ComponentDimension[j]; //dim M_sourcedegree
			targetdegree := sourcedegree+M`AlgebraDegrees[i];
			p := Position(M`Support,targetdegree);
			sourceoffset := ArraySum([ M`ComponentDimension[l] : l in [1..j-1] ] : ZeroElement:=Zero(Integers()));  // basis vectors of M_d then have numbers offset+d..offset+sourcedim
			targetoffset := ArraySum([ M`ComponentDimension[l] : l in [1..p-1] ] : ZeroElement:=Zero(Integers())); // basis vectors of M_{imagedegree} then have numbers offset+d..offset+dim(M_{imagedegree})
			for k:=1 to sourcedim do
				v := Zero(V); //will be image of action of A.i on k-th basis vector of M_d
				if IsDefined(M`Matrices[i],sourcedegree) then
					w := M`Matrices[i][sourcedegree][k];
					for l in Support(w) do
						v[targetoffset+l] := w[l];
					end for;
				end if;
				if Rep eq "Dense" then
					Q[i][sourceoffset+k] := v;
				else
					for j in Support(v) do
						Q[i][sourceoffset+k][j] := v[j];
					end for;
				end if;
			end for;
		end for;
	end for;

	return Q;

end intrinsic;

intrinsic RModule(M::ModGr) -> ModRng
{}

	return RModule(ActionGenerators(M) : Rep:="Dense");

end intrinsic;

//==============================================================================
intrinsic NumberOfNonZeroEntries(M::ModGr) -> RngIntElt
{}

	N := 0;
	for i:=1 to #M`AlgebraDegrees do
		for d in M`Support do
			if IsDefined(M`Matrices[i],d) then
				N +:= NumberOfNonZeroEntries(M`Matrices[i][d]);
			end if;
		end for;
	end for;

	return N;

end intrinsic;

//==============================================================================
intrinsic ChangeRing(M::ModGr, theta::Map : Rep:="Parent") -> ModGr
{}

	if Rep eq "Parent" then
		Rep := M`Rep;
	end if;

	Msupp := Support(M);
	Mcomponentdims := [ <d,Dimension(M`ModuleComponents[d])> : d in Msupp ];
	matrices := [* *];
	for i:=1 to #M`AlgebraDegrees do
		for d in Msupp do
			if IsDefined(M`Matrices[i], d) then
				m := ChangeRing(M`Matrices[i][d],theta);
				if Rep eq "Sparse" and Type(m) ne MtrxSprs then
					m := SparseMatrix(m);
				elif Rep eq "Dense" and Type(m) ne Mtrx then
					m := Matrix(m);
				end if;
				Append(~matrices, <i,d,m>);
			end if;
		end for;
	end for;

	return GradedModule(Codomain(theta), M`AlgebraDegrees, Mcomponentdims, matrices : Rep:=Rep);

end intrinsic;

//==============================================================================
intrinsic GradedModuleOld(M::ModGr : Rep:="Parent") -> ModGrOld
{}

	if Rep eq "Parent" then
		Rep := M`Rep;
	end if;

	Q := ActionGenerators(M : Rep:=Rep);
	Qdegs := [ Degree(M.i) : i in [1..M`Dimension] ];
	Algdegs := M`AlgebraDegrees;

	return GradedModuleOld(Q, Qdegs, Algdegs);

end intrinsic;
