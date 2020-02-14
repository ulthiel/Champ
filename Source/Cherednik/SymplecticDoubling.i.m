/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Intrinsics for symplectic doublings of matrix groups and their invariant rings.
*/

declare attributes GrpMat:
    SymplecticDoubling, // The symplectic doubling as a matrix group.
    SymplecticDoublingCoordinateAlgebra, //The coordinate algebra K[V^*+V] of the symplectic doubling (note the order V^*+V !), so first y, then x.
    SymplecticDoublingPrimaryInvariants, //Primary invariants of the symplectic doubling.
	SymplecticDoublingSecondaryInvariants, //Secondary invariants of the symplectic doubling.
    SymplecticDoublingFundamentalInvariants,
    SympecticDoublingCoordinateAlgebraHilbertSeries,
    CoordinateAlgebraToSymplecticDoublingCoordinateAlgebraXpart, //K[V] to K[V] \subset K[V^*+V]
	CoordinateAlgebraToSymplecticDoublingCoordinateAlgebraYpart; //K[V] to K[V^*] \subset K[V^*+V]

//the following additional attributes help for computations (we assume that these have been assigned)
declare attributes RngMPol:
	xPart, //we assume that P=K[y_1,...,y_n,x_1,...,x_n], then xPart is K[x_1,...,x_n]
	yPart, // this is K[y_1,...,y_n] ; this distinction is not really necessary here, in fact we take yPart just to be xPart but it becomes necessary for quotients.
	xEmbedding, //the embedding xPart -> K[x_1,...,x_n] \subset P
	yEmbedding;	//the embedding yPart -> K[y_1,...,y_n] \subset P

declare attributes RngMPolRes:
	xPart, //we assume that P=K[y_1,...,y_n,x_1,...,x_n], then xPart is K[x_1,...,x_n]
	yPart, // this is K[y_1,...,y_n] ; this distinction is not really necessary here, in fact we take yPart just to be xPart but it becomes necessary for quotients.
	xEmbedding, //the embedding xPart -> K[x_1,...,x_n] \subset P
	yEmbedding;	//the embedding yPart -> K[y_1,...,y_n] \subset P

//===========================================================================
intrinsic SymplecticDoubling(g::GrpMatElt) -> AlgMatElt
{
	If g acts on the vector space V, this returns the matrix describing the action of g on V \oplus V^*. This is simply the block diagonal matrix with g in the (1,1)-block and the transpose-inverse of g in the (2,2)-block.
}

	return BlockMatrix( 2, 2, [Matrix(g), ZeroMatrix(BaseRing(g), Nrows(g), Ncols(g)), ZeroMatrix(BaseRing(g), Nrows(g), Ncols(g)), Transpose(Matrix(g))^-1]);

end intrinsic;

//===========================================================================
intrinsic SymplecticDoubling(~G::GrpMat)
{
	If G acts on V, attaches/returns the matrix group describing the action of G on V \oplus V^*.
}

    if assigned G`SymplecticDoubling then
        return;
    end if;

    DualGroup(~G);

    mats := [ SymplecticDoubling(G.i) : i in [1..Ngens(G)] ];

    G`SymplecticDoubling := MatrixGroup<2*Dimension(G), BaseRing(G)|mats>;

    G`SymplecticDoubling`InvariantRing := InvariantRing(G`SymplecticDoubling); //this is correct; no dual group here!

    G`SymplecticDoublingCoordinateAlgebra := PolynomialRing(G`SymplecticDoubling`InvariantRing);
    names := [ "y"*Sprint(i) : i in [1..Dimension(G)]] cat [ "x"*Sprint(i) : i in [1..Dimension(G)]];
    AssignNames(~G`SymplecticDoublingCoordinateAlgebra, names);

    Generators(~G`SymplecticDoublingCoordinateAlgebra);
    G`SymplecticDoublingCoordinateAlgebra`One := One(G`SymplecticDoublingCoordinateAlgebra);
	G`SymplecticDoublingCoordinateAlgebra`Zero := Zero(G`SymplecticDoublingCoordinateAlgebra);

    CoordinateAlgebra(~G);

    G`Dimension := Dimension(G);
    n := G`Dimension;

	G`SymplecticDoublingCoordinateAlgebra`xPart := G`CoordinateAlgebra;
	G`SymplecticDoublingCoordinateAlgebra`yPart := G`DualCoordinateAlgebra;
	G`SymplecticDoublingCoordinateAlgebra`xEmbedding :=
	hom<G`CoordinateAlgebra -> G`SymplecticDoublingCoordinateAlgebra | [G`SymplecticDoublingCoordinateAlgebra`Generators[i] : i in [n+1..2*n]]>;
	G`SymplecticDoublingCoordinateAlgebra`yEmbedding := hom<G`DualCoordinateAlgebra -> G`SymplecticDoublingCoordinateAlgebra | [G`SymplecticDoublingCoordinateAlgebra`Generators[i] : i in [1..n]]>;

end intrinsic;

//===========================================================================
intrinsic SymplecticDoubling(G::GrpMat) -> GrpMat
{}

    SymplecticDoubling(~G);
    return G`SymplecticDoubling;

end intrinsic;


//===========================================================================
intrinsic SymplecticDoublingPrimaryInvariants(~G::GrpMat)
{Attaches/computes primary invariants of the symplectic doubling.}

    if assigned G`SymplecticDoublingPrimaryInvariants then
        return;
    end if;

    SymplecticDoubling(~G);
    G`SymplecticDoublingPrimaryInvariants := PrimaryInvariants(G`SymplecticDoubling`InvariantRing);

end intrinsic;

//===========================================================================
intrinsic SymplecticDoublingPrimaryInvariants(G::GrpMat) -> SeqEnum
{}

    SymplecticDoublingPrimaryInvariants(~G);
    return G`SymplecticDoublingPrimaryInvariants;

end intrinsic;

//===========================================================================
intrinsic SymplecticDoublingSecondaryInvariants(~G::GrpMat)
{}

    if assigned G`SymplecticDoublingSecondaryInvariants then
        return;
    end if;

    SymplecticDoubling(~G);
    G`SymplecticDoublingSecondaryInvariants := SecondaryInvariants(G`SymplecticDoubling`InvariantRing);

end intrinsic;

//===========================================================================
intrinsic SymplecticDoublingSecondaryInvariants(G::GrpMat) -> SeqEnum
{}

    SymplecticDoublingSecondaryInvariants(~G);
    return G`SymplecticDoublingSecondaryInvariants;

end intrinsic;

//===========================================================================
intrinsic SymplecticDoublingFundamentalInvariants(~G::GrpMat : UseDB:=true, SaveToDB:=false)
{
	Attaches/computes fundamental invariants of the symplectic doubling.
}

    if assigned G`SymplecticDoublingFundamentalInvariants and not SaveToDB then
        return;
    end if;

	SymplecticDoubling(~G);

	if not assigned G`SymplecticDoublingFundamentalInvariants then
		if not SaveToDB and UseDB and assigned G`DBDir and CHAMP_ExistsInDB(G`DBDir, "Invariants/SymplecticDoublingFundamentalInvariants") then
			G`SymplecticDoublingFundamentalInvariants := CHAMP_GetFromDB(G`DBDir, "Invariants/SymplecticDoublingFundamentalInvariants");
			print "Fetched from DB";
		else
    		G`SymplecticDoublingFundamentalInvariants := FundamentalInvariants(G`SymplecticDoubling`InvariantRing);
    	end if;
    end if;

    if SaveToDB then
    	if not assigned G`DBDir then
    		error "No database directory assigned to group.";
    	end if;
    	str := "/*\n\tSymplectic doubling fundamental invariants\n\tVersion: "*CHAMP_GetVersion()*"\n\tDate: "*Date()*"\n*/\n";
		str *:= "P<";
		for j:=1 to Dimension(G) do
			str *:= "y"*Sprint(j);
			if j lt Dimension(G) then
				str *:= ",";
			end if;
		end for;
		str *:= ",";
		for j:=1 to Dimension(G) do
			str *:= "x"*Sprint(j);
			if j lt Dimension(G) then
				str *:= ",";
			end if;
		end for;

    	str *:= "> := PolynomialRing(";
    	K := BaseRing(G);
    	if Type(K) eq FldRat then
    		str *:= "Rationals(),";
    	elif Type(K) eq FldCyc then
    		str *:= "CyclotomicField("*Sprint(CyclotomicOrder(K))*"),";
    	else
    		error "No routine for this base ring yet.";
    	end if;
    	str *:= Sprint(2*Dimension(G))*");\n";
    	if Type(K) eq FldCyc then
    		str *:= Sprint(K.1)*" := RootOfUnity("*Sprint(CyclotomicOrder(K))*");\n";
    	end if;
    	for i:=1 to #G`SymplecticDoublingFundamentalInvariants do
    		str *:= "f"*Sprint(i)*" := "*Sprint(G`SymplecticDoublingFundamentalInvariants[i])*";\n";
    	end for;
    	str *:= "return [";
    	for i:=1 to #G`SymplecticDoublingFundamentalInvariants do
    		str *:= "f"*Sprint(i);
    		if i lt #G`SymplecticDoublingFundamentalInvariants then
    			str *:= ",";
    		end if;
    	end for;
    	str *:= "]";
    	CHAMP_SaveToDB(str, G`DBDir, "SymplecticDoublingFundamentalInvariants");
    	print "Saved";
    end if;

end intrinsic;

//===========================================================================
intrinsic SymplecticDoublingFundamentalInvariants(G::GrpMat) -> SeqEnum
{}

    SymplecticDoublingFundamentalInvariants(~G);
    return G`SymplecticDoublingFundamentalInvariants;

end intrinsic;

//===========================================================================
intrinsic Bidegree(f::RngMPolElt) -> Tup
{
	Suppose f is an element of a polynomial ring R of even rank r. We can refine the natural N-grading on R into an (N,N)-grading, where the first component is the degree of the first r/2 variables and the second component is the degree of the last r/2 variables. We assume here that f is bi-homogeneous.
}

    mon := Monomials(f)[1];
    N := Rank(Parent(f));
    if N mod 2 ne 0 then
        error "Polynomial ring has to be of even rank.";
    end if;
    n := N div 2;
    exp := Exponents(mon);
    return <&+[exp[i] : i in [1..n]], &+[exp[i] : i in [n+1..N]]>;

end intrinsic;

//===========================================================================
intrinsic Bidegree(f::RngMPolResElt) -> Tup
{

}

    mon := Monomials(f)[1];
    N := Rank(Parent(f));
    if N mod 2 ne 0 then
        error "Polynomial ring has to be of even rank.";
    end if;
    n := N div 2;
    exp := Exponents(mon);
    return <&+[exp[i] : i in [1..n]], &+[exp[i] : i in [n+1..N]]>;

end intrinsic;


//=========================================================================
intrinsic SymplecticDoublingAction(f::RngMPolElt, g::GrpMatElt) -> RngMPolElt
{
	Action of g on an element f of the symplectic doubling coordinate algebra K[V^* + V]. Note the order V^* + V.
}

	G := Parent(g);
	n := Ncols(g);
	P := Parent(f);

	//for some reason the following method is the quickest, I don't know why
	//assume that P`xPart, P`yPart, P`xEmbedding, P`yEmbedding are properly assigned (see beginning of file)
	res := P`Zero;
	for mon in Monomials(f) do
		coeff := MonomialCoefficient(f, mon);
		exp := Exponents(mon);
		ypart := ArrayProduct([(P`yPart`Generators[i]^exp[i])^g : i in [1..n] | exp[i] ne 0] : OneElement:=P`yPart`One);
		xpart := ArrayProduct([(P`xPart`Generators[i-n]^exp[i])^Transpose(g^-1) : i in [n+1..2*n] | exp[i] ne 0] : OneElement:=P`xPart`One);
		mong := P`yEmbedding(ypart)*P`xEmbedding(xpart);
		res +:= coeff*mong;
	end for;

	return res;

	//split method
/*
	C := P`xPart; //see beginning of file for description
	res := Zero(P);
	for mon in Monomials(f) do
		exp := Exponents(mon);
		coeff := MonomialCoefficient(f,mon);
		ypart := Monomial(C, exp[1..n]);
		xpart := Monomial(C, exp[n+1..2*n]);
		ypartg := Action(ypart, g : Dual:=false);
		xpartg := Action(xpart, g : Dual:=true);
		ypartg := P`yEmbedding(ypartg);
		xpartg := P`xEmbedding(xpartg);
		mong := ypartg*xpartg;
		res +:= coeff*mong;
	end for;

	return res;
*/

	//direct method
/*
	return f^SymplecticDoubling(g);
*/

end intrinsic;

//=========================================================================
intrinsic SymplecticDoublingAction(f::RngMPolResElt, g::GrpMatElt) -> RngMPolElt
{
	Symplectic doubling action on a quotient of the symplectic doubling coordinate algebra (by a non-mixed ideal).
}

	G := Parent(g);
	n := Ncols(g);
	P := Parent(f);
	res := Zero(P);
	for mon in Monomials(f) do
		exp := Exponents(mon);
		coeff := MonomialCoefficient(f,mon);
		ypart := Monomial(P`yPart, exp[1..n]);
		xpart := Monomial(P`xPart, exp[n+1..2*n]);
		ypartg := Action(ypart, g : Dual:=false);
		xpartg := Action(xpart, g : Dual:=true);
		ypartg := P`yEmbedding(ypartg);
		xpartg := P`xEmbedding(xpartg);
		mong := ypartg*xpartg;
		res +:= coeff*mong;
	end for;

	return res;

end intrinsic;

//============================================================================
intrinsic SympecticDoublingCoordinateAlgebraHilbertSeries(~G::GrpMat)
{}

	if assigned G`SympecticDoublingCoordinateAlgebraHilbertSeries then
        return;
    end if;

   	Classes(~G);
   	P<t,u> := RationalFunctionField(BaseRing(G),2);
   	ser := Zero(P);
   	I := IdentityMatrix(P, Dimension(G));

   	for i:=1 to #G`Classes do
   		g := ClassRepresentative(G,i);
   		len := G`Classes[i][2];
   		ser +:= len*1/(Determinant(I-ChangeRing(g, P)*t)*Determinant(I-ChangeRing(g^-1,P)*u	));
	end for;

	ser *:= 1/Order(G);

	G`SympecticDoublingCoordinateAlgebraHilbertSeries := ser;

end intrinsic;
