/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Compute a presentation of the center of the rational Cherednik algebra at t=0.

	Joint work with Cédric Bonnafé (Montpellier).
*/

declare attributes AlgChe:
	CenterSpace, //The polynomial ring R[z_1,...,z_N] where R is the base ring of H and the z_i are the center generators of H. The Calogero-Moser space lives in the associated affine space.
	CenterPresentation;

//============================================================================
intrinsic ChangeRing(f::RngMPolElt, phi::Map) -> RngMPolElt
{Change the coefficient ring of the polynomial ring of f along phi.}

	P := PolynomialRing(Codomain(phi), Ngens(Parent(f)));
	AssignNames(~P, Names(Parent(f)));
	P := ChangeOrder(P, MonomialOrder(Parent(f)));
	g := Zero(P);
	for mon in Monomials(f) do
		coeff := MonomialCoefficient(f, mon);
		exp := Exponents(mon);
		g +:= phi(coeff)*Monomial(P,exp);
	end for;

	return g;

end intrinsic;

//============================================================================
intrinsic ChangeRing(I::RngMPol, phi::Map) -> RngMPol
{Change the coefficient ring of the polynomial ring of I along phi.}

	P := PolynomialRing(Codomain(phi), Ngens(Parent(Basis(I)[1])));
	P := ChangeOrder(P, MonomialOrder(Parent(Basis(I)[1])));
	AssignNames(~P, Names(Parent(Basis(I)[1])));
	return ideal<P|[ P!ChangeRing(f,phi) : f in Basis(I) ]>;

end intrinsic;

//============================================================================
intrinsic ChangeRing(M::AlgMatElt[RngMPol], phi::Map) -> AlgMatElt
{Change the coefficient ring of the polynomial ring of the base ring of M along phi.}

	//get the codomain of ring change
	R := Parent(ChangeRing(Zero(BaseRing(M)), phi));
	Mnew := ZeroMatrix(R,Nrows(M),Ncols(M));
	for i:=1 to Nrows(M) do
		for j:=1 to Ncols(M) do
			Mnew[i,j] := R!ChangeRing(M[i,j],phi);
		end for;
	end for;
	return Mnew;

end intrinsic;

//============================================================================
intrinsic CenterSpace(~H::AlgChe)
{
	The polynomial ring R[z_1,...,z_N] where R is the base ring of H and the z_i are the center generators of H. The Calogero-Moser space lives in the associated affine space.
}
	if assigned H`CenterSpace then
		return;
	end if;

	W := H`Group;
	SymplecticDoublingFundamentalInvariants(~W);
	N := #W`SymplecticDoublingFundamentalInvariants;
	Q := PolynomialRing(H`BaseRing, N);
	AssignNames(~Q, ["z"*Sprint(i) : i in [1..N]]);
	H`CenterSpace := Q;

end intrinsic;

//============================================================================
intrinsic Preimage(H::AlgChe, z::AlgCheElt) -> RngMPolElt
{
	The preimage of a central element z in a rational Cherednik algebra. Only works if the parameters of H can be specialized to 0, i.e., it must be possible to specialize to the undeformed situation.
}

	CenterGenerators(~H);
	CenterSpace(~H);
	W := H`Group;
	P := Parent(W`SymplecticDoublingFundamentalInvariants[1]); //K[V^* + V]

	//create map to kill the parameters c_s to get to undeformed situation
	if Type(H`BaseRing) ne RngMPol and Type(H`BaseRing) ne RngUPol then
		error "Need parameters which may be specialized to 0";
	end if;

	psi := hom<H`BaseRing -> BaseRing(W) | [0 : i in [1..Rank(H`BaseRing)]]>;

	for i:=1 to #W`ReflectionClasses do
		if psi(H`cParameter(i)) ne 0 then
			error "Need parameters which may be specialized to 0";
		end if;
	end for;

	res := Zero(H`CenterSpace);
	if z eq Zero(H) then
		return res;
	else
		coeff := Coefficient(z`Element, Identity(W));

		k := Minimum({Valuation(c) : c in Coefficients(coeff) | c ne 0});
		coeffk := &+[ Projection(MonomialCoefficient(coeff,mon),k)*mon : mon in Monomials(coeff) ];
		znew := z;
		for monc in MonomialsOfDegree(H`BaseRing, k) do
			hmonc := P!ChangeRing(&+[ MonomialCoefficient(MonomialCoefficient(coeffk,mon),monc)*mon : mon in Monomials(coeffk)], psi); //kill variables
			r,fmonc:=HomogeneousModuleTest(W`SymplecticDoublingFundamentalInvariants,[P!1],hmonc);

			res +:= monc*H`CenterSpace!fmonc[1];

			znew := znew - monc*Evaluate(fmonc[1], H`CenterGenerators);

		end for;

		return res + Preimage(H,znew);
	end if;

end intrinsic;

//============================================================================
intrinsic Evaluate(f::RngMPolElt, elts::List) -> AlgCheElt
{Evaluate the multivariate polynomial in a list of ring elements.}

	H := elts[1]`Parent;
	h := Zero(H);
	for mon in Monomials(f) do
		coeff := MonomialCoefficient(f,mon);
		exp := Exponents(mon);
		prod := One(H);
		for i:=1 to #exp do
			prod := prod*elts[i]^exp[i];
		end for;
		prod := prod*coeff;
		h := h + prod;
	end for;

	return h;

end intrinsic;

//============================================================================
intrinsic Valuation(f::RngMPolElt) -> RngIntElt
{The valuation of f. If f is zero, then -1 is returned.}

	if f eq 0 then
		return -1;
	else
		return Minimum({Degree(mon) : mon in Monomials(f)});
	end if;

end intrinsic;

//============================================================================
intrinsic Projection(f::RngMPolElt, k::RngIntElt) -> RngMPolElt
{The projection of f to the homogeneous degree k polynomials.}

	return ArraySum([ MonomialCoefficient(f,mon)*mon : mon in Monomials(f) | Degree(mon) eq k] : ZeroElement:=Zero(Parent(f)));

end intrinsic;

//============================================================================
intrinsic Projection(f::RngUPolElt, k::RngIntElt) -> RngUPolElt
{The projection of f to the homogeneous degree k polynomials.}

	return ArraySum([ MonomialCoefficient(f,mon)*mon : mon in Monomials(f) | Degree(mon) eq k] : ZeroElement:=Zero(Parent(f)));

end intrinsic;

//===========================================================================
intrinsic MonomialsOfDegree(R::RngUPol, k::RngIntElt) -> SeqEnum
{}

	return [R.1^k];

end intrinsic;

//===========================================================================
intrinsic Presentation(R::RngInvar) -> SeqEnum
{A presentation of the invariant ring R (the relations are with respect to the
fundamental invariants of R).}

	fund := FundamentalInvariants(R);
	prim := PrimaryInvariants(R);
	sec := IrreducibleSecondaryInvariants(R);
	invar := prim cat sec;
	P := PolynomialRing(BaseRing(R), #fund);
	A := Algebra(R);

	invarpres := [];
	for f in invar do
		b,g := HomogeneousModuleTest(fund,[R!1],f);
		Append(~invarpres, g[1]);
	end for;

	rel := RelationIdeal(R);

	phi:=hom<A->P|invarpres>;

	return ideal<P|[phi(r) : r in Basis(rel)]>;

end intrinsic;

//============================================================================
intrinsic CenterPresentation(~H::AlgChe : Weights:=false, SaveToDB := false, UseDB:=true)
{The presentation of the center of the generic rational Cherednik algebra H.}

	if assigned H`CenterPresentation then
		return;
	end if;

	W := H`Group;
	SymplecticDoublingFundamentalInvariants(~W);
	R := W`SymplecticDoubling`InvariantRing;
	print "Computing presentation of invariant ring.";
	pres := Presentation(R);
	CenterGenerators(~H);
	rels := [];
	P := Parent(Basis(pres)[1]);
	CenterSpace(~H);
	phi:=hom<P->H`CenterSpace | [H`CenterSpace.i : i in [1..Ngens(H`CenterSpace)] ]>;
	count:=0;
	for f in Basis(pres) do
		count +:= 1;
		print "Deforming relation "*Sprint(count)*" of "*Sprint(#Basis(pres));
		z := Evaluate(f,H`CenterGenerators);
		print "Evaluation done. Computing preimage.";
		zpre := Preimage(H,z);
		rel := phi(f)-zpre;
		Append(~rels, rel);
	end for;

	H`CenterPresentation := rels;

	if SaveToDB then
		K := BaseRing(W);
		R := H`BaseRing;
		P := H`CenterSpace;

		str := "K := ";
		fieldstr := "";
		if Type(K) eq FldRat then
			fieldstr := "RationalField()";
			str *:= fieldstr*";\n";
		elif Type(K) eq FldCyc then
			fieldstr := "CyclotomicField("*Sprint(CyclotomicOrder(K))*")";
			str *:= fieldstr*";\n";
			str *:= Sprint(K.1)*" := RootOfUnity("*Sprint(CyclotomicOrder(K))*");\n";
		else
			error "Not yet implemented for this type of base ring.";
		end if;

		if Type(R) eq RngMPol then
			if BaseRing(R) ne K then
				error "Not yet implemented for this type of base ring.";
			end if;
			rstr := "PolynomialRing("*fieldstr*", "*Sprint(Ngens(R))*")";
			str *:= "R := "*rstr*";\n";
			for i:=1 to Ngens(R) do
				str *:= Sprint(R.i)*" := R."*Sprint(i)*";\n";
			end for;
		elif Type(R) eq Fld then
			if R ne K then
				error "Not yet implemented for this type of base ring.";
			end if;
			rstr := fieldstr;
			str *:= "R := "*rstr*";\n";
		end if;

		str *:= "P := PolynomialRing(R,"*Sprint(Ngens(P))*");\n";
		for i:=1 to Ngens(P) do
			str *:= Sprint(P.i)*" := P."*Sprint(i)*";\n";
		end for;

		for i:=1 to #rels do
			str *:= "rel"*Sprint(i)*" := "*Sprint(rels[i])*";\n";
		end for;

		str *:= "\n";
		str *:= "return [";
		for i:=1 to #rels do
			str *:= "f"*Sprint(i);
			if i lt #rels then
				str *:= ",";
			end if;
		end for;
		str *:= "]";

		CHAMP_SaveToDB(str, H`DBDir, "CenterPresentation");

	end if;

end intrinsic;

intrinsic CenterPresentation(H::AlgChe : Weights:=false) -> SeqEnum
{}

	CenterPresentation(~H : Weights:=Weights);
	return H`CenterPresentation;

end intrinsic;

//============================================================================
intrinsic FixedPoints(H::AlgChe, d::RngIntElt) -> SeqEnum
{The mu_d fixed points of the center of the generic rational Cherednik algebra H.}

	W := H`Group;
	pres := CenterPresentation(H);
	funddegs := [ -Bidegree(f)[1]+Bidegree(f)[2] : f in W`SymplecticDoublingFundamentalInvariants ];
	posd := [ i : i in [1..#funddegs] | IsDivisibleBy(funddegs[i], d) ];
	A := Parent(pres[1]);
	P := PolynomialRing(H`BaseRing, #posd);
	morph := [ Zero(P) : i in [1..#funddegs] ];
	for j:=1 to #posd do
		morph[posd[j]] := P.j;
	end for;

	phi := hom<A->P | morph >;

	return [phi(f) : f in pres];

end intrinsic;

//============================================================================
intrinsic Homogenize(I::RngMPol) -> RngMPol
{The homogenization of I.}

	P := Parent(Basis(I)[1]);
	Q := PolynomialRing(BaseRing(P), Ngens(P)+1);
	AssignNames(~Q, Names(P) cat ["t"]);
	t := Q.Ngens(Q);
	phi := hom<P->Q | [Q.i : i in [1..Ngens(P)]]>;

	J := [];
	for f in Basis(I) do
		fhomog := Zero(Q);
		d := Degree(f);
		for mon in Monomials(f) do
			fhomog +:= fhomog+MonomialCoefficient(f,mon)*t^(d-Degree(mon))*phi(mon);
		end for;
		Append(~J, fhomog);
	end for;

	return ideal<Q|J>;

end intrinsic;

//============================================================================
intrinsic Scheme(I::RngMPol) -> Sch
{The subscheme of the parent polynomial ring of I defined by I.}

	P := Parent(Basis(I)[1]);
	A := AffineSpace(P);
	return Scheme(A,I);

end intrinsic;

//============================================================================
intrinsic Specialize(I::RngMPol, X::SeqEnum) -> RngMPol
{If I is an ideal of a polynomial over a base ring involving base rings, then create the ideal obtained by specializing the parameters of the base ring to X.}

	P := Parent(Basis(I)[1]);
	k := BaseRing(P);
	l := BaseRing(k);
	phi := hom<k->l | X >;
	return ChangeRing(I,phi);

end intrinsic;

//============================================================================
intrinsic BettiNumbers(S::Sch) -> SeqEnum
{}

	I := Ideal(S);
	if not IsHomogeneous(I) then
		error "Scheme has to be defined by homogenous ideal.";
	end if;
	return BettiNumbers(GradedModule(I));

end intrinsic;
