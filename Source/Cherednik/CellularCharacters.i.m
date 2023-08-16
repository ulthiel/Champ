/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Joint work with Cédric Bonnafé (Montpellier).
*/

//============================================================================
intrinsic CommonBaseRingForCherednikStuff(G::GrpMat, c::Map, reps::List) -> Rng
{A common base ring to do Cherednik type computations with G, c, and specified representations.}

	L := BaseRing(G);

	for i:=1 to #reps do
		L := CommonOverfield(L, BaseRing(Codomain(reps[i])));
	end for;

	if Type(Codomain(c)) eq FldRat or Type(Codomain(c)) eq FldCyc then
		L := CommonOverfield(L, Codomain(c));
	elif Type(Codomain(c)) eq RngMPol or Type(Codomain(c)) eq FldFunRat then
		L := ChangeRing(Codomain(c), L);
	end if;

	return L;

end intrinsic;

//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map, y::ModTupFldElt) -> AlgMatElt
{The action of the Gaudin operator at y on the full group algebra.}

	L := CommonBaseRingForCherednikStuff(G,c,[**]);

    NumberingMap(~G);
    ReflectionLibrary(~G);
    K := RationalFunctionField(Codomain(c), Dimension(G)); //k(V), k=Codomain(c)
    AssignNames(~K, ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
    D:=ZeroMatrix(K, #G, #G);
    E:=VectorSpace(K, #G);
    for i:=1 to #G do
    	w := G`InverseNumberingMap(i);
        //image of e_w under D^{c,v,vstar}_y
        img := Zero(E);
        for s in G`ReflectionLibraryFlat do
        	coroot := &+[s`Coroot[i]*K.i : i in [1..Ngens(K)]];
            img +:= s`Eigenvalue/(s`Eigenvalue-1)*c(s`ReflectionClass)*CanonicalPairing(y,s`Coroot)/coroot*E.G`NumberingMap(s`Element*w);
        end for;
        //print img;
        D[i]:=img;
    end for;

    return D;

end intrinsic;

//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map) -> AlgMatElt
{The full Gaudin operator.}

	L := CommonBaseRingForCherednikStuff(G,c,[**]);

	K := RationalFunctionField(L, Dimension(G));
	S := RationalFunctionField(L, 2*Dimension(G));
	phi := hom<K->S | [S.(i+Dimension(G)) : i in [1..Dimension(G)]]>;
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]] cat ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(ChangeRing(GaudinOperator(G,c,V.i),K),phi) : i in [1..Dimension(G)]];

end intrinsic;


//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map, y::ModTupFldElt, rho::Map) -> AlgMatElt
{The action of the Gaudin operator at y on rho.}

	L := CommonBaseRingForCherednikStuff(G,c,[*rho*]);

    NumberingMap(~G);
    ReflectionLibrary(~G);
    K := RationalFunctionField(L, Dimension(G)); //k(V), k=Codomain(c)
    AssignNames(~K, ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
    D:=ZeroMatrix(K, Dimension(Codomain(rho)), Dimension(Codomain(rho)));
    E:=VectorSpace(K, Dimension(Codomain(rho)));

    D := ZeroMatrix(K, Degree(Codomain(rho)), Degree(Codomain(rho)));
    for s in G`ReflectionLibraryFlat do
     	coroot := &+[s`Coroot[i]*K.i : i in [1..Ngens(K)]];
        D +:=  (K!s`Eigenvalue/(K!s`Eigenvalue-1))*(K!c(s`ReflectionClass))*CanonicalPairing(y,s`Coroot)/coroot*ChangeRing(rho(s`Element), K);
    end for;

    return D;

end intrinsic;


//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map, rho::Map) -> AlgMatElt
{The Gaudin operator acting on rho.}

	L := CommonBaseRingForCherednikStuff(G,c,[*rho*]);

	K := RationalFunctionField(L, Dimension(G));
	S := RationalFunctionField(L, 2*Dimension(G));
	phi := hom<K->S | [S.(i+Dimension(G)) : i in [1..Dimension(G)]]>;
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]] cat ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(ChangeRing(GaudinOperator(G,c,V.i,rho),K),phi) : i in [1..Dimension(G)]];

end intrinsic;

//============================================================================
intrinsic GaudinOperators(G::GrpMat, c::Map, fam::SetIndx) -> List
{}

	Representations(~G);

	L := CommonBaseRingForCherednikStuff(G,c,[*G`Representations[0][i] : i in fam*]);

	gaudins := [**];

	for i in fam do
		print i;
		rho := G`Representations[0][i];
		D := GaudinOperator(G,c,rho);
		R := BaseRing(D);
		S := ChangeRing(R, L);
		D := ChangeRing(D, S);
		Append(~gaudins, D);
	end for;

	return gaudins;

end intrinsic;

//============================================================================
intrinsic GaudinOperatorSpecialized(G::GrpMat, c::Map, y::ModTupFldElt, vreg::ModTupFldElt, rho::Map) -> AlgMatElt
{The action of the Gaudin operator at y specialized in the regular vector rho on rho.}

	K := CommonBaseRingForCherednikStuff(G,c,[*rho*]);

    NumberingMap(~G);
    ReflectionLibrary(~G);

    D:=ZeroMatrix(K, Degree(Codomain(rho)), Degree(Codomain(rho)));
    for s in G`ReflectionLibraryFlat do
     	//coroot := &+[s`Coroot[i]*K.i : i in [1..Ngens(K)]];
        D +:=  (K!s`Eigenvalue/(K!s`Eigenvalue-1))*(K!c(s`ReflectionClass))*CanonicalPairing(y,s`Coroot)/CanonicalPairing(vreg,s`Coroot)*ChangeRing(rho(s`Element), K);
    end for;

    return D;

end intrinsic;

//============================================================================
intrinsic GaudinOperatorSpecialized(G::GrpMat, c::Map, vreg::ModTupFldElt, rho::Map) -> AlgMatElt
{The Gaudin operator specialized in the regular vector vreg acting on rho.}
	L := CommonBaseRingForCherednikStuff(G,c,[*rho*]);

	S := PolynomialRing(L, Dimension(G) : Global:=true);
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(GaudinOperatorSpecialized(G,c,V.i,vreg,rho),L) : i in [1..Dimension(G)]];

end intrinsic;

//============================================================================
intrinsic GaudinOperatorsSpecialized(G::GrpMat, c::Map, vreg::ModTupFldElt, fam::SetIndx) -> AlgMatElt
{The Gaudin operator specialized in the regular vector vreg acting on rho.}

	Representations(~G);

	L := CommonBaseRingForCherednikStuff(G,c,[*G`Representations[0][i] : i in fam*]);

	gaudins := [**];

	for i in fam do
		print i;
		rho := G`Representations[0][i];
		D := GaudinOperatorSpecialized(G,c,vreg,rho);
		R := BaseRing(D);
		S := ChangeRing(R, L);
		D := ChangeRing(D, S);
		Append(~gaudins, D);
	end for;

	return gaudins;

end intrinsic;


intrinsic DualRepresentation(f::Map) -> Map
{}

	G := Domain(f);
	codom := Codomain(f);
	K := BaseRing(codom);
	n := Dimension(codom);
	gl := GL(n, K);

	mats := [ Transpose(f(g^-1)) : g in Generators(G) ];

	return hom<G->gl | mats>;

end intrinsic;

//============================================================================
intrinsic SemisimplePart(f::RngUPolElt) -> RngUPolElt
{The semisimple part of f.}

	g := Gcd(f, Derivative(f));
	u := LeadingCoefficient(g);

	return f div g/u;

end intrinsic;

//============================================================================
intrinsic NonZeroPoint(f::RngUPolElt) -> RngElt
{Returns a point x such that f(x) is nonzero.}

	R := BaseRing(Parent(f));

	x := R!0;
	while Evaluate(f,x) eq 0 do
		x +:= 1;
	end while;
	return x;

end intrinsic;

//============================================================================
intrinsic NonZeroPoint(f::FldFunRatMElt) -> RngElt
{Returns a point x such that f(x) is nonzero.}

	R := BaseRing(Parent(f));
	point := [];

	for i:=1 to Ngens(Parent(f)) do
		x := R!0;
		while Evaluate(Denominator(f),i,x) eq 0 or Evaluate(Numerator(f),i,x) eq 0 do
			x +:= 1;
		end while;
		Append(~point, x);
		f := Evaluate(f,i,x);
	end for;

	return point;

end intrinsic;

//============================================================================
intrinsic NonZeroPoint(f::RngMPolElt) -> RngElt
{Returns a point x such that f(x) is nonzero.}

	R := BaseRing(Parent(f));
	point := [];

	for i:=1 to Ngens(Parent(f)) do
		x := R!0;
		while Evaluate(f,i,x) eq 0 do
			x +:= 1;
		end while;
		Append(~point, x);
		f := Evaluate(f,i,x);
	end for;

	return point;

end intrinsic;

//============================================================================
intrinsic NonZeroPoint(F::SeqEnum[FldFunRatMElt]) -> RngElt
{Returns a point x such that f(x) is nonzero for all f in F.}

	R := BaseRing(Parent(F[1]));
	N := Ngens(Parent(F[1]));
	point := [];

	for i:=1 to N do
		x := R!0;
		while exists{f: f in F | Evaluate(Denominator(f),i,x) eq 0 or Evaluate(Numerator(f),i,x) eq 0} do
			x +:= 1;
		end while;
		Append(~point, x);
		F := [ Evaluate(f,i,x) : f in F];
	end for;

	return point;

end intrinsic;


//============================================================================
intrinsic CalogeroMoserCellularCharacters(W::GrpMat, c::Map, fam::SetIndx : vreg:=0, xreg:=0) -> AlgMatElt
{The decomposition matrix of the Calogero-Moser c-cellular characters for W in the family fam of characters of W (fam should be a union of CM families).}

	Representations(~W);
	ReflectionLibrary(~W);

	if #fam eq 1 then
		mults := IdentityMatrix(Integers(), 1);
		return mults;
	end if;

	//If vreg is not provided (vreg=0), then we compute generic Gaudin
	//operators. Otherwise, we compute them specialized them in vreg
	//already.
	print "Computing Gaudin operators";
	if vreg eq 0 then
		gaudins := GaudinOperators(W,c,fam);
	else
		gaudins := GaudinOperatorsSpecialized(W,c,vreg,fam);
	end if;

	// count := 0;
	// for i in fam do
	// 	print i;
	// 	count +:= 1;
	// 	t := Cputime();
	// 	if vreg eq 0 then
	// 		Append(~gaudins, GaudinOperator(W,c,W`Representations[0][i]));
	// 	else
	// 		Append(~gaudins, GaudinOperatorSpecialized(W,c,vreg,W`Representations[0][i]));
	// 	end if;
	// 	print Sprint(Cputime(t))*" seconds";
	// 	//PrintPercentage(count, #fam);
	// end for;

	//Now, we need to find xreg such that the discriminant is non-zero
	//If xreg is provided, we just specialize the Gaudins.
	//Otherwise, the computation is much more complicated.
	//Note: if xreg is provided, the cellular characters cannot be proven
	//to be correct. It's just a heuristic.
	if xreg eq 0 then

		print "Computing characteristic polynomials";
		charpols := [ ];
		count := 0;
		for D in gaudins do
			count +:= 1;
			t := Cputime();
			print "Characteristic polynomial "*Sprint(fam[count]);
			if vreg eq 0 then
				Append(~charpols, CharacteristicPolynomialNaive(D));
			else
				Append(~charpols, CharacteristicPolynomial(D));
			end if;
			print Sprint(Cputime(t))*" seconds";
			PrintPercentage(count, #fam);
		end for;

		print "Computing semisimple parts";
		charpolsss := [];
		count := 0;
		for f in charpols do
			count +:= 1;
			t := Cputime();
			print "Semisimple part "*Sprint(fam[count]);
			Append(~charpolsss, SemisimplePart(f));
			print Sprint(Cputime(t))*" seconds";
			PrintPercentage(count, #fam);
		end for;

		print "Computing semisimple part of product";
		//return charpols,charpolsss; //for debugging

		charpolssssprodss := charpolsss[1];
		for i:=2 to #charpolsss do
			charpolssssprodss *:= charpolsss[i];
			charpolssssprodss := SemisimplePart(charpolssssprodss);
			PrintPercentage(i, #fam-1);
		end for;

		print "Computing discriminant";
		//print charpolssssprodss;
		t := Cputime();
		disc := Discriminant(charpolssssprodss);
		print Sprint(Cputime(t))*" seconds";


		//Commented code
		print "Searching for nonzero point";
		R := BaseRing(gaudins[1]);
		if vreg eq 0 then
			point := NonZeroPoint([disc] cat [&*[ &+[s`Coroot[i]*R.(Dimension(W)+i) : i in [1..Dimension(W)]] : s in W`ReflectionLibraryFlat ]]);
		else
			point := NonZeroPoint(disc);
		end if;

		print point;

	else

		point := Eltseq(xreg);

	end if;

	print "Specializing Gaudin operators";
	gaudinsspec := [* Evaluate(D, point) : D in gaudins *];

	print "Computing characteristic polynomials";
	speccharpols := [];
	count := 0;
	for D in gaudinsspec do
		count +:= 1;
		print "Characteristic polynomial "*Sprint(fam[count]);
		Append(~speccharpols, CharacteristicPolynomial(D));
		PrintPercentage(count, #W`Representations[0]);
	end for;

	//UT: Aug 16, 2023: Here were some unstaged changes.
	//I think this should be correct now.
	//print "Computing semisimple parts";
	//speccharpolsss := [ SemisimplePart(f) : f in speccharpols ];

	print "Factorizing characteristic polynomials";
	//factorizations := [ Factorization(f) : f in speccharpolsss];
	factorizations := [ Factorization(f) : f in speccharpols];

	//print factorizations;

	factors := {};
	for f in factorizations do
		for g in f do
			factors join:={g[1]};
		end for;
	end for;
	factors := SetToSequence(factors);

	print "Computing multiplicities";
	mults := ZeroMatrix(Integers(), #factors, #fam);

	for l:=1 to #factors do
		f := factors[l];
		for i:=1 to #fam do
			fpos := Position([factorizations[i][j][1] : j in [1..#factorizations[i]]], f);
			if fpos eq 0 then
				continue;
			end if;
			mults[l][i] := Dimension(Kernel(Evaluate(f^factorizations[i][fpos][2], gaudinsspec[i])))/Degree(f);
		end for;
	end for;

	return mults;

end intrinsic;

intrinsic CalogeroMoserCellularCharacters(W::GrpMat, c::Map : vreg:=0, xreg:=0) -> List
{The decomposition matrix of the Calogero-Moser c-cellular characters for W, computed per Euler family. The output is a pair, consisting of an Euler family and the decomposition matrix.}

	eulerfam := EulerFamilies(W,c);
	eulerfam := {@ fam[1] : fam in eulerfam @};
	mults := [* <fam,CalogeroMoserCellularCharacters(W,c,fam : vreg:=vreg, xreg:=xreg)> : fam in eulerfam *];
	return mults;

end intrinsic;

//============================================================================
intrinsic CharacteristicPolynomialNaive(M::Mtrx) -> RngUPolElt
{Computes the characteristic polynomial of M in a naive way. Might beform much better for complicated matrices (entries are e.g. large polynomials) of small size.}

	P := PolynomialRing(BaseRing(M));
	t := P.1;
	zero := Zero(P);
	assert Ncols(M) eq Nrows(M);
	n := Ncols(M);
	chi := zero;
	count := 0;
	Sigma := SymmetricGroup(n);
	tI := t*IdentityMatrix(P, n);
	one := One(P);
	N := #Sigma;
	//step := (#Sigma div 12) + 1;
	step := 1;

	//even more naive
	for sigma in Sigma do
		count +:= 1;
		sigmaeltseq := Eltseq(sigma);
		prod := one;
		for i:=1 to n do
			prod *:= tI[i][sigmaeltseq[i]] - M[i][sigmaeltseq[i]];
		end for;
		if Sign(sigma) eq 1 then
			chi +:= prod;
		else
			chi -:= prod;
		end if;
		if count mod step eq 0 then
			PrintPercentage(count, N);
		end if;
	end for;

	return chi;

end intrinsic;

//============================================================================
intrinsic DeterminantNaive(M::Mtrx) -> RngUPolElt
{Computes the characteristic polynomial of M in a naive way. Might beform much better for complicated matrices (entries are e.g. large polynomials) of small size.}

	zero := Zero(BaseRing(M));
	one := One(BaseRing(M));
	assert Ncols(M) eq Nrows(M);
	n := Ncols(M);
	count := 0;
	Sigma := SymmetricGroup(n);
	N := #Sigma;
	//step := (#Sigma div 12) + 1;
	step := 1;

	//even more naive
	det := zero;
	for sigma in Sigma do
		count +:= 1;
		sigmaeltseq := Eltseq(sigma);
		prod := one;
		for i:=1 to n do
			prod *:= M[i][sigmaeltseq[i]];
		end for;
		if Sign(sigma) eq 1 then
			det +:= prod;
		else
			det -:= prod;
		end if;
		if count mod step eq 0 then
			PrintPercentage(count, N);
		end if;
	end for;

	return det;

end intrinsic;

//============================================================================
intrinsic DiscriminantNaive(f::RngUPolElt) -> RngElt
{Computes the discriminant of f in a naive way. Might perform better for small degree polynomials with complicated coefficients.}

	d := Degree(f);
	X := BaseRing(f).1;
	if d eq 4 then
		a := Coefficients(f)[5];
		b := Coefficients(f)[4];
		c := Coefficients(f)[3];
		d := Coefficients(f)[2];
		e := Coefficients(f)[1];
		disc := 256*a^3*e^3 - 192*a^2*b*d*e^2 - 128*a^2*c^2*e^2 + 144*a^2*c*d^2*e - 27*a^2*d^4 + 144*a*b^2*c*e^2 - 6*a*b^2*d^2*e - 80*a*b*c^2*d*e + 18*a*b*c*d^3 + 16*a*c^4*e - 4*a*c^3*d^2 - 27*b^4*e^2 + 18*b^3*c*d*e - 4*b^3*d^3 - 4*b^2*c^3*e + b^2*c^2*d^2;
		return disc;
	end if;

end intrinsic;

//============================================================================
intrinsic ChangeRing(f::RngUPolElt, phi::Map) -> RngUPolElt
{}

	P := PolynomialRing(Codomain(phi));
	return &+[ phi(Coefficients(f)[i])*P.1^(i-1) : i in [1..#Coefficients(f) ]];

end intrinsic;

//============================================================================
intrinsic IsRegular(W::GrpMat, v::ModTupFldElt) -> BoolElt
{Returns true if v is a regular vector, i.e. v is not fixed by any reflection in W.}

	Reflections(~W);
	for s in W`Reflections do
		if v*s eq v then
			return false;
		end if;
	end for;

	return true;

end intrinsic;

intrinsic RandomRegularVector(W::GrpMat) -> ModTupFldElt
{Returns a random regular vector (tries to keep the entries small).}

	N := 10;
	V := VectorSpace(W);
	n := Dimension(V);
	while true do
		i := 1;
		S := {-i..i};
		for j:=1 to N do
			v := V![ Random(S) : k in [1..n] ];
			if IsRegular(W,v) then
				return v;
			end if;
		end for;
		i +:= 1;
	end while;

end intrinsic;


intrinsic GaudinOperatorsForJulia(W::GrpMat, c::Map, fam::SetIndx : file:="") -> MonStgElt
{}

	str := "";
	gaudins := GaudinOperators(W,c,fam);
	R := BaseRing(gaudins[1]); //a rational function field
	Rb := BaseRing(R);
	if Type(Rb) eq FldRat or Type(Rb) eq FldCyc then
		k := Rb;
	elif Type(Rb) eq RngMPol or Type(Rb) eq FldFunRat then
		k := BaseRing(Rb);
	else
		error "Not implemented for this type of base ring";
	end if;

	if Type(k) eq FldRat then
		str *:= "k = QQ";
	elif Type(k) eq FldCyc then
		n := CyclotomicOrder(k);
		str *:= Sprintf("k,z = CyclotomicField(%o, \"z\")", n);
	end if;
	str *:= "\n";

	if Type(Rb) eq FldRat or Type(Rb) eq FldCyc then
		str *:= "P = k";
	elif Type(Rb) eq RngMPol or Type(Rb) eq FldFunRat then
		str *:= "P, ";
		if Rank(Rb) gt 1 then
			str *:= "(";
		end if;
		for i:=1 to Rank(Rb) do
			str *:= Sprint(Rb.i);
			if i lt Rank(Rb) then
				str *:= ",";
			end if;
		end for;
		if Rank(Rb) gt 1 then
			str *:= ")";
		end if;
		str *:= " = PolynomialRing(k, ";
		if Rank(Rb) gt 1 then
			str *:= "[";
		end if;
		for i:=1 to Rank(Rb) do
			str *:= "\""*Sprint(Rb.i)*"\"";
			if i lt Rank(Rb) then
				str *:= ",";
			end if;
		end for;
		if Rank(Rb) gt 1 then
			str *:= "]";
		end if;
		str *:= ")";
	end if;
	str *:= "\n";

	if Type(Rb) eq FldFunRat then
		str *:= "P = FractionField(P)";
		str *:= "\n";
		for i:=1 to Rank(Rb) do
			str *:= Sprint(Rb.i);
			if i lt Rank(Rb) then
				str *:= ",";
			end if;
		end for;
		str *:= " = map(P, gens(base_ring(P)))";
		str *:= "\n";
	end if;

	str *:= "K, (";
	for i:=1 to Rank(R) do
		str *:= Sprint(R.i);
		if i lt Rank(R) then
			str *:= ",";
		end if;
	end for;
	str *:= ") = PolynomialRing(P, [";
	for i:=1 to Rank(R) do
		str *:= "\""*Sprint(R.i)*"\"";
		if i lt Rank(R) then
			str *:= ",";
		end if;
	end for;
	str *:= "])";
	str *:= "\n";

	str *:= "K = FractionField(K)";
	str *:= "\n";
	for i:=1 to Rank(R) do
		str *:= Sprint(R.i);
		if i lt Rank(R) then
			str *:= ",";
		end if;
	end for;
	str *:= " = map(K, gens(base_ring(K)))";
	str *:= "\n";

	if Type(Rb) eq RngMPol or Type(Rb) eq FldFunRat then
		for i:=1 to Rank(Rb) do
			str *:= Sprint(Rb.i);
			if i lt Rank(Rb) then
				str *:= ",";
			end if;
		end for;
		str *:= " = map(K,[";
		for i:=1 to Rank(Rb) do
			str *:= Sprint(Rb.i);
			if i lt Rank(Rb) then
				str *:= ",";
			end if;
		end for;
		str *:= "])";
		str *:= "\n";
	end if;

	for i:=1 to #gaudins do
		N := Ncols(gaudins[i]);
		Dstr := Sprint(Eltseq(gaudins[i]));
		Dstr := Replace(Dstr, "\\\\/", "\\\\/\\\\\/");
		if Type(k) eq FldCyc then
			Dstr := Replace(Dstr, Sprint(k.1), "z");
		end if;
		str *:= Sprintf("D%o = matrix(K, %o, %o, %o)", fam[i], N, N, Dstr);
		str *:= "\n";
	end for;

	str *:= "D = [";
	for i:=1 to #gaudins do
		str *:= Sprintf("D%o", fam[i]);
		if i lt #gaudins then
			str *:= ",";
		end if;
	end for;
	str *:= "]";
	str *:= "\n";

	coroots_str := "coroots = [";
	for j:=1 to #W`ReflectionLibraryFlat do
		s := W`ReflectionLibraryFlat[j];
		coroots_str *:= Sprint(&*[ &+[s`Coroot[i]*R.(Dimension(W)+i) : i in [1..Dimension(W)]]]);
		if j lt #W`ReflectionLibraryFlat then
			coroots_str *:= ",";
		end if;
	end for;
	coroots_str *:= "]";
	coroots_str := Replace(coroots_str, "\\\\/", "\\\\/\\\\\/");
	if Type(k) eq FldCyc then
		coroots_str := Replace(coroots_str, Sprint(k.1), "z");
	end if;
	str *:= coroots_str;

	if file ne "" then
		Write(file, str : Overwrite:=true);
	end if;

	return str;


end intrinsic;
