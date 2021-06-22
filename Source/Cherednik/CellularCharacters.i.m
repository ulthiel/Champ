/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Joint work with Cédric Bonnafé (Montpellier).
*/



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
intrinsic GaudinOperator(G::GrpMat, c::Map, y::ModTupFldElt) -> AlgMatElt
{The action of the Gaudin operator at y on the full group algebra.}

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
intrinsic GaudinOperator(G::GrpMat, c::Map, y::ModTupFldElt, rho::Map) -> AlgMatElt
{The action of the Gaudin operator at y on rho.}


    NumberingMap(~G);
    ReflectionLibrary(~G);
    K := RationalFunctionField(Codomain(c), Dimension(G)); //k(V), k=Codomain(c)
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
intrinsic GaudinOperatorSpecialized(G::GrpMat, c::Map, y::ModTupFldElt, vreg::ModTupFldElt, rho::Map) -> AlgMatElt
{The action of the Gaudin operator at y specialized in the regular vector rho on rho.}

    NumberingMap(~G);
    ReflectionLibrary(~G);
    K := Codomain(c);
    D:=ZeroMatrix(K, Degree(Codomain(rho)), Degree(Codomain(rho)));
    for s in G`ReflectionLibraryFlat do
     	//coroot := &+[s`Coroot[i]*K.i : i in [1..Ngens(K)]];
        D +:=  s`Eigenvalue/(s`Eigenvalue-1)*c(s`ReflectionClass)*CanonicalPairing(y,s`Coroot)/CanonicalPairing(vreg,s`Coroot)*ChangeRing(rho(s`Element), K);
    end for;

    return D;

end intrinsic;

/*
//============================================================================
intrinsic GaudinOperatorModified(G::GrpMat, c::Map, y::ModTupFldElt, rho::Map) -> AlgMatElt
{The action of the Gaudin operator at y on rho.}

    NumberingMap(~G);
    ReflectionLibrary(~G);
    K := PolynomialRing(Codomain(c), #G`ReflectionLibraryFlat); //k(V), k=Codomain(c)
    AssignNames(~K, ["C"*Sprint(i) : i in [1..#G`ReflectionLibraryFlat]]);
    D:=ZeroMatrix(K, Degree(Codomain(rho)), Degree(Codomain(rho)));

    for i:=1 to #G`ReflectionLibraryFlat do
     	s := G`ReflectionLibraryFlat[i];
     	//coroot := &+[s`Coroot[i]*K.i : i in [1..Ngens(K)]];
        D +:=  s`Eigenvalue*c(s`ReflectionClass)*CanonicalPairing(y,s`Coroot)*K.i*ChangeRing(rho(s`Element), K);
    end for;

    return D;

end intrinsic;
*/

//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map) -> AlgMatElt
{The full Gaudin operator.}
	K := RationalFunctionField(Codomain(c), Dimension(G));
	S := RationalFunctionField(Codomain(c), 2*Dimension(G));
	phi := hom<K->S | [S.(i+Dimension(G)) : i in [1..Dimension(G)]]>;
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]] cat ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(ChangeRing(GaudinOperator(G,c,V.i),K),phi) : i in [1..Dimension(G)]];

end intrinsic;

//============================================================================
intrinsic GaudinOperator(G::GrpMat, c::Map, rho::Map) -> AlgMatElt
{The Gaudin operator acting on rho.}
	K := RationalFunctionField(Codomain(c), Dimension(G));
	S := RationalFunctionField(Codomain(c), 2*Dimension(G));
	phi := hom<K->S | [S.(i+Dimension(G)) : i in [1..Dimension(G)]]>;
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]] cat ["y"*Sprint(i) : i in [1..Dimension(G)] ]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(ChangeRing(GaudinOperator(G,c,V.i,rho),K),phi) : i in [1..Dimension(G)]];

end intrinsic;

//============================================================================
intrinsic GaudinOperatorSpecialized(G::GrpMat, c::Map, vreg::ModTupFldElt, rho::HomGrp) -> AlgMatElt
{The Gaudin operator specialized in the regular vector vreg acting on rho.}
	K := Codomain(c);
	S := PolynomialRing(Codomain(c), Dimension(G) : Global:=true);
	AssignNames(~S, ["X"*Sprint(i) : i in [1..Dimension(G)]]);
	V := VectorSpace(G);

	return &+[ S.i*ChangeRing(GaudinOperatorSpecialized(G,c,V.i,vreg,rho),K) : i in [1..Dimension(G)]];

end intrinsic;

//============================================================================
intrinsic CalogeroMoserCellularCharacters(W::GrpMat, c::Map, fam::Setq : vreg:=0, xreg:=0) -> AlgMatElt
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
	gaudins := [* *];
	count := 0;
	for i in fam do
		print i;
		count +:= 1;
		t := Cputime();
		if vreg eq 0 then
			Append(~gaudins, GaudinOperator(W,c,W`Representations[0][i]));
		else
			Append(~gaudins, GaudinOperatorSpecialized(W,c,vreg,W`Representations[0][i]));
		end if;
		print Sprint(Cputime(t))*" seconds";
		PrintPercentage(count, #fam);
	end for;

	//No, we need to find yreg such that the discriminant is non-zero
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
			print "Characteristic polynomial "*Sprint(count);
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
			print "Semisimple part "*Sprint(count);
			Append(~charpolsss, SemisimplePart(f));
			print Sprint(Cputime(t))*" seconds";
			PrintPercentage(count, #fam);
		end for;

		//charpolsssprodss := SemisimplePart(&*charpolsss);
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
		print "Characteristic polynomial "*Sprint(count);
		Append(~speccharpols, CharacteristicPolynomialNaive(D));
		PrintPercentage(count, #W`Representations[0]);
	end for;

	print "Computing semisimple parts";
	speccharpolsss := [ SemisimplePart(f) : f in speccharpols ];

	print "Factorizing semisimple parts";
	factorizations := [ Factorization(f) : f in speccharpolsss];

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
