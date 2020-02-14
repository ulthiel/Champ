/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Some tools for Euler families
*/


//=============================================================================
intrinsic FakeDegreeResidues(~G::GrpMat)
/*
    History:
        Sunday, October 20, 2013 16:19:51: Ported from old project.
*/
{}
	if assigned G`FakeDegreeResidues then
        return;
    end if;

    FakeDegrees(~G);
    CoinvariantAlgebraPoincareSeries(~G);

    G`FakeDegreeResidues := [];
    P<t> := PolynomialRing(Integers():Global:=false);

    p := Characteristic(BaseRing(G));
    Modules(~G,p);

    for i:=1 to #G`Modules[p] do
        num := (Integers()!Dimension(G`Modules[p][i]))*t^Valuation(G`FakeDegrees[i])*G`CoinvariantAlgebraPoincareSeries;
        res := num mod G`FakeDegrees[i];
        Append(~G`FakeDegreeResidues, res);
    end for;

end intrinsic;

//==============================================================================
intrinsic SupersingularRepresentations(~G::GrpMat)
/*
    History:
        Sunday, October 20, 2013 16:33:12: Initial.
*/
{}

    if assigned G`SupersingularRepresentations then
        return;
    end if;

    FakeDegreeResidues(~G);
    p := Characteristic(BaseRing(G));

    G`SupersingularRepresentations := [i : i in [1..#G`Representations[p]] | G`FakeDegreeResidues[i] ne 0 ];

end intrinsic;


//==============================================================================
intrinsic IsSupersingular(f::RngUPolElt[RngInt], g::RngUPolElt[RngInt]) -> SetEnum
/*
    History:
        Saturday, August 24, 2013 21:20:24: Initial.
*/
{True iff (f,g) is supersingular.}

    fact := Factorization(f);
    mults := [ fact[i][2] : i in [1..#fact] ];
    dom := CartesianProduct([ [0..mults[i]] : i in [1..#fact] ]);

    for d in dom do
        q := &*[ fact[i][1]^d[i] : i in [1..#fact] ];
        print q;
    end for;

    print fact;
    print dom;
    return {};

end intrinsic;


//==============================================================================
intrinsic EulerScalar(G::GrpMat, c::Map, rho::Map : Reversed:=false) -> RngElt
/*
    Intrinsic: EulerScalar

    History:
        Monday, September 23, 2013 22:16:33: Initial.
*/
{}
    //first have to find common ring for computations
    R := Codomain(c);
    K := CommonOverfield( BaseRing(G), BaseRing(Codomain(rho)) );
    if Type(R) eq RngMPol or Type(R) eq FldFunRat then
        if BaseRing(R) ne K then
            R := ChangeRing(R, K);
        end if;
    end if;

    val := Zero(R);
    rhodim := Ncols(rho(G.1));

    for s in G`ReflectionLibraryFlat do
        eps := s`Eigenvalue;
        term := R!(K!(1/(eps-1))*K!rho(s`Element)[1][1])*c(s`ReflectionClass);
        if not Reversed then
            term *:= eps;
        end if;
        val +:= term;
    end for;

    return val;

end intrinsic;


//==============================================================================
intrinsic EulerScalar(G::GrpMat, c::Map, chi::AlgChtrElt : Reversed:=false) -> RngElt
/*
    History:
        Monday, September 23, 2013 22:16:33: Initial.
*/
{}

    //first have to find common ring for computations
    R := Codomain(c);
    K := CommonOverfield( BaseRing(G), BaseRing(chi) );
    if Type(R) eq RngMPol or Type(R) eq FldFunRat then
        if BaseRing(R) ne K then
            R := ChangeRing(R, K);
        end if;
    end if;



    val := Zero(R);
    rhodim := chi(1);
    for s in G`ReflectionLibraryFlat do
        eps := s`Eigenvalue;
        term := R!(K!(1/(eps-1))*K!chi(s`Element))*R!c(s`ReflectionClass);
        if not Reversed then
            term *:= eps;
        end if;
        val +:= term;
    end for;

    return R!(K!(1/rhodim))*val;

end intrinsic;

//==============================================================================
intrinsic EulerFamilies(G::GrpMat, c::Map : Print:=false, Reversed:=false) -> SetEnum
/*
    Intrinsic: EulerFamilies

    History:
        Monday, September 23, 2013 22:24:05: Initial.
*/
{}

    p := Characteristic(BaseRing(G));
    if p eq 0 then
        CharacterTable(~G);

        //first have to find common ring for computations
        R := Codomain(c);
        K := CommonOverfield(BaseRing(G), CommonOverfield([* BaseRing(chi) : chi in G`CharacterTable *] ));
        if Type(R) eq RngMPol or Type(R) eq FldFunRat then
            if BaseRing(R) ne K then
                R := ChangeRing(R, K);
            end if;
        end if;

        eulerscalarsarray := [ R!EulerScalar(G, c, G`CharacterTable[i] : Reversed:=Reversed) : i in [1..#G`CharacterTable] ];
    else
        Representations(~G, p);
        eulerscalarsarray := [ EulerScalar(G, c, G`Representations[p][i] : Reversed:=Reversed) : i in [1..#G`Representations[p]] ];
    end if;
    eulerscalarsset := SequenceToSet(eulerscalarsarray);
    eulerfamilies := {@@};
    chars := {@ i : i in [1..#eulerscalarsarray] @};
    for x in eulerscalarsset do
        eulerscalarsset diff:={x};
        fam := {@ i : i in chars | eulerscalarsarray[i] eq x @};
        Diff(~chars, fam);
        eulerfamilies join:={@ < fam, x> @};
    end for;

    if Print then
        str := "";
        for j:=1 to #eulerfamilies do
            fam := eulerfamilies[j][1];
            str *:= "\\lbrace";
            for i:=1 to #fam do
                if p ne 0 then
                    str *:= Sprint(G`RepresentationNames[p][fam[i]]);
                else
                    str *:= Sprint(G`CharacterNames[fam[i]]);
                end if;
                if i lt #fam then
                    str *:= ", ";
                end if;
            end for;
            str *:= "\\rbrace";
            if j lt #eulerfamilies then
                str *:= ", ";
            end if;
        end for;
        print str;
    end if;

    return eulerfamilies;

end intrinsic;

//==============================================================================
intrinsic EulerFamilies(G::GrpMat) -> SetEnum
{}

    c := CherednikParameter(G);
    return EulerFamilies(G,c);

end intrinsic;

//==============================================================================
intrinsic GoodAndBadEulerFamilies(G::GrpMat, c::Map : Reversed:=false) -> SetEnum, SetEnum
/*
    History:
        Sunday, October 20, 2013 16:43:34: Initial.
*/
{}
    p := Characteristic(BaseRing(G));
    Modules(~G,p);

    EU := EulerFamilies(G,c : Reversed:=Reversed);
    eulerfams := {@ x[1] : x in EU @};
    SupersingularRepresentations(~G);

    supersings := SequenceToSet(G`SupersingularRepresentations);

    good := {@ x : x in eulerfams | #x eq 1 or (#x eq 2 and #(x meet supersings) ge 1) or (#x eq 3 and #(x meet supersings) eq 3) @};
    bad := eulerfams diff good;
    return good,bad;

end intrinsic;

//==============================================================================
intrinsic PrintGoodAndBadEulerFamilies(G::GrpMat,c::Map : Reversed:=false)
/*
    History:
        Sunday, October 20, 2013 16:49:47: Initial.
*/
{}

    p := Characteristic(BaseRing(G));
    good,bad := GoodAndBadEulerFamilies(G,c : Reversed:=Reversed);

    print "Good Euler c-families:";
    str := "";
    for j:=1 to #good do
        x := good[j];
        if assigned G`RepresentationNames and IsDefined(G`RepresentationNames, p) then
            str *:= "\\lbrace ";
            for i:=1 to #x do
                str *:= G`RepresentationNames[p][x[i]];
                if i lt #x then
                    str *:= ", ";
                end if;
            end for;
            str *:= " \\rbrace";
        else
            str *:= Sprint(x);
        end if;
        if j lt #good then
            str *:= ", ";
        end if;
    end for;
    print str;

    print "Bad Euler c-families:";
    str := "";
    for j:=1 to #bad do
        x := bad[j];
        if assigned G`RepresentationNames and IsDefined(G`RepresentationNames, p) then
            str *:= "\\lbrace ";
            for i:=1 to #x do
                str *:= G`RepresentationNames[p][x[i]];
                if i lt #x then
                    str *:= ", ";
                end if;
            end for;
            str *:= " \\rbrace";
        else
            str *:= Sprint(x);
        end if;
        if j lt #bad then
            str *:= ", ";
        end if;
    end for;
    print str;

end intrinsic;

//==============================================================================
intrinsic EulerScheme(G::GrpMat, c::Map : Reversed:=false) -> SeqEnum
/*
    Intrinsic: EulerScheme

    History:
        * Wednesday, April 2, 2014 at 13:59:44: Optimized; not using RadicalDecomposition any more.

        * Tuesday, September 24, 2013 19:45:39: Initial.
*/
{}

    p := Characteristic(BaseRing(G));
    if p eq 0 then
        CharacterTable(~G);
        //first have to find common ring for computations
        R := Codomain(c);
        K := CommonOverfield(BaseRing(G), CommonOverfield([* BaseRing(chi) : chi in G`CharacterTable *] ));
        if Type(R) eq RngMPol or Type(R) eq FldFunRat then
            if BaseRing(R) ne K then
                R := ChangeRing(R, K);
            end if;
        end if;
        eulerscalarsset := { R!EulerScalar(G, c, G`CharacterTable[i] : Reversed:=Reversed) : i in [1..#G`CharacterTable] };
    else
        Representations(~G, p);
        eulerscalarsset := { EulerScalar(G, c, G`Representations[p][i] : Reversed:=Reversed) : i in [1..#G`Representations[p]] };
    end if;

    planes := { x-y : x,y in eulerscalarsset | x ne y };

    P := Universe(planes);

    //normalize
    normplanes := {P | };
    for x in planes do
        Coeffs := Coefficients(x);
        Mons := Monomials(x);
        ZCoeffs := [Integers()!Coeffs[i] : i in [1..#Coeffs]];
        LCF := GreatestCommonDivisor(ZCoeffs)*Sign(Integers()!LeadingCoefficient(x));
        for i:=1 to #Coeffs do
            Coeffs[i] := Coeffs[i]/LCF;
        end for;
        normplanes join:={Polynomial(Coeffs,Mons)};
    end for;

    return Sort(SetToSequence(normplanes));

end intrinsic;

//==============================================================================
intrinsic EulerScheme(G::GrpMat : Reversed:=false) -> RngMPol
/*
    History:
        Tuesday, September 24, 2013 19:45:33: Initial.
*/
{}

    c := CherednikParameter(G:Rational:=false, Type:="GGOR");
    CherednikParameterSpace(~G);
    return ChangeUniverse(EulerScheme(G,c : Reversed:=Reversed), G`CherednikParameterSpace);

end intrinsic;

//==============================================================================
intrinsic SymmetricReflectionGroupGenericEulerFamilies(n::RngIntElt) -> SetIndx
/*
    History:
        Tuesday, September 24, 2013 22:00:24: Initial.
*/
{The generic Euler familie for a non-modular symmetric reflection group.}
    parts := Partitions(n);
    eulerscalarsarray := [ &+[ lambda[j]*(lambda[j] - 2*j + 1) : j in [1..#lambda]] : lambda in parts ];
    eulerscalarsset := SequenceToSet(eulerscalarsarray);
    eulerfamilies := {@@};
    chars := {@ i : i in [1..#parts] @};
    for x in eulerscalarsset do
        eulerscalarsset diff:={x};
        fam := {@ i : i in chars | eulerscalarsarray[i] eq x @};
        Diff(~chars, fam);
        eulerfamilies join:={@ fam @};
    end for;

    return eulerfamilies;

end intrinsic;
