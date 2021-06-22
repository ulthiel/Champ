/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Intrinsics around group recognition.
*/

forward GroupRecognitionGL;
forward GroupRecognitionSL;
forward GroupRecognitionPGL;
forward GroupRecognitionPSL;
forward GroupRecognitionGU;
forward GroupRecognitionSU;
forward GroupRecognitionPGU;
forward GroupRecognitionPSU;
forward GroupRecognitionSp;
forward GroupRecognitionPSp;
forward GroupRecognitionGOPlus;
forward GroupRecognitionGOMinus;
forward GroupRecognitionGO;
forward GroupRecognitionSOPlus;
forward GroupRecognitionSOMinus;
forward GroupRecognitionSO;
forward GroupRecognitionOmegaPlus;
forward GroupRecognitionOmegaMinus;
forward GroupRecognitionOmega;
forward GroupRecognitionPGOPlus;
forward GroupRecognitionPGOMinus;
forward GroupRecognitionPGO;
forward GroupRecognitionPSOPlus;
forward GroupRecognitionPSOMinus;
forward GroupRecognitionPSO;
forward GroupRecognitionPOmegaPlus;
forward GroupRecognitionPOmegaMinus;
forward GroupRecognitionPOmega;
forward GroupRecognitionAlternatingGroup;
forward GroupRecognitionSymmetricGroup;
forward GroupRecognitionDihedralGroup;
forward GroupRecognitionQuaternionGroup;

declare verbose GrpRec, 5;

//==============================================================================
intrinsic GroupRecognition(G::Grp : Check:=["All"]) -> SeqEnum
/*
    History:
        Friday, June 21, 2013 17:32:28: Initial.
*/
{}

    result := "";

    if Check eq ["All"] then
        Check := ["ClassicalGroups", "ProjectiveClassicalGroups", "AlternatingGroup", "SymmetricGroup", "ShephardTodd"];
    end if;
    if "ClassicalGroups" in Check then
        Check cat:=["GL", "SL", "GU", "SU", "Sp", "GO", "GO+", "GO-", "SO", "SO+", "SO-", "Omega", "Omega+", "Omega-"];
    end if;
    if "ProjectiveClassicalGroups" in Check then
        Check cat:=["PGL", "PSL", "PGU", "PSU", "PSp", "PGO", "PGO+", "PGO-", "PSO", "PSO+", "PSO-", "POmega", "POmega+", "POmega-"];
    end if;


    if "ShephardTodd" in Check then
        result *:= GroupRecognitionShephardTodd(G);
    end if;
    if "GL" in Check then
        vprint GrpRec, 5: "Testing GL.";
        result *:= GroupRecognitionGL(G);
    end if;
    if "SL" in Check then
        vprint GrpRec, 5: "Testing SL.";
        result *:= GroupRecognitionSL(G);
    end if;
    if "PGL" in Check then
        vprint GrpRec, 5: "Testing PGL.";
        result *:= GroupRecognitionPGL(G);
    end if;
    if "GU" in Check then
        vprint GrpRec, 5: "Testing GU.";
        result *:= GroupRecognitionGU(G);
    end if;
    if "SU" in Check then
        vprint GrpRec, 5: "Testing SU.";
        result *:= GroupRecognitionSU(G);
    end if;
    if "PGU" in Check then
        vprint GrpRec, 5: "Testing PGU.";
        result *:= GroupRecognitionPGU(G);
    end if;
    if "PSU" in Check then
        vprint GrpRec, 5: "Testing PSU.";
        result *:= GroupRecognitionPSU(G);
    end if;
    if "Sp" in Check then
        vprint GrpRec, 5: "Testing Sp.";
        result *:= GroupRecognitionSp(G);
    end if;
    if "PSp" in Check then
        vprint GrpRec, 5: "Testing PSp.";
        result *:= GroupRecognitionPSp(G);
    end if;
    if "GO+" in Check then
        vprint GrpRec, 5: "Testing GO+.";
        result *:= GroupRecognitionGOPlus(G);
    end if;
    if "GO-" in Check then
        vprint GrpRec, 5: "Testing GO-.";
        result *:= GroupRecognitionGOMinus(G);
    end if;
    if "GO" in Check then
        vprint GrpRec, 5: "Testing GO.";
        result *:= GroupRecognitionGO(G);
    end if;
    if "SO+" in Check then
        vprint GrpRec, 5: "Testing SO+.";
        result *:= GroupRecognitionSOPlus(G);
    end if;
    if "SO-" in Check then
        vprint GrpRec, 5: "Testing SO-.";
        result *:= GroupRecognitionSOMinus(G);
    end if;
    if "SO" in Check then
        vprint GrpRec, 5: "Testing SO.";
        result *:= GroupRecognitionSO(G);
    end if;
    if "Omega+" in Check then
        vprint GrpRec, 5: "Testing Omega+.";
        result *:= GroupRecognitionOmegaPlus(G);
    end if;
    if "Omega-" in Check then
        vprint GrpRec, 5: "Testing Omega-.";
        result *:= GroupRecognitionOmegaMinus(G);
    end if;
    if "Omega" in Check then
        vprint GrpRec, 5: "Testing Omega.";
        result *:= GroupRecognitionOmega(G);
    end if;
    if "PGO+" in Check then
        vprint GrpRec, 5: "Testing PGO+.";
        result *:= GroupRecognitionPGOPlus(G);
    end if;
    if "PGO-" in Check then
        vprint GrpRec, 5: "Testing PGO-.";
        result *:= GroupRecognitionPGOMinus(G);
    end if;
    if "PGO" in Check then
        vprint GrpRec, 5: "Testing PGO.";
        result *:= GroupRecognitionPGO(G);
    end if;
    if "PSO+" in Check then
        vprint GrpRec, 5: "Testing PSO+.";
        result *:= GroupRecognitionPSOPlus(G);
    end if;
    if "PSO-" in Check then
        vprint GrpRec, 5: "Testing PSO-.";
        result *:= GroupRecognitionPSOMinus(G);
    end if;
    if "PSO" in Check then
        vprint GrpRec, 5: "Testing PSO.";
        result *:= GroupRecognitionPSO(G);
    end if;
    if "PSL" in Check then
        vprint GrpRec, 5: "Testing PSL.";
        result *:= GroupRecognitionPSL(G);
    end if;
    if "POmega+" in Check then
        vprint GrpRec, 5: "Testing POmega+.";
        result *:= GroupRecognitionPOmegaPlus(G);
    end if;
    if "POmega-" in Check then
        vprint GrpRec, 5: "Testing POmega-.";
        result *:= GroupRecognitionPOmegaMinus(G);
    end if;
    if "POmega" in Check then
        vprint GrpRec, 5: "Testing POmega.";
        result *:= GroupRecognitionPOmega(G);
    end if;
    if "AlternatingGroup" in Check then
        vprint GrpRec, 5: "Testing Alternating Group.";
        result *:= GroupRecognitionAlternatingGroup(G);
    end if;
    if "SymmetricGroup" in Check then
        vprint GrpRec, 5: "Testing Symmetric Group.";
        result *:= GroupRecognitionSymmetricGroup(G);
    end if;
    if "DihedralGroup" in Check then
        result *:= GroupRecognitionDihedralGroup(G);
    end if;
    if "QuaternionGroup" in Check then
        result *:= GroupRecognitionQuaternionGroup(G);
    end if;

    return result;

end intrinsic;


//==============================================================================
/*
    Checks if G is isomorphic to GL(n,q).
*/
function GroupRecognitionGL(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=1;
    while 1 gt 0 do
        f := (&*[ (q^n - q^t) : t in [0..n-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        //magma is too stupid: just by order PGL(3,13) and GL(1,810534817) could theoretically be isomorphic and magma computes for ages. But of course PGL is not abelian, so they can't be isomorphic!
        if #{IsAbelian(G), IsAbelian(GL(n,q))} eq 2 then
            continue;
        end if;

        if IsIsomorphic(G,GL(n,q)) then
            result *:= ">> Group isomorphic to general linear group GL("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to SL(n,q) and if G is a subgroup of GL(n,q), checks if it is conjugate to SL(n,q).
*/
function GroupRecognitionSL(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := (&*[ (q^n - q^t) : t in [0..n-1] ]) div (q-1);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,SL(n,q)) then
            result *:= ">> Group isomorphic to special linear group SL("*Sprint(n)*","*Sprint(q)*").\n";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,SL(n,q)) then
                    result *:= ">> Matrix group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to special linear group SL("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PGL(n,q).
*/
function GroupRecognitionPGL(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := (&*[ (q^n - q^t) : t in [0..n-1] ]) div (q-1);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,PGL(n,q)) then
            result *:= ">> Group isomorphic to projective general linear group PGL("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;




//==============================================================================
/*
    Checks if G is isomorphic to PSL(n,q).
*/
function GroupRecognitionPSL(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Rationals());
    n:=2;
    toobig := false;
    while n gt 0 do
        for d in Divisors(n) do
            f := (1/d)*q^( (n*(n-1)) div 2)*ArrayProduct([q^i-1 : i in [2..n]]);
            if Evaluate(f,2) gt order then
                toobig := true;
                break;
            end if;
            roots := Roots(f-order);
            for r in roots do
                if not IsIntegral(r[1]) then
                    continue;
                end if;
                if r[1] lt 2 then
                    continue;
                end if;
                if IsPrimePower(Integers()!r[1]) then
                    candidates join:={<n,Integers()!r[1]>};
                end if;
            end for;
        end for;
        if toobig then
            break;
        end if;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,PSL(n,q)) then
            result *:= ">> Group isomorphic to projective special linear group PSL("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to PSp(n,q).
*/
function GroupRecognitionPSp(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Rationals());
    n:=2;
    toobig := false;
    while n gt 0 do
        m := n div 2;
        for d in Divisors(2) do
            f := (1/d)*q^(m^2)*ArrayProduct([q^(2*i)-1 : i in [1..m]]);
            if Evaluate(f,2) gt order then
                toobig := true;
                break;
            end if;
            roots := Roots(f-order);
            for r in roots do
                if not IsIntegral(r[1]) then
                    continue;
                end if;
                if r[1] lt 2 then
                    continue;
                end if;
                if IsPrimePower(Integers()!r[1]) then
                    candidates join:={<n,Integers()!r[1]>};
                end if;
            end for;
        end for;
        if toobig then
            break;
        end if;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,PSp(n,q)) then
            result *:= ">> Group isomorphic to projective symplectic group PSp("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to GU(n,q).
*/
function GroupRecognitionGU(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := q^( (n^2-n) div 2)*(&*[ (q^i - (-1)^i) : i in [1..n] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:=1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,GU(n,q)) then
            result *:= ">> Group isomorphic to general unitary group GU("*Sprint(n)*","*Sprint(q)*").\n";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q^2),ChangeRing(G,GF(q^2)),GU(n,q)) then
                    result *:= ">> Matrix group conjugate in GL("*Sprint(n)*","*Sprint(q^2)*") to general unitary group GU("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to SU(n,q).
*/
function GroupRecognitionSU(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := q^( (n^2-n) div 2)*ArrayProduct([ (q^i - (-1)^i) : i in [1..n] ]);
        f := f div (q+1);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:=1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,SU(n,q)) then
            result *:= ">> Group isomorphic to special unitary group SU("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q^2),ChangeRing(G,GF(q^2)),SU(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q^2)*") to special unitary group ("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to PGU(n,q).
*/
function GroupRecognitionPGU(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := (q^( (n^2-n) div 2)*(&*[ (q^i - (-1)^i) : i in [1..n] ])) div (q+1);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,PGU(n,q)) then
            result *:= ">> Group isomorphic to projective general unitary group PGU("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PSU(n,q).
*/
function GroupRecognitionPSU(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Rationals());
    n:=2;
    toobig := false;
    while n gt 0 do
        for d in Divisors(n) do
            f := (1/d)*q^( (n^2-n) div 2)*ArrayProduct([ (q^i - (-1)^i) : i in [1..n] ]);
            f := f div (q+1);
            if Evaluate(f,2) gt order then
                toobig := true;
                break;
            end if;
            roots := Roots(f-order);
            for r in roots do
                if not IsIntegral(r[1]) then
                    continue;
                end if;
                if r[1] lt 2 then
                    continue;
                end if;
                if IsPrimePower(Integers()!r[1]) then
                    candidates join:={<n,Integers()!r[1]>};
                end if;
            end for;
        end for;
        if toobig then
            break;
        end if;
        n +:= 1;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,PSU(n,q)) then
            result *:= ">> Group isomorphic to projective special unitary group PSU("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to Sp(n,q).
*/
function GroupRecognitionSp(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Integers());
    n:=2; //n=1 is not really an interesting case, its the trivial group.
    while 1 gt 0 do
        f := &*[ (q^(2*i)-1)*q^(2*i-1) : i in [1..(n div 2)] ];
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(r[1]) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,Sp(n,q)) then
            result *:= ">> Group isomorphic to symplectic group Sp("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,Sp(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to symplectic group Sp("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to GO+(n,q).
*/
function GroupRecognitionGOPlus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        //TODO: magma has some problems checking isomorphy in even characteristic! annoying!
        power,p,k := IsPrimePower(q);
        if IsEven(p) then
            continue;
        end if;

        if IsIsomorphic(G,GeneralOrthogonalGroupPlus(n,q)) then
            result *:= ">> Group isomorphic to general orthogonal group of plus type GO+("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,GeneralOrthogonalGroupPlus(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to general orthogonal group of plus type GO+("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to GO-(n,q).
*/
function GroupRecognitionGOMinus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,GeneralOrthogonalGroupMinus(n,q)) then
            result *:= ">> Group isomorphic to general orthogonal group of minus type GO-("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,GeneralOrthogonalGroupMinus(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to general orthogonal group of minus type GO-("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to GO(n,q).
*/
function GroupRecognitionGO(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    //for odd characteristic
    n:=3;
    while true do
        m := (n-1) div 2;
        f := 2*q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=3;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsEven(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,GeneralOrthogonalGroup(n,q)) then
            result *:= ">> Group isomorphic to general orthogonal group GO("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,GeneralOrthogonalGroup(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to general orthogonal group GO("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to SO+(n,q).
*/
function GroupRecognitionSOPlus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic (I guessed the order, but it should be true)
    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsEven(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        //MAGMA IS SO STUPID. THERE IS AN INTERNAL ERROR WITH THE ORDER OF SOPLUS HERE. WE FIX IT.
        T :=  SpecialOrthogonalGroupPlus(n,q);
        Torder := 1;
        power,p,k := IsPrimePower(q);
        if IsOdd(p) then
            m := n div 2;
            Torder := q^(m*(m-1))*(q^m-1);
            Torder *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        else
            Torder := 2*q^(m*(m-1))*(q^m-1);
            Torder *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        end if;
        AssertAttribute(T,"Order",Torder);

        if IsIsomorphic(G,T) then
            result *:= ">> Group isomorphic to special orthogonal group of plus type SO+("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,T) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to special orthogonal group of plus type SO+("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to SO-(n,q).
*/
function GroupRecognitionSOMinus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic (I guessed the order, but it should be true)
    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsEven(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,SpecialOrthogonalGroupMinus(n,q)) then
            result *:= ">> Group isomorphic to special orthogonal group of minus type SO-("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,SpecialOrthogonalGroupMinus(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to special orthogonal group of minus type SO-("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to SO(n,q).
*/
function GroupRecognitionSO(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,SpecialOrthogonalGroup(n,q)) then
            result *:= ">> Group isomorphic to special orthogonal group SO("*Sprint(n)*","*Sprint(q)*").\n";

            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,SpecialOrthogonalGroup(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to special orthogonal group SO("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to Omega+(n,q).
*/
function GroupRecognitionOmegaPlus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //order of Omega+(4,2) is 36; not covered by formula above.
    if order eq 36 then
        candidates join:={<4,2>};
    end if;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,OmegaPlus(n,q)) then
            result *:= ">> Group isomorphic to Omega+("*Sprint(n)*","*Sprint(q)*").\n";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,OmegaPlus(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to Omega+("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to Omega-(n,q).
*/
function GroupRecognitionOmegaMinus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,OmegaMinus(n,q)) then
            result *:= ">> Group isomorphic to Omega-("*Sprint(n)*","*Sprint(q)*").\n";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,OmegaMinus(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to Omega-("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to Omega(n,q).
*/
function GroupRecognitionOmega(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //for odd characteristic
    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,Omega(n,q)) then
            result *:= ">> Group isomorphic to Omega("*Sprint(n)*","*Sprint(q)*").\n";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)),Fld) and Characteristic(BaseRing(G)) ne 0 and #BaseRing(G) eq q and Dimension(G) eq n then
                if IsConjugate(GL(n,q),G,Omega(n,q)) then
                    result *:= ">> Group conjugate in GL("*Sprint(n)*","*Sprint(q)*") to Omega("*Sprint(n)*","*Sprint(q)*").\n";
                end if;
            end if;
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to PGO+(n,q).
*/
function GroupRecognitionPGOPlus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsEven(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PGOPlus(n,q)) then
            result *:= ">> Group isomorphic to projective general orthogonal group of plus type PGO+("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PGO-(n,q).
*/
function GroupRecognitionPGOMinus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := 2*q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power and IsEven(p) then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //some fixes for small orders not covered by formula:
    if order eq 2 then
        candidates join:={<2,3>};
    end if;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PGOMinus(n,q)) then
            result *:= ">> Group isomorphic to projective general orthogonal group of minus type PGO-("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PGO(n,q).
*/
function GroupRecognitionPGO(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PGO(n,q)) then
            result *:= ">> Group isomorphic to projective general orthogonal group PGO("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PSO+(n,q).
*/
function GroupRecognitionPSOPlus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f *:= 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PSOPlus(n,q)) then
            result *:= ">> Group isomorphic to projective special orthogonal group of plus type PSO+("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to PSO-(n,q).
*/
function GroupRecognitionPSOMinus(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //for odd characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //for even characteristic
    n:=2;
    while true do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f *:= 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PSOMinus(n,q)) then
            result *:= ">> Group isomorphic to projective special orthogonal group of minus type PSO-("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to PSO(n,q).
*/
function GroupRecognitionPSO(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Integers());

    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(r[1]);
            if power then
                candidates join:={<n,r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,PSO(n,q)) then
            result *:= ">> Group isomorphic to projective special orthogonal group PSO("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;


//==============================================================================
/*
    Checks if G is isomorphic to POmegaPlus(n,q).
*/
function GroupRecognitionPOmegaPlus(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Rationals());

    //odd characteristic
    n:=2;
    toobig := false;
    while n gt 0 do
        m := n div 2;
        for d in Divisors(4) do
            f := (1/d)*q^(m*(m-1))*(q^m-1);
            f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
            if Evaluate(f,2) gt order then
                toobig := true;
                break;
            end if;
            roots := Roots(f-order);
            for r in roots do
                if not IsIntegral(r[1]) then
                    continue;
                end if;
                if r[1] lt 2 then
                    continue;
                end if;
                if IsPrimePower(Integers()!r[1]) and IsOdd(Integers()!r[1]) then
                    candidates join:={<n,Integers()!r[1]>};
                end if;
            end for;
        end for;
        if toobig then
            break;
        end if;
        n +:= 2;
    end while;

    //even characteristic
    n:=2;
    toobig := false;
    while n gt 0 do
        m := n div 2;
        f := q^(m*(m-1))*(q^m-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
           toobig := true;
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(Integers()!r[1]) and IsEven(Integers()!r[1]) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        if toobig then
            break;
        end if;
        n +:= 2;
    end while;

    //one fix
    if order eq 1 then
        candidates join:={<2,2>};
    end if;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,POmegaPlus(n,q)) then
            result *:= ">> Group isomorphic to POmega+("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to POmegaMinus(n,q).
*/
function GroupRecognitionPOmegaMinus(G)

    result := "";
    order := #G;

    candidates := {};
    P<q> := PolynomialRing(Rationals());

    //odd characteristic
    n:=2;
    toobig := false;
    while n gt 0 do
        m := n div 2;
        for d in Divisors(4) do
            f := (1/d)*q^(m*(m-1))*(q^m+1);
            f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
            if Evaluate(f,2) gt order then
                toobig := true;
                break;
            end if;
            roots := Roots(f-order);
            for r in roots do
                if not IsIntegral(r[1]) then
                    continue;
                end if;
                if r[1] lt 2 then
                    continue;
                end if;
                if IsPrimePower(Integers()!r[1]) and IsOdd(Integers()!r[1]) then
                    candidates join:={<n,Integers()!r[1]>};
                end if;
            end for;
        end for;
        if toobig then
            break;
        end if;
        n +:= 2;
    end while;

    //even characteristic
    n:=2;
    toobig := false;
    while n gt 0 do
        m := n div 2;
        f := q^(m*(m-1))*(q^m+1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
           toobig := true;
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            if IsPrimePower(Integers()!r[1]) and IsEven(Integers()!r[1]) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        if toobig then
            break;
        end if;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];
        if IsIsomorphic(G,POmegaMinus(n,q)) then
            result *:= ">> Group isomorphic to POmega-("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to POmega(n,q).
*/
function GroupRecognitionPOmega(G)

    result := "";
    order := #G;

    candidates := {};

    P<q> := PolynomialRing(Rationals());

    //odd characteristic
    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        f := f div 2;
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsOdd(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    //even characteristic
    n:=1;
    while true do
        m := (n-1) div 2;
        f := q^(m^2)*(q^(2*m)-1);
        f *:= ArrayProduct([ q^(2*i)-1 : i in [1..m-1] ]);
        if Evaluate(f,2) gt order then
            break;
        end if;
        roots := Roots(f-order);
        for r in roots do
            if not IsIntegral(r[1]) then
                continue;
            end if;
            if r[1] lt 2 then
                continue;
            end if;
            power,p,k := IsPrimePower(Integers()!r[1]);
            if power and IsEven(p) then
                candidates join:={<n,Integers()!r[1]>};
            end if;
        end for;
        n +:= 2;
    end while;

    for cand in candidates do
        n := cand[1];
        q := cand[2];

        if IsIsomorphic(G,POmega(n,q)) then
            result *:= ">> Group isomorphic to POmega("*Sprint(n)*","*Sprint(q)*").\n";
        end if;
    end for;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to an alternating group.
*/
function GroupRecognitionAlternatingGroup(G)

    result := "";
    order := #G;

    n:=2;
    while 1 gt 0 do
        o := Factorial(n) div 2;
        if o ge order then
            if IsIsomorphic(G,AlternatingGroup(n)) then
                result *:= ">> Group isomorphic to alternating group of degree "*Sprint(n)*".\n";
            end if;
            break;
        else
            n+:=1;
        end if;
    end while;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to a symmetric group.
*/
function GroupRecognitionSymmetricGroup(G)

    result := "";
    order := #G;

    n:=2;
    while 1 gt 0 do
        o := Factorial(n);
        if o ge order then
            if IsIsomorphic(G,SymmetricGroup(n)) then
                result *:= ">> Group isomorphic to symmetric group of degree "*Sprint(n)*".\n";
            end if;
            break;
        else
            n+:=1;
        end if;
    end while;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to a dihedral group.
*/
function GroupRecognitionDihedralGroup(G)

    result := "";
    order := #G;

    if order mod 2 ne 0 then
        return "";
    end if;

    n := order div 2;

    if n in {1,2} then
        if #Center(G) ne 2*n then
            return "";
        end if;
    else
        if IsOdd(n) then
            if #Center(G) ne 1 then
                return "";
            end if;
        else
            if #Center(G) ne 2 then
                return "";
            end if;
        end if;
    end if;

    if IsIsomorphic(G,DihedralGroup(n)) then
        result cat:= ">> Group isomorphic to dihedral group of degree "*Sprint(n)*".\n";
    end if;

    return result;

end function;

//==============================================================================
/*
    Checks if G is isomorphic to a generalized quaternion group.
*/
function GroupRecognitionQuaternionGroup(G)

    result := "";
    order := #G;

    if order mod 4 ne 0 then
        return [];
    end if;

    n := order div 4;
    Q := Group<x,y | x^(2*n), y^4, x^n = y^2, y^-1*x*y=x^-1>;
    //the order of the center of the generalized quaternion group is equal to 2.
    if #Center(G) ne 2 then
        return "";
    end if;

    if IsIsomorphic(G,Q) then
        result *:= ">> Group isomorphic to generalized quaternion group of degree "*Sprint(n)*".\n";
    end if;

    return result;

end function;

//==============================================================================
intrinsic GroupRecognitionShephardTodd(G::Grp) -> List
/*
    History:
        Monday, July 01, 2013 09:41:19: Initial.
*/
{}
    order := #G;

    result := "";

    if order eq 24 then
        iso, f:= IsIsomorphic(G, ShephardTodd(4));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G4.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G4.\n";
            end if;
        end if;
    elif order eq 48 then
        iso, f:= IsIsomorphic(G, ShephardTodd(6));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G6.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G6.";
            end if;
        end if;
        iso, f:= IsIsomorphic(G, ShephardTodd(12));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G12.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G12.";
            end if;
        end if;
    elif order eq 72 then
        iso, f:= IsIsomorphic(G, ShephardTodd(5));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G5.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G5.";
            end if;
        end if;
    elif order eq 96 then
         iso, f:= IsIsomorphic(G, ShephardTodd(8));
         if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G8.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G8.";
            end if;
        end if;
        iso, f:= IsIsomorphic(G, ShephardTodd(13));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G13.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G13.";
            end if;
        end if;
    elif order eq 120 then
         iso, f:= IsIsomorphic(G, ShephardTodd(23));
         if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G23.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G23.";
            end if;
        end if;
    elif order eq 144 then
        iso, f:= IsIsomorphic(G, ShephardTodd(7));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G7.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G7.";
            end if;
        end if;
        iso, f:= IsIsomorphic(G, ShephardTodd(14));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G14.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G14.";
            end if;
        end if;
    elif order eq 192 then
        iso, f:= IsIsomorphic(G, ShephardTodd(9));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G9.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G9.";
            end if;
        end if;
    elif order eq 240 then
        iso, f:= IsIsomorphic(G, ShephardTodd(22));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G22.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G22.";
            end if;
        end if;
    elif order eq 288 then
        iso, f:= IsIsomorphic(G, ShephardTodd(10));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G10.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G10.";
            end if;
        end if;
        iso, f:= IsIsomorphic(G, ShephardTodd(15));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G15.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G15.";
            end if;
        end if;
    elif order eq 336 then
        iso, f:= IsIsomorphic(G, ShephardTodd(24));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G24.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G24.";
            end if;
        end if;
    elif order eq 360 then
        iso, f:= IsIsomorphic(G, ShephardTodd(20));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G20.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G20.";
            end if;
        end if;
    elif order eq 576 then
        iso, f:= IsIsomorphic(G, ShephardTodd(11));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G11.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G11.";
            end if;
        end if;
    elif order eq 600 then
        iso, f:= IsIsomorphic(G, ShephardTodd(16));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G16.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G16.";
            end if;
        end if;
    elif order eq 648 then
        iso, f:= IsIsomorphic(G, ShephardTodd(25));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G25.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G25.";
            end if;
        end if;
    elif order eq 720 then
        iso, f:= IsIsomorphic(G, ShephardTodd(21));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G21.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G21.";
            end if;
        end if;
    elif order eq 1152 then
        iso, f:= IsIsomorphic(G, ShephardTodd(28));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G28.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G28.";
            end if;
        end if;
    elif order eq 1200 then
        iso, f:= IsIsomorphic(G, ShephardTodd(17));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G17.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G17.";
            end if;
        end if;
    elif order eq 1296 then
        iso, f:= IsIsomorphic(G, ShephardTodd(26));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G26.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G26.";
            end if;
        end if;
    elif order eq 1800 then
        iso, f:= IsIsomorphic(G, ShephardTodd(18));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G18.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G18.";
            end if;
        end if;
    elif order eq 2160 then
        iso, f:= IsIsomorphic(G, ShephardTodd(27));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G27.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G27.";
            end if;
        end if;
    elif order eq 3600 then
        iso, f:= IsIsomorphic(G, ShephardTodd(19));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G19.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G19.";
            end if;
        end if;
    elif order eq 7680 then
        iso, f:= IsIsomorphic(G, ShephardTodd(29));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G29.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G29.";
            end if;
        end if;
    elif order eq 14400 then
        iso, f:= IsIsomorphic(G, ShephardTodd(30));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G30.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G30.";
            end if;
        end if;
    elif order eq 46080 then
        iso, f:= IsIsomorphic(G, ShephardTodd(31));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G31.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G31.";
            end if;
        end if;
    elif order eq 51840 then
        iso, f:= IsIsomorphic(G, ShephardTodd(33));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G33.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G33.";
            end if;
        end if;
        iso, f:= IsIsomorphic(G, ShephardTodd(35));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G35.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G35.";
            end if;
        end if;
    elif order eq 155520 then
        iso, f:= IsIsomorphic(G, ShephardTodd(32));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G32.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G32.";
            end if;
        end if;
    elif order eq 2903040 then
        iso, f:= IsIsomorphic(G, ShephardTodd(36));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G36.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G36.";
            end if;
        end if;
    elif order eq 39191040 then
        iso, f:= IsIsomorphic(G, ShephardTodd(34));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G34.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G34.";
            end if;
        end if;
    elif order eq 696729600 then
        iso, f:= IsIsomorphic(G, ShephardTodd(37));
        if iso then
            result *:= "Group isomorphic to exceptional complex reflection group of Shephard-Todd type G37.";
            if Type(G) eq GrpMat and ISA(Type(BaseRing(G)), Fld) and Characteristic(BaseRing(G)) eq 0 then
                result *:= "Matrix group equivalent to exceptional complex reflection group of Shephard-Todd type G37.";
            end if;
        end if;
    end if;

    return result;

end intrinsic;
