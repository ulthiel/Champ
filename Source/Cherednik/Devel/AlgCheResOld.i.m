/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Intrinsics for rational Cherednik algebras
*/

forward ApplyRuleToTerm;
forward RewriteMonomial;
forward Normalize;
forward NormalizeTerm;
forward FinalTermNormalization;
forward FinalNormalization;
forward VectorSpaceMap;

//==============================================================================
intrinsic RestrictedRationalCherednikAlgebra(G::GrpMat, c::Map) -> AlgCheRes
/*
    History:
        Sunday, October 27, 2013 00:03:54: Initial.
*/
{}
    H := New(AlgCheRes);
    H`ProductsDomain := [];
    H`ProductsCodomain := [**];
    H`Generators := [**];
    H`KeepProducts := true;
    H`RewriteRulesDomain := [];
    H`RewriteRulesCodomain := [**];
    H`RewriteRulesCounter := [];

    //ALGEBRA SPECIFIC FROM HERE
    ReflectionLibrary(~G);
    StrongGenerators(~G);
    CoinvariantAlgebra(~G);
    CoinvariantAlgebra(~G`DualGroup);
    d := Dimension(G);
    G`Generators := [G.i : i in [1..Ngens(G)]];
    N := 2*d + #G`Generators;
    R := Codomain(c);
    H`CoveringFreeAlgebra := FreeAlgebra(R, N);
    Generators(~H`CoveringFreeAlgebra);
    One(~H`CoveringFreeAlgebra);
    Zero(~H`CoveringFreeAlgebra);
    H`NumberOfRewriteRules:=0;
    H`GeneratorDegrees := [];
    WordEmbedding(~G : NoInverse:=true);
    G`NumberOfGenerators := #Generators(G);

    //we have ygx since we are working in the OPPOSITE algebra!
    names := [];
    for i in Reverse([1..d]) do
        Append(~names, "y"*Sprint(i));
        Append(~H`GeneratorDegrees, -1);
    end for;
    for i in Reverse([1..#Generators(G)]) do
        Append(~names, "g"*Sprint(i));
        Append(~H`GeneratorDegrees, 0);
    end for;
    for i in Reverse([1..d]) do
        Append(~names, "x"*Sprint(i));
        Append(~H`GeneratorDegrees, 1);
    end for;

    AssignNames(~H`CoveringFreeAlgebra, names);

    for i:=1 to N do
        x := New(AlgCheResElt);
        x`Parent := H;
        x`CoveringFreeAlgebraElement := H`CoveringFreeAlgebra.i;
        x`IsInNormalForm := true;
        Append(~H`Generators,x);
    end for;

    H`BaseRing := R;
    H`Group := G;
    H`cParameter := c;

    H`NumberOfStrongGenerators := #H`Group`StrongGenerators;

    //one element
    x := New(AlgCheResElt);
    x`Parent := H;
    x`CoveringFreeAlgebraElement := H`CoveringFreeAlgebra`One;
    H`One := x;

    //set up rules

    //NOTE: WE ARE WORKING IN THE OPPOSITE ALGEBRA! SO EVERYTHING IS REVERSED!

    //x: Elements of V^*
    //y: Elements of V
    //g: Elements of G

    //x_jx_i -> x_ix_j for i>j
    for j:=1 to d do
        for i:=j+1 to d do
            rew := H`CoveringFreeAlgebra.GetVariableNumber(H,"x"*Sprint(i))*H`CoveringFreeAlgebra.GetVariableNumber(H,"x"*Sprint(j));
            AddRewriteRule(~H, [GetVariableNumber(H, "x"*Sprint(j)), GetVariableNumber(H, "x"*Sprint(i)) ], rew);
        end for;
    end for;

    //x_i g_j -> g_j (^g_j x_i)
    //note: action on x is dual action!
    for i:=1 to d do
        for j:=1 to #G`Generators do
            action := Eltseq( Transpose((G`Generators[j])^-1)[i] ); //need transpose inverse of g_j because of dual action
            actionelt := ArraySum([ action[k]*H`CoveringFreeAlgebra.GetVariableNumber(H, "x"*Sprint(k)) : k in [1..d] ] : ZeroElement:=Zero(H)); // ^g_j x_i
            rew := H`CoveringFreeAlgebra.GetVariableNumber(H, "g"*Sprint(j))*actionelt;
            AddRewriteRule(~H, [GetVariableNumber(H, "x"*Sprint(i)), GetVariableNumber(H, "g"*Sprint(j))], rew);
        end for;
    end for;

    //y_jy_i -> y_iy_j for i>j
    for j:=1 to d do
        for i:=j+1 to d do
            rew := H`CoveringFreeAlgebra.GetVariableNumber(H,"y"*Sprint(i))*H`CoveringFreeAlgebra.GetVariableNumber(H,"y"*Sprint(j));
            AddRewriteRule(~H, [GetVariableNumber(H, "y"*Sprint(j)), GetVariableNumber(H, "y"*Sprint(i)) ], rew);
        end for;
    end for;


    //g_jy_i -> ^(g_j^-1) y_i g_j
    for i:=1 to d do
        for j:=1 to #G`Generators do
            action := Eltseq((G`Generators[j]^-1)[i]);
            actionelt := ArraySum([ action[k]*H`CoveringFreeAlgebra.GetVariableNumber(H, "y"*Sprint(k)) : k in [1..d] ] : ZeroElement:=Zero(H)); // ^(g_j^-1) y_j (note the doube inverse!)
            rew := actionelt*H`CoveringFreeAlgebra.GetVariableNumber(H, "g"*Sprint(j));
            AddRewriteRule(~H, [GetVariableNumber(H, "g"*Sprint(j)), GetVariableNumber(H, "y"*Sprint(i))], rew);
        end for;
    end for;

    //x_iy_j -> y_jx_i + [y_j,x_i] (non-opp commutator) //NOTE: WE ARE WORKING IN THE OPPOSITE ALGEBRA!
    V := VectorSpace(G);
    for i:=1 to d do
        for j:=1 to d do
            //non-opposite comutator [y_j,x_i]
            comm := Zero(H`CoveringFreeAlgebra);
            for s in G`ReflectionLibraryFlat do
                comm := comm+CherednikCoefficient(V.j,V.i,s)*c(s`ReflectionClass)*(H!s`Element)`CoveringFreeAlgebraElement;
            end for;
            rew := H`CoveringFreeAlgebra.GetVariableNumber(H, "y"*Sprint(j))*H`CoveringFreeAlgebra.GetVariableNumber(H, "x"*Sprint(i)) + comm;
            AddRewriteRule(~H, [GetVariableNumber(H, "x"*Sprint(i)), GetVariableNumber(H, "y"*Sprint(j))], rew);
        end for;
    end for;

    H`xProjection := hom<H`CoveringFreeAlgebra->H`CoveringFreeAlgebra | [ One(H`CoveringFreeAlgebra) : i in [1..d+#H`Group`Generators]] cat [ H`CoveringFreeAlgebra.i : i in [d+#H`Group`Generators+1..#H`CoveringFreeAlgebra`Generators]]>;

    H`gProjection := hom<H`CoveringFreeAlgebra->H`CoveringFreeAlgebra | [ One(H`CoveringFreeAlgebra) : i in [1..d]] cat [ H`CoveringFreeAlgebra.i : i in [d+1..d+#H`Group`Generators]] cat [ One(H`CoveringFreeAlgebra) : i in [d+#H`Group`Generators+1..#H`CoveringFreeAlgebra`Generators]] >;

    H`yProjection := hom<H`CoveringFreeAlgebra->H`CoveringFreeAlgebra | [ H`CoveringFreeAlgebra.i : i in [1..d]] cat [ One(H`CoveringFreeAlgebra) : i in [d+1..#H`CoveringFreeAlgebra`Generators]]>;

    return H;

end intrinsic;


//==============================================================================
intrinsic RestrictedRationalCherednikAlgebra(G::GrpMat : Type:="EG") -> AlgCheRes
/*
    History:
        Sunday, October 27, 2013 16:51:05: Initial.
*/
{}
    c := CherednikParameter(G:Type:=Type);
    return RestrictedRationalCherednikAlgebra(G, c);

end intrinsic;

//==============================================================================
intrinsic Print(H::AlgCheRes)
/*
    History:
        Friday, September 20, 2013 16:27:50: Initial.
*/
{}

    printf "Restricted rational Cherednik algebra\n", H`BaseRing;
    printf "Generators:\n";
    IndentPush();
    str := "";
    for i:=1 to Ngens(H) do
        str *:= Sprint(H.i);
        if i lt Ngens(H) then
            str *:= ", ";
        end if;
    end for;
    printf "%o\n", str;
    IndentPop();
    printf "Generator degrees:\n";
    IndentPush();
    str := "";
    for i:=1 to Ngens(H) do
        str *:= Sprint(H`GeneratorDegrees[i]);
        if i lt Ngens(H) then
            str *:= ", ";
        end if;
    end for;
    printf "%o\n", str;
    IndentPop();
    printf "Base ring:\n";
    IndentPush();
    printf "%o\n", H`BaseRing;
    IndentPop();
    printf "Number of elementary rewrite rules:\n";
    IndentPush();
    printf "%o\n", #H`RewriteRulesDomain;
    IndentPop();
    printf "Group:\n";
    IndentPush();
    printf "%o\n", H`Group;
    IndentPop();
    printf "c-parameter:\n";
    IndentPush();
    printf "%o", H`cParameter;
    IndentPop();

end intrinsic;

//==============================================================================
intrinsic Print(x::AlgCheResElt)
/*
    History:
        Friday, September 20, 2013 16:27:50: Initial.
*/
{}

    printf "%o", x`CoveringFreeAlgebraElement;

end intrinsic;

//==============================================================================
/*
    History:
        * Friday, February 21, 2014 at 16:48:41: Moved old code to IsCoercible to be able to use !.
        * Sunday, October 27, 2013 17:07:21: Initial.
*/
intrinsic IsCoercible(H::AlgCheRes, g::GrpMatElt) -> BoolElt
{}
    if H`Group ne Parent(g) then
        return false;
    end if;

    w := H`Group`WordEmbedding(g);
    d := Dimension(H`Group);
    ginalg := H`CoveringFreeAlgebra`One;
    for i:=1 to #w do
        ginalg := ginalg * H`CoveringFreeAlgebra`Generators[d+H`Group`NumberOfGenerators+1-w[i]];
    end for;

    g := New(AlgCheResElt);
    g`CoveringFreeAlgebraElement := ginalg;
    g`Parent := H;
    return true, g;

end intrinsic;

//==============================================================================
intrinsic GetVariableNumber(H::AlgCheRes, var::MonStgElt) -> RngIntElt
{}

    part := var[1];
    num := StringToInteger(var[2]);

    if part eq "x" then
        return Dimension(H`Group)+H`Group`NumberOfGenerators+Dimension(H`Group)-num+1;
    elif part eq "g" then
        return Dimension(H`Group)+H`Group`NumberOfGenerators-num+1;
    elif part eq "y" then
        return Dimension(H`Group)-num+1;
    end if;

end intrinsic;

//==============================================================================
intrinsic GetVariable(H::AlgCheRes, var::MonStgElt) -> RngIntElt
{}

    return H.GetVariableNumber(H,var);

end intrinsic;



//==============================================================================
intrinsic Generators(H::AlgCheRes) -> SeqEnum
{}

    return H`Generators;

end intrinsic;

//==============================================================================
intrinsic AddRewriteRule(~H::AlgCheRes, code::SeqEnum, x::AlgFrElt)
/*
    History:
        Sunday, October 27, 2013 16:43:48: Initial.
*/
{}

    Append(~H`RewriteRulesDomain, code);
    Append(~H`RewriteRulesCodomain, x);
    H`NumberOfRewriteRules+:=1;

end intrinsic;

//==============================================================================
intrinsic BaseRing(H::AlgCheRes) -> Rng
{}

    return BaseRing(H`CoveringFreeAlgebra);

end intrinsic;

//==============================================================================
intrinsic Ngens(H::AlgCheRes) -> RngIntElt
/*
    History:
        Friday, September 20, 2013 16:27:31: Initial.
*/
{}

    return Ngens(H`CoveringFreeAlgebra);

end intrinsic;

//==============================================================================
intrinsic Zero(H::AlgCheRes) -> AlgCheResElt
/*
    History:
        Sunday, October 27, 2013 16:14:05: Initial.
*/
{}

    x := New(AlgCheResElt);
    x`Parent := H;
    x`CoveringFreeAlgebraElement := Zero(H`CoveringFreeAlgebra);

    return x;

end intrinsic;

//==============================================================================
intrinsic One(H::AlgCheRes) -> AlgCheResElt
/*
    History:
        Sunday, October 27, 2013 16:14:05: Initial.
*/
{}

    return H`One;

end intrinsic;

//==============================================================================
intrinsic Parent(x::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:29:07: Initial.
*/
{}

    return x`Parent;

end intrinsic;

//==============================================================================
intrinsic '.'(H::AlgCheRes, i::RngIntElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:30:08: Initial.
*/
{}

    return H`Generators[i];

end intrinsic;

//==============================================================================
intrinsic '*'(x::AlgCheResElt, y::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}
    H := x`Parent;

    pos := Position(H`ProductsDomain, <x`CoveringFreeAlgebraElement,y`CoveringFreeAlgebraElement>);
    if pos ne 0 then
        return H`ProductsCodomain[pos];
    end if;
    prod := Normalize(H,x`CoveringFreeAlgebraElement*y`CoveringFreeAlgebraElement);
    if H`KeepProducts then
        Append(~H`ProductsDomain, <x`CoveringFreeAlgebraElement,y`CoveringFreeAlgebraElement>);
        Append(~H`ProductsCodomain, prod);
    end if;
    return prod;

end intrinsic;

//==============================================================================
intrinsic '*'(x::RngElt, y::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := y`Parent;
    z`CoveringFreeAlgebraElement := x*y`CoveringFreeAlgebraElement;
    return z;

end intrinsic;

//==============================================================================
intrinsic '*'(x::AlgCheResElt, y::RngElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement*y;
    return z;

end intrinsic;

//==============================================================================
intrinsic '+'(x::AlgCheResElt, y::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y`CoveringFreeAlgebraElement;
    return z;

end intrinsic;


//==============================================================================
intrinsic '+'(x::AlgCheResElt, y::RngElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y;
    return z;

end intrinsic;

//==============================================================================
intrinsic '+'(y::RngElt, x::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgCheResElt, y::RngElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(y::RngElt, x::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgCheResElt, y::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y`CoveringFreeAlgebraElement;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgCheResElt) -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgCheResElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := -x`CoveringFreeAlgebraElement;
    return z;

end intrinsic;

//==============================================================================
intrinsic '^'(x::AlgCheResElt, i::RngIntElt : Method:="Quick") -> AlgCheResElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}
    H := Parent(x);

    //in rewrite systems the following is probalby much quicker
    if Method eq "Quick" then
        pow := H`One;
        for j:=1 to i do
            pow := pow*x;
        end for;
        return pow;
    else
        //this is the standard power
        pow := Normalize(H,x`CoveringFreeAlgebraElement^i);
        return pow;
    end if;

end intrinsic;


//==============================================================================
intrinsic Monomials(x::AlgCheResElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:38:22: Initial.
*/
{}

    mons := Monomials(x`CoveringFreeAlgebraElement);
    monscorrect := [* *];
    for mon in mons do
        y := New(AlgCheResElt);
        y`Parent := x`Parent;
        y`CoveringFreeAlgebraElement := mon;
        Append(~monscorrect, y);
    end for;

    return monscorrect;

end intrinsic;

//==============================================================================
intrinsic Terms(x::AlgCheResElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:38:22: Initial.
*/
{}

    terms := Terms(x`CoveringFreeAlgebraElement);
    termscorrect := [* *];
    for term in terms do
        y := New(AlgCheResElt);
        y`Parent := x`Parent;
        y`CoveringFreeAlgebraElement := term;
        Append(~termscorrect, y);
    end for;

    return termscorrect;

end intrinsic;

//==============================================================================
intrinsic Coefficients(x::AlgCheResElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:47:12: Initial
*/
{}

    return Coefficients(x`CoveringFreeAlgebraElement);

end intrinsic;


//==============================================================================
function ApplyRuleToTerm(H, x)
/*
    History:
        Friday, September 20, 2013 16:53:13: Initial.
*/

    mon := LeadingMonomial(x);
    coeff := LeadingCoefficient(x);
    code := Eltseq(mon);

    //try to find an elementary rule
    foundrule := false;
    if not foundrule then
        for i:=1 to H`NumberOfRewriteRules do
            subs, startcode := IsSubsequenceExtended(H`RewriteRulesDomain[i], code);
            if subs then
                endcode := startcode+#H`RewriteRulesDomain[i]-1;
                rew := H`RewriteRulesCodomain[i];
                foundrule:=true;
                break;
            end if;
        end for;
    end if;

    if foundrule then
        return coeff*RewriteMonomial(mon,startcode,endcode,rew), false;
    else
        return x, true; //element is in normal form here
    end if;

end function;

//==============================================================================
/*
    Description:
        Replace the part in mon starting at startrew and ending at endrew by rew
*/
function RewriteMonomial(mon, startrew, endrew, rew)

    return Subword(mon, 1, startrew-1)*rew*Subword(mon, endrew+1, Length(mon)-endrew);

end function;

//==============================================================================
/*
    History:
        * Saturday, March 1, 2014 at 18:01:29: Initial
*/
function NormalizeTerm(H, x)


    badterms := [ x ];
    y := Zero(H`CoveringFreeAlgebra);

    while not IsEmpty(badterms) do
        t := badterms[1];
        Remove(~badterms, 1);
        trewrite, normal := ApplyRuleToTerm(H,t);
        if normal then
            y := y+trewrite;
        else
            badterms cat:=Terms(trewrite);
        end if;
    end while;

    return FinalNormalization(H,y);

end function;


//==============================================================================
function Normalize(H, x)
/*
    History:
        Sunday, October 27, 2013 01:43:36: Initial.
*/

    y := Zero(H`CoveringFreeAlgebra);

    for t in Terms(x) do
        y := y+NormalizeTerm(H,t);
    end for;

    z := New(AlgCheResElt);
    z`CoveringFreeAlgebraElement := y;
    z`Parent := H;

    return z;

end function;

//==============================================================================
/*
    Algebra specific final normalization
*/
function FinalTermNormalization(H, y)

    //FINAL NORMALIZATION
    //rewrite group algebra part
    coeff := LeadingCoefficient(y);
    mon := LeadingMonomial(y);

    d := Dimension(H`Group);
    //find start and end of group algebra part of mon. if there is no group algebra part, there is nothing to do.
    //group algebra part is between start and ende
    code := Eltseq(mon);
    codelength:=#code;

    //part in y
    leftpart := H`yProjection(y);

    //construct leftpart in dual coinvarint algebra
    leftpartcode := Eltseq(leftpart);
    leftpartproj := One(H`Group`DualGroup`CoinvariantAlgebra);
    for i:=1 to #leftpartcode do
        leftpartproj := leftpartproj*H`Group`DualGroup`CoinvariantAlgebra.(d-leftpartcode[i]+1);
    end for;

    //reconstruct leftpartproj in H (note that leftpartproj does not have to be a monomial any more!!)
    leftpartrew := H`CoveringFreeAlgebra`Zero;
    for t in Terms(leftpartproj) do
        tcode := Exponents(t);
        trew := H`CoveringFreeAlgebra`One;
        for i:=#tcode to 1 by -1 do
            trew := trew*H`CoveringFreeAlgebra`Generators[d-i+1]^tcode[i];
        end for;
        leftpartrew := leftpartrew + (H`BaseRing!LeadingCoefficient(t))*trew;
    end for;

   // print leftpart, leftpartproj, leftpartrew;
   // print "";

    //part in x
    rightpart := H`xProjection(y);

    //construct right part in coinvarint algebra
    rightpartcode := Eltseq(rightpart);
    rightpartproj := One(H`Group`CoinvariantAlgebra);
    for i:=1 to #rightpartcode do
        rightpartproj := rightpartproj*H`Group`CoinvariantAlgebra.(d+H`Group`NumberOfGenerators+d-rightpartcode[i]+1);
    end for;

    //reconstruct rightpartproj in H (note that rightpartproj does not have to be a monomial any more!!)
    rightpartrew := H`CoveringFreeAlgebra`Zero;
    for t in Terms(rightpartproj) do
        tcode := Exponents(t);
        trew := H`CoveringFreeAlgebra`One;
        for i:=#tcode to 1 by -1 do
            trew := trew*H`CoveringFreeAlgebra`Generators[d+H`Group`NumberOfGenerators+d-i+1]^tcode[i];
        end for;
        rightpartrew := rightpartrew + (H`BaseRing!LeadingCoefficient(t))*trew;
    end for;

    //print rightpart, rightpartproj, rightpartrew;
    //print "";

    //part in g
    //NOTE THAT ALGEBRA GENERATORS ARE REVERSED!!!
    grouppart := H`gProjection(y);
    grouppartcode := Eltseq(grouppart);
    grouppartproj := Identity(H`Group);
    offset := H`Group`NumberOfGenerators+d+1;   //for computing map between H generators to strong perm generators.
    for i:=1 to #grouppartcode do
        grouppartproj := grouppartproj*H`Group`Generators[offset-grouppartcode[i]];
    end for;

    grouppartrew := H`Group`WordEmbedding(grouppartproj); //rewrites grouppart

    //same code as in coercion for group elements into AlgCheRes
    midpart := H`CoveringFreeAlgebra`One;
    for i:=1 to #grouppartrew do
        midpart := midpart*H`CoveringFreeAlgebra`Generators[offset-grouppartrew[i]];
    end for;

    return coeff*leftpartrew*midpart*rightpartrew;

end function;

//==============================================================================
function FinalNormalization(H,y)

    z := Zero(H`CoveringFreeAlgebra);
    for t in Terms(y) do
        z := z+FinalTermNormalization(H,t);
    end for;
    return z;

end function;

//==============================================================================
intrinsic Commutator(x::AlgCheResElt, y::AlgCheResElt) -> AlgCheResElt
{}

    return x*y-y*x;

end intrinsic;


//==============================================================================
intrinsic 'eq'(x::AlgCheResElt, y::AlgCheResElt) -> BoolElt
/*
    History:
        Sunday, October 27, 2013 21:37:27: Initial.
*/
{}

    return x`CoveringFreeAlgebraElement eq y`CoveringFreeAlgebraElement;

end intrinsic;

//==============================================================================
intrinsic IsCentral(x::AlgCheResElt) -> BoolElt
/*
    History:
        Sunday, October 27, 2013 21:36:17: Initial.
*/
{}
    H := Parent(x);
    for i:=1 to Ngens(H) do
        for j:=1 to Ngens(H) do
            if not H.i*x eq x*H.i then
                return false;
            end if;
        end for;
    end for;

    return true;

end intrinsic;

//==============================================================================
intrinsic EulerElement(H::AlgCheRes) -> AlgCheResElt
/*
    As in [BR13], 4.4

    History:
        Sunday, October 27, 2013 19:14:32: Initial
*/
{}
    d := Dimension(H`Group);

    eu := Zero(H);

    for i:=1 to d do
        eu := eu+H.GetVariableNumber(H, "x"*Sprint(i))*H.GetVariableNumber(H, "y"*Sprint(i));
    end for;


    for s in H`Group`ReflectionLibraryFlat do
        eu := eu+H`BaseRing!(1/(Determinant(s`Element)-1))*H`cParameter(s`ReflectionClass)*(H!s`Element);
    end for;

    return eu;

end intrinsic;




//==============================================================================
intrinsic Basis(~H::AlgCheRes)
{}
    CoinvariantAlgebra(~H`Group);
    Basis(~H`Group`CoinvariantAlgebra);
    CoinvariantAlgebra(~H`Group`DualGroup);
    Basis(~H`Group`DualGroup`CoinvariantAlgebra);
    NumberingMap(~H`Group);
    H`Basis := [* *];
    d := Dimension(H`Group);
    for i:=1 to #H`Group`CoinvariantAlgebra`Basis do
        xcode := Exponents(H`Group`CoinvariantAlgebra`Basis[i]);
        Hx := &*[ H.GetVariableNumber(H, "x"*Sprint(l))^xcode[l] : l in [1..d]];
        for j:=1 to #H`Group do
            g := H`Group`InverseNumberingMap(j);
            Hg := H!g;
            for k:=1 to #H`Group`DualGroup`CoinvariantAlgebra`Basis do
                ycode := Exponents(H`Group`DualGroup`CoinvariantAlgebra`Basis[k]);
                Hy := &*[ H.GetVariableNumber(H, "y"*Sprint(l))^ycode[l] : l in [1..d]];
                Append(~H`Basis, Hy*Hg*Hx);
            end for;
        end for;
    end for;

    H`VectorSpace := VectorSpace(BaseRing(H), #H`Basis);
    H`VectorSpaceMap := func<x|VectorSpaceMap(H,x)>;

end intrinsic;

//==============================================================================
function VectorSpaceMap(H,x)

    v:=Zero(H`VectorSpace);
    for t in Terms(x`CoveringFreeAlgebraElement) do
        mon := LeadingMonomial(t);
        coeff := LeadingCoefficient(t);
        monH := New(AlgCheResElt);
        monH`Parent := H;
        monH`CoveringFreeAlgebraElement := mon;
        pos := Position(H`Basis, monH);
        v := v + coeff*H`VectorSpace.pos;
    end for;

    return v;

end function;

//==============================================================================
intrinsic StructureConstantAlgebra(~H::AlgCheRes)
{}

    if assigned H`StructureConstantAlgebra then
        return;
    end if;

    Basis(~H);

    consts := [];
    for i:=1 to #H`Basis do
        PrintPercentage(i, #H`Basis);
        for j:=1 to #H`Basis do
            prod := H`Basis[i]*H`Basis[j];
            prodv := H`VectorSpaceMap(prod);
            for k in Support(prodv) do
                Append(~consts, <i,j,k,prodv[k]>);
            end for;
        end for;
    end for;

    H`StructureConstantAlgebra := AssociativeAlgebra<BaseRing(H), #H`Basis | consts : Rep:="Sparse", Check:=false>;

end intrinsic;

//==============================================================================
intrinsic MatrixAlgebra(~H::AlgCheRes)
{}

    Basis(~H);
    mats := [ ZeroMatrix(BaseRing(H),#H`Basis,#H`Basis) : i in [1..#Generators(H)]] ;
    for i:=1 to #Generators(H) do
        print "Computing action of generator "*Sprint(i);
        for j:=1 to #H`Basis do
            PrintPercentage(j, #H`Basis);
            v := H`VectorSpaceMap(H`Basis[j]*H`Generators[i]);
            mats[i][j] := v;
        end for;
    end for;
    H`MatrixAlgebra := MatrixAlgebra<BaseRing(H),#H`Basis | mats>;

end intrinsic;



//==============================================================================
//==============================================================================


declare type AlgCheRes[AlgCheResElt];

declare attributes AlgCheRes:
    BaseRing,
    CoveringFreeAlgebra,
    Generators,
    RewriteRulesDomain,
    RewriteRulesCodomain,
    Group,
    cParameter,
    KeepProducts,
    ProductsDomain,
    ProductsCodomain,
    Basis,
    VectorSpace,
    VectorSpaceMap,
    One,
    xProjection,
    gProjection,
    yProjection,
    RewriteRulesCounter,
    NumberOfStrongGenerators,
    NumberOfRewriteRules,
    StructureConstantAlgebra,
    MatrixAlgebra,
    GeneratorDegrees;

declare attributes AlgCheResElt:
    Parent,
    CoveringFreeAlgebraElement,
    IsInNormalForm;
