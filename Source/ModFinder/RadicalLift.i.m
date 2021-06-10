/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Intrinsics aroud computing heads of local modules as described in my thesis.

    History:
        Saturday, October 05, 2013 10:19:19: Ported.
*/

declare verbose RadicalLift, 5;

//==============================================================================
intrinsic Specialize(M::ModGrOld, P::Tup : Rep:="Dense") -> ModGrOld
/*
    History:
        Monday, September 23, 2013 14:38:35: Initial.
*/
{Specialize M in the specialization described by P.}

    R := BaseRing(M);

    //the specialization will be theta with codomain L.
    //we construct it now.
    if #P eq 2 then
        //in this case M lives over a rational function field over a number field; P[1] are values for the indeterminates; P[2] is a prime ideal of the number field.
        K := BaseRing(R);
        theta1 := hom<R->K | P[1]>;
        L,theta2 := ResidueClassFieldGF(P[2]);
        theta := theta1*theta2;
    elif #P eq 1 then
        //in this case P[1] are values for the indeterminates or it is a prime ideal of a number field.
        if Type(P[1]) eq SeqEnum then
            L := BaseRing(R);
            theta := hom<R->L | P[1]>;
        elif Type(P[1]) eq RngOrdIdl then
            L, theta := ResidueClassFieldGF(P[1]);
        elif Type(P[1]) eq RngIntElt then
            K := BaseRing(M);
            if Type(K) eq FldRat then
                L, theta := ResidueClassFieldGF(P[1]);
            else
                m,p := RandomPrimeIdeal(K, P[1]); //random prime ideal over P[1];
                L, theta := ResidueClassFieldGF(m);
            end if;
        end if;
    end if;

    if Rep eq "Preserve" then
        Rep := M`Rep;
    end if;
    if Rep eq "Sparse" then
        Mspec := GradedModuleOld([ SparseMatrix(L, Dimension(M), Dimension(M)) : i in [1..Ngens(M)] ], M`RowDegrees, M`MatrixDegrees);
    else
        Mspec := GradedModuleOld([ ZeroMatrix(L, Dimension(M), Dimension(M)) : i in [1..Ngens(M)] ], M`RowDegrees, M`MatrixDegrees);
    end if;

    for i:=1 to Ngens(M) do
        for j:=1 to Dimension(M) do
            for k in Support(M`Matrices[i][j]) do
                Mspec`Matrices[i][j][k] := theta(M`Matrices[i][j][k]);
            end for;
        end for;
    end for;

    return Mspec;

end intrinsic;


//==============================================================================
intrinsic RandomFiniteFieldSpecialization(M::ModGrOld : ParameterRange:="Automatic", pRange:="Automatic", PreferTotallySplit:=true, pExclude:={}) -> ModGrOld, Tup
/*
    History:
        Monday, September 23, 2013 14:25:45: Initial.
*/
{Pick a random datum defining a specialization from the base ring of M to a finite field.}


    if Type(ParameterRange) eq MonStgElt and ParameterRange eq "Automatic" then
        N := NumberOfNonZeroEntries(M);
        //ParameterRange := {N..Maximum(10^7, Round(N^(1.5)))}; //the is no evidence that this is a good choice...
        ParameterRange := {N..Maximum(3*N, Minimum(10^6, Round(N^(1.5))))}; //the is no evidence that this is a good choice...
    end if;

    if Type(pRange) eq MonStgElt and pRange eq "Automatic" then
        N := NumberOfNonZeroEntries(M);
        pRange := {N..Maximum(3*N, Minimum(10^6, Round(N^(1.5))))}; //the is no evidence that this is a good choice...
    end if;

    R := BaseRing(M);
    if Type(R) eq FldFunRat then
        K := BaseRing(R);
        if Type(K) eq FldNum or Type(K) eq FldCyc or Type(K) eq FldRat then
            if Type(K) eq FldRat then
                m := RandomPrime(pRange diff pExclude);
            elif Type(K) eq FldNum or Type(K) eq FldCyc then
                m,p := RandomPrimeIdeal(K : pRange:=pRange diff pExclude, PreferTotallySplit:=PreferTotallySplit);
            end if;

            P := [ Random(ParameterRange) : i in [1..Ngens(R)] ];

            return <P,m>;

        elif Type(K) eq FldFin then
            P := [ Random(K) : i in [1..Ngens(R)] ];
            return <P>;

        else
            error "No method implemented for this type of base field.";
        end if;

    elif Type(R) eq FldNum or Type(R) eq FldCyc or Type(R) eq FldRat then

        if Type(R) eq FldRat then
            m := RandomPrime(pRange diff pExclude);
        elif Type(R) eq FldNum or Type(R) eq FldCyc then
            m,p := RandomPrimeIdeal(R : pRange:=pRange diff pExclude, PreferTotallySplit:=PreferTotallySplit);
        end if;
        return <m>;

    elif Type(R) eq FldFin then
        return <>;

    else
        error "No method implemented for this type of base field.";
    end if;

end intrinsic;

//==============================================================================
intrinsic IsRadicalLinearlyLiftable(M::ModGrOld, Mspec::ModRng : GeneratorSets:=[{}], QuotientAlg:="Own",removeredundant:=true,sortdatums:=true,updateinterval:="Automatic") -> BoolElt, ModGrOld, ModRng
/*
    History:
        Monday, September 23, 2013 14:20:34: Initial.
*/
{Computes the abstract structure of the radical of Mspec and tries to find a submodule of M with this abstract structure. Returns true, the quotient by the lifted radical (which does not have to coincide with the radical!), the lifted radical and the radical of Mspec if the method works; otherwise false is returned.}

    vprint RadicalLift, 5: "Computing radical of specialized module.";
    t := Cputime();
    Jspec := JacobsonRadical(Mspec);
    vprint RadicalLift, 5: "Time: "*Sprint(Cputime(t));
    vprint RadicalLift, 5: "Dimension: "*Sprint(Dimension(Jspec))*".";
    if Dimension(Jspec) eq 0 then
        vprint RadicalLift, 5: "Module is already irreducible.";
        return true, [**], Jspec;
    end if;
    vprint RadicalLift, 5: "Computing abstract structure of radical of specialized module.";
    A := AbstractSubspace(Jspec, Mspec);
    vprint RadicalLift, 5: "Number of variables: "*Sprint(#A`Variables)*".";
    vprint RadicalLift, 5: "Initializing submodule finder.";
    F := InitializeModuleFinder(A, M);
    vprint RadicalLift, 5: "Starting submodule finder.";
    IndentPush();
    if GeneratorSets eq [{}] then
        GeneratorSets := [ {1..#M`Matrices} ];
    end if;
    for gens in GeneratorSets do
        print "Using generators "*Sprint(gens);
        FindSubmoduleWithAbstractStructure(~F : gens:=gens,sortdatums:=sortdatums,removeredundant:=removeredundant,updateinterval:=updateinterval);
        if #F`DeterminedVariables eq #A`Variables then
            break;
        end if;
    end for;
    IndentPop();
    if #F`DeterminedVariables ne #A`Variables then
        vprint RadicalLift, 5: "Could not find a submodule with the given abstract structure.";
        return false,_,_;
    end if;
    vprint RadicalLift, 5: "Could find a subspace with the given abstract structure.";

    vprint RadicalLift, 5: "Spinning this subspace.";
    J:=Spin(M,Concretize(F));
    if #J gt Dimension(Jspec) then
        vprint RadicalLift, 5: "Subspace is not a submodule.";
        return false,_,_;
    end if;
    vprint RadicalLift, 5: "Subspace is a submodule.";
    //Q:=Quotient(M, J : Alg:=QuotientAlg, PerformSpinning:=false);
    //vprint RadicalLift, 5: "Quotient dimension: "*Sprint(Dimension(Q));

    return true, J, Jspec;

end intrinsic;

//==============================================================================
intrinsic IsRadicalLinearlyLiftable(M::ModGrOld: Rounds:=2, ParameterRange:="Automatic", pRange:="Automatic", PreferTotallySplit:=true,GeneratorSets:=[{}], QuotientAlg:="Own",removeredundant:=true,sortdatums:=true,updateinterval:="Automatic",pExclude:={}) -> BoolElt, ModGrOld, ModRng, ModRng, Tup
/*
    History:
        Monday, September 23, 2013 15:17:23: Initial.
*/
{Try to lift the radical from using IsRadicalLinearlyLiftable and a random finite field specialization.}

    for rounds:=1 to Rounds do
        P := RandomFiniteFieldSpecialization(M:ParameterRange:=ParameterRange, pRange:=pRange, pExclude:=pExclude, PreferTotallySplit:=true);
        if GetVerbose("ModGrOld") ge 5 then
            print "Trying to lift radical using finite field specialization:";
            IndentPush();
            print P;
            IndentPop();
        end if;
        IndentPush();
        Mspec := RModule(Specialize(M, P : Rep:="Dense"));
        res, J, Jspec := IsRadicalLinearlyLiftable(M, Mspec : GeneratorSets:=GeneratorSets, QuotientAlg:=QuotientAlg,removeredundant:=removeredundant,sortdatums:=sortdatums,updateinterval:=updateinterval);
        IndentPop();
        if res eq false then
            vprint RadicalLift, 5: "Attempt was not successful.";
            continue;
        else
            return true, J, Mspec, Jspec, P;
        end if;
    end for;

    vprint RadicalLift, 5: "Round limit reached without success.";
    return false,_,_,_,_;

end intrinsic;

//==============================================================================
intrinsic IsSpecializationIrreducible(M::ModGrOld, P::Tup) -> BoolElt
/*
    History:
        Monday, September 23, 2013 15:23:43: Initial.
*/
{Computes the specialization of M in the finite field specialization P and checks if M is irreducible.}

    Mspec := Specialize(M, P : Rep:="Dense");

    d := Dimension(JacobsonRadical(RModule(Mspec)));

    if d eq 0 then
        return true;
    else
        return false;
    end if;

end intrinsic;

//==============================================================================
intrinsic IsSpecializationIrreducible(M::ModGrOld : Rounds:=2, ParameterRange:="Automatic", pRange:="Automatic", pExclude:={}, PreferTotallySplit:=true) -> BoolElt
/*
    History:
        Monday, September 23, 2013 15:25:44: Initial.
*/
{}

    for rounds:=1 to Rounds do
        P := RandomFiniteFieldSpecialization(M:ParameterRange:=ParameterRange, pRange:=pRange, pExclude:=pExclude, PreferTotallySplit:=true);
        if GetVerbose("ModGrOld") ge 5 then
            print "Checking irreducibility of specialization in finite field specialization:";
            IndentPush();
            print P;
            IndentPop();
        end if;
        if IsSpecializationIrreducible(M, P) then
            return true;
        else
            continue;
        end if;
    end for;

    vprint RadicalLift, 5: "Round limit reached without success.";
    return false;


end intrinsic;

//==============================================================================
intrinsic HeadOfLocalModule(M::ModGrOld, P::Tup : GeneratorSets:=[{}],removeredundant:=true,sortdatums:=true, updateinterval:="Automatic", QuotientAlg:="Own", TopChop:=true) -> BoolElt, ModGrOld, ModGrOld
{}

    if GetVerbose("ModGrOld") ge 5 then
        print "Trying to lift radical using finite field specialization:";
        IndentPush();
        print P;
        IndentPop();
    end if;
    IndentPush();
    Mspec := RModule(Specialize(M, P : Rep:="Dense"));
    res, J, Jspec := IsRadicalLinearlyLiftable(M, Mspec : GeneratorSets:=GeneratorSets, QuotientAlg:=QuotientAlg, removeredundant:=removeredundant, sortdatums:=sortdatums,updateinterval:=updateinterval);
    IndentPop();
    if res eq false then
        if TopChop then
            vprint RadicalLift, 5: "Performing top chop.";
            U := Spin(M, [* M.Dimension(M) *] );    //works for Verma's
            if #U ne Dimension(M) then
                Q := Quotient(M,U : PerformSpinning:=true);
                res, Q, Qspec := HeadOfLocalModule(Q, P : GeneratorSets:=GeneratorSets, QuotientAlg:=QuotientAlg, removeredundant:=removeredundant, sortdatums:=sortdatums,updateinterval:=updateinterval, TopChop:=true);
                if res then
                    return true, Q, Qspec; //stupid problem. would have to lift J back to M again; I don't want to do this now...
                else
                    return false, _, _;
                end if;
            end if;
        else
            vprint RadicalLift, 5: "Attempt was not successful.";
            return false, _, _;
        end if;
    else
        Q := Quotient(M,J : PerformSpinning:=true);
        vprint RadicalLift, 5: "Specializing quotient.";
        Qspec := Specialize(Q, P : Rep:="Dense");
        RModule(~Qspec);
        vprint RadicalLift, 5: "Checking if specialized quotient is irreducible.";
        if Dimension(Qspec) ne 0 and IsIrreducible(Qspec`RModule) then //Magma has a problem if dim=0?
            vprint RadicalLift, 5: "Algorithm successful.";
            return true, Q, Qspec;
        else
            return false, _, _;
        end if;
    end if;

end intrinsic;

//==============================================================================
intrinsic HeadOfLocalModule(M::ModGrOld: Rounds:=2, ParameterRange:="Automatic", pRange:="Automatic", PreferTotallySplit:=true,GeneratorSets:=[{}],removeredundant:=true,sortdatums:=true, updateinterval:="Automatic", QuotientAlg:="Own", pExclude:={}, TopChop:=true) -> BoolElt, ModGrOld, ModGrOld, Tup
/*
    History:
        Saturday, October 05, 2013 10:36:56: Initial.
*/
{If M is local, try to compute its head (and thus its radical) using a finite field specialization}

    for rounds:=1 to Rounds do
        P := RandomFiniteFieldSpecialization(M:ParameterRange:=ParameterRange, pRange:=pRange, pExclude:=pExclude, PreferTotallySplit:=true);
        if GetVerbose("ModGrOld") ge 5 then
            print "Trying to lift radical using finite field specialization:";
            IndentPush();
            print P;
            IndentPop();
        end if;
        IndentPush();
        res, Q, Qspec := HeadOfLocalModule(M, P : GeneratorSets:=GeneratorSets, QuotientAlg:=QuotientAlg, removeredundant:=removeredundant, sortdatums:=sortdatums,updateinterval:=updateinterval, TopChop:=TopChop);
        IndentPop();
        if res eq true then
            return res, Q, Qspec, P;
        else
            vprint RadicalLift, 5: "Attempt was not successful.";
        end if;
    end for;

    vprint RadicalLift, 5: "Round limit reached without success.";
    return false,_,_,_;

end intrinsic;


//==============================================================================
intrinsic HeadsOfLocalModules(M::List : Rounds:=2, ParameterRange:="Automatic", pRange:="Automatic", PreferTotallySplit:=true, GeneratorSets:=[{}],updateinterval:="Automatic",sortdatums:=true,removeredundant:=true, L:=[**],pExclude:={}, TopChop:=true) -> BoolElt, List, Mtrx, Tup
/*
    History:
        Monday, September 23, 2013 15:47:45: Initial.

        Saturday, October 05, 2013 11:01:17: Extended.
*/
{Tries to compute the heads and decompositions of a constituent closed family of local modules having pairwise non-isomorphic heads using finite field specializations.}

    for rounds:=1 to Rounds do
        P := RandomFiniteFieldSpecialization(M[1]:ParameterRange:=ParameterRange, pRange:=pRange, pExclude:=pExclude, PreferTotallySplit:=true);
        D := ZeroMatrix(Integers(), #M, #M);
        goodP := true;
        if GetVerbose("ModGrOld") ge 5 then
            print "";
            print "Trying finite field specialization:";
            IndentPush();
            print P;
            IndentPop();
        end if;
        Mspecs := [* *];
        Lspecs := [* *];

        //compute simples if not assigned
        if L eq [* *] then
            for i:=1 to #M do
                if GetVerbose("ModGrOld") ge 5 then
                    print "";
                    print "Module "*Sprint(i)*".";
                    IndentPush();
                end if;
                res, Q, Qspec := HeadOfLocalModule(M[i], P : GeneratorSets:=GeneratorSets,sortdatums:=sortdatums,removeredundant:=removeredundant,updateinterval:=updateinterval, TopChop:=TopChop);
                if res eq false then
                    L:=[**];
                    Mspecs := [**];
                    Lspecs := [**];
                    goodP := false;
                    break;
                end if;
                Append(~Mspecs, RModule(Specialize(M[i],P : Rep:="Dense")));
                Append(~Lspecs, RModule(Qspec)); //Q is now indeed the head of M[i] and its specialization is the head of Mspec!
                Append(~L, Q);

                if GetVerbose("ModGrOld") ge 5 then
                    IndentPop();
                end if;

            end for;
        else
            Lspecs := [* *];
            for i:=1 to #M do
                Append(~Mspecs, RModule(Specialize(M[i],P : Rep:="Dense")));
                Append(~Lspecs, RModule(Specialize(L[i],P : Rep:="Dense")));
            end for;
        end if;

        if goodP then
            vprint RadicalLift, 5: "Trying to identify constituents.";
            IndentPush();
            for i:=1 to #M do
                if GetVerbose("RadicalLift") ge 5 then
                    PrintPercentage(i, #M);
                end if;
                const := ConstituentsWithMultiplicities(Mspecs[i]);
                //print const;
                cands := {1..#Lspecs};
                reconstructedconst := [];
                for j:=1 to #const do
                    jcands := [];
                    for k in cands do
                        if IsIsomorphic(Lspecs[k], const[j][1]) then
                            Append(~jcands, k);
                            cands diff:={k};
                        end if;
                    end for;
                    if #jcands gt 1 or #jcands eq 0 then
                        goodspec := false;
                        break i;
                    else
                        Append(~reconstructedconst, jcands[1]);
                    end if;
                end for;

                for j:=1 to #const do
                    D[i][reconstructedconst[j]] := const[j][2];
                end for;
            end for;
            IndentPop();
        end if;

        if goodP eq false then
            continue;
        else
            break;
        end if;

    end for;

    if not goodP then
        vprint RadicalLift, 5: "Round limit reached without success.";
        return false, _, _;
    end if;

    return true, L, D, P;

end intrinsic;
