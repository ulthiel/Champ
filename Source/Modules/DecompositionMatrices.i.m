/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Intrinsics for finite-dimensional modules over finite-dimensional algebras over fields.
*/

//==============================================================================
/*
    Intrinsic:
        DecompositionMatrix

    Description:
        Computes the decomposition of the family of modules M over some algebra A into the simple modules of A.
*/
intrinsic DecompositionMatrix(M::SeqEnum) -> Mtrx, SeqEnum
{}

    //first collect the simples
    print "Computing constituents";
    Msimples := [];
    for i:=1 to #M do
        consts := Constituents(M[i]);
        for j:=1 to #consts do
            identifiedas:=0;
            for k:=1 to #Msimples do
                if IsIsomorphicFixed(Msimples[k],consts[j]) then
                    identifiedas := k;
                    break;
                end if;
            end for;
            if identifiedas eq 0 then
                Append(~Msimples,consts[j]);
                identifiedas:=#Msimples;
            end if;
        end for;
        PrintPercentage(i, #M);
    end for;

    //now collect multiplicities
    print "Computing multiplicities";
    decs := [];
    for i:=1 to #M do
        for j:=1 to #Msimples do
            Append(~decs, <i,j,Multiplicity(M[i],Msimples[j])>);
        end for;
        PrintPercentage(i, #M);
    end for;

/*
    algebrasimples := [];
    decs:=[];
    for i:=1 to #M do
        consts := ConstituentsWithMultiplicities(M[i]);
        for j:=1 to #consts do
            identifiedas:=0;
            for k:=1 to #algebrasimples do
                if IsIsomorphicFixed(algebrasimples[k],consts[j][1]) then
                    identifiedas := k;
                    break;
                end if;
            end for;
            if identifiedas eq 0 then
                Append(~algebrasimples,consts[j][1]);
                identifiedas:=#algebrasimples;
            end if;
            //print i, identifiedas;
            Append(~decs, <i,identifiedas,consts[j][2]>); //identifiedas occurs in i with multipliciy consts[j][2]
        end for;
        PrintPercentage(i, #M);
    end for;
*/

    D := ZeroMatrix(Integers(), #M, #Msimples);

    for d in decs do
        D[d[1]][d[2]] := d[3];
    end for;

    return D, Msimples;

end intrinsic;


//==============================================================================
intrinsic DecompositionMatrix(M::SeqEnum, P::RngOrdIdl) -> AlgMatElt, SeqEnum, SeqEnum
/*
    Intrinsic:
        DecompositionMatrix

    Description:
        Let M be a family of modules over an algebra A over a Dedekind domain O and let P be a prime ideal of O. This function computes the constituents of the modules of M, the constituents of the P-modular reductions of the modules in M, and the decomposition matrix from the simples to the mod-P simples.
*/
{}
    Msimples := [];
    for i:=1 to #M do
        consts := Constituents(M[i]);
        for j:=1 to #consts do
            identifiedas:=0;
            for k:=1 to #Msimples do
                if IsIsomorphicFixed(Msimples[k],consts[j]) then
                    identifiedas := k;
                    break;
                end if;
            end for;
            if identifiedas eq 0 then
                Append(~Msimples,consts[j]);
                identifiedas:=#Msimples;
            end if;
        end for;
    end for;

    MP := [ Specialize(M[i],P) : i in [1..#M] ];
    MPsimples := [];
    for i:=1 to #MP do
        consts := Constituents(MP[i]);
        for j:=1 to #consts do
            identifiedas:=0;
            for k:=1 to #MPsimples do
                if IsIsomorphicFixed(MPsimples[k],consts[j]) then
                    identifiedas := k;
                    break;
                end if;
            end for;
            if identifiedas eq 0 then
                Append(~MPsimples,consts[j]);
                identifiedas:=#MPsimples;
            end if;
        end for;
    end for;

    decs := [];

    for i:=1 to #Msimples do
        simplespec := Specialize(Msimples[i],P);
        consts := ConstituentsWithMultiplicities(simplespec);
        cands := {1..#MPsimples};
        for j:=1 to #consts do
            identifiedas := 0;
            for k in cands do
                if IsIsomorphicFixed(consts[j][1], MPsimples[k]) then
                    identifiedas := k;
                    cands diff:={k};
                    break;
                end if;
            end for;
            if identifiedas eq 0 then
                error "Could not identify constituent.";
            end if;
            Append(~decs, <i,identifiedas, consts[j][2]>);
        end for;
    end for;

    D := ZeroMatrix(Integers(), #Msimples, #MPsimples);

    for d in decs do
        D[d[1]][d[2]] := d[3];
    end for;

    return D, Msimples, MPsimples;

end intrinsic;



//============================================================================
intrinsic IsIsomorphicFixed(M::ModRng, N::ModRng) -> BoolElt
{
	There seems to be a bug with IsIsomorphic for modules in characteristic zero (tested in 2.19-8). This was fixed in 2.20-6.
}

    res := IsIsomorphic(M,N);

	a,b,c := GetVersion();
	if [a,b,c] ge [2,20,6] then
		return res;
	end if;

    if Characteristic(BaseRing(M)) ne 0 or res eq true then
        return res;
    else
        if Dimension(M) ne Dimension(N) then
            return false;
        end if;
        if IsIrreducible(M) and IsIrreducible(N) then
            H := AHom(M,N);
            return not Dimension(H) eq 0;
        else
            error "Cannot decide";
        end if;
    end if;

end intrinsic;

//==============================================================================
intrinsic Multiplicity(M::ModRng, S::ModRng) -> RngIntElt
/*
    There seems to be a bug with IsIsomorphic for modules in characteristic zero (tested in 2.19-8).
*/
{The multipliciy of the simple module S in M.}

    return Integers()!(Dimension(AHom(S,M))/Dimension(AHom(S,S)));

end intrinsic;
