/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


/*
	Intrinsics for computing decomposition matrices (always just in split case as usual).
*/


declare attributes Grp:
	DecompositionMatrix;

declare attributes ModGrp:
	Character;	//the character of a G-module

//===================================================================
intrinsic Decomposition(chi::AlgChtrElt : StringOutput:=false) -> ModTupRng
{}

	G := Group(Parent(chi));
	CharacterTable(~G);
	dec := Decomposition(G`CharacterTable, chi);
	if StringOutput then
		str := "";
		supp := Sort(SetToSequence({i : i in [1..#G`CharacterTable] | dec[i] ne 0}));
		if IsEmpty(supp) then
			str := "0";
		else
			for i:=1 to #supp do
				str *:= Sprint(dec[supp[i]])*"*"*G`CharacterNames[supp[i]];
				if i lt #supp then
					str *:= " + ";
				end if;
			end for;
		end if;
		print str;
	end if;

	return dec;

end intrinsic;

//============================================================================
/*
    Intrinsic: DecompositionMatrix

    Split p-modular Decomposition matrix of a group

    Declaration:
        :intrinsic DecompositionMatrix(~G::Grp, p::RngIntElt)
        :intrinsic DecompositionMatrix(G::Grp, p::RngIntElt) -> Mtrx

    Parameters:
       G - A group.
       p - The characteristic.

    Description:
        Computes the split +p+-modular decomposition matrix of +G+ and assigns it to the corresponding attribute of +G+. The +i+-th row of the split +p+-modular decomposition matrix consists of the multiplicities of the absolutely irreducible representations of +G+ in characteristic +p+ in the +p+-modular reduction of the +i+-th absolutely irreducible representations of +G+ in characteristic 0.

    History:
        * Tuesday, April 1, 2014 at 17:54:34: Using DecompositionInGrothendieckGroup now.
        * Saturday, September 21, 2013 15:08:31: Initial.

*/
intrinsic DecompositionMatrix(~G::Grp, p::RngIntElt)
{The split decomposition matrix of G in characteristic p.}

    if not assigned G`DecompositionMatrix then
        G`DecompositionMatrix := AssociativeArray(Integers()); //indexed by characteristic
    end if;

    Modules(~G, 0);
    Modules(~G, p);

    A := ZeroMatrix(Integers(), #G`Modules[0], #G`Modules[p]);

    for i:=1 to #G`Modules[0] do
        Vspec := Specialize(G`Modules[0][i], p); //the p-modular reduction of G`Modules[0][i]
        A[i] := DecompositionInGrothendieckGroup(Vspec);
    end for;

    G`DecompositionMatrix[p] := A;

end intrinsic;

//============================================================================
intrinsic DecompositionMatrix(G::GrpMat, p::RngIntElt) -> ModMatRngElt
{}

    DecompositionMatrix(~G,p);
    return G`DecompositionMatrix[p];

end intrinsic;


//============================================================================
intrinsic DecompositionInGrothendieckGroup(M::ModGrp : UseCharacters:=true) -> ModTupRngElt
{
	The decomposition of M in the split Grothendieck group of the underlying group (in the characteristic of the base ring of M). In characteristic zero, it is more efficient to use the character of M.
}

    G:=Group(M);
    p:=Characteristic(BaseRing(M));
    Modules(~G,p);
    V := RSpace(Integers(), #G`Modules[p]);

    if p ne 0 or (p eq 0 and not UseCharacters) then
        dec := Zero(V);
        K := SplittingField(G,p);
        K:=CommonOverfield([*K*] cat [*BaseRing(M)*]);
        C:=ConstituentsWithMultiplicities(ChangeRing(M,K));
        K:=CommonOverfield([*K*] cat [*BaseRing(C[i][1]) : i in [1..#C]*] cat [*BaseRing(G`Modules[p][i]) : i in [1..#G`Modules[p]]*]);
        cands := {1..#G`Modules[p]};
        for i:=1 to #C do
            identified:=false;
            for j in cands do
                if IsIsomorphic(ChangeRing(C[i][1],K), ChangeRing(G`Modules[p][j], K)) then
                    dec[j] := C[i][2];
                    cands diff:={j};
                    identified:=true;
                    break;
                end if;
            end for;
            if not identified then
                error "Cannot identify some constituent.";
            end if;
        end for;
        return dec;
    else
        CharacterTable(~G);
        Character(~M);
        D:=Decomposition(G`CharacterTable, M`Character);
        return V!D;
    end if;

end intrinsic;


//============================================================================
intrinsic CharacterFixed(M::ModGrp) -> AlgChtrElt
{
	Magma's Character intrinsic is buggy over big fields (MakeCyclotmic errors). This workaround should fix this.
}

    if Characteristic(BaseRing(M)) ne 0 then
        error "Base ring of module has to be of characteristic zero.";
    end if;
    T := Type(BaseRing(M));
    if T eq FldRat or T eq FldCyc or T eq FldNum then
        return Character(M);
    end if;
    G := Group(M);
    Classes(~G);
    CharacterRing(~G);
    K := SplittingField(G,0);
    chi := [ Zero(K) : i in [1..#G`Classes]];
    for i:=1 to #G`Classes do
        w := ElementToWord(ClassRepresentative(G,i) : NoInverse:=true);
        matprod := IdentityMatrix(BaseRing(M), Dimension(M));
        for j:=1 to #w do
            matprod := matprod * ActionGenerator(M,w[j]);
        end for;
        chi[i] := K!Trace(matprod);
    end for;

    return G`CharacterRing!chi;

end intrinsic;

//============================================================================
intrinsic CharacterFixed(~M::ModGrp)
{Attaches the character.}

	if assigned M`Character then
		return;
	end if;

	M`Character := CharacterFixed(M);

end intrinsic;

//============================================================================
intrinsic Character(~M::ModGrp)
{Attaches the character.}

	CharacterFixed(~M);

end intrinsic;
