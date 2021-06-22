/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions to write elements as words in the generators.
*/

//===========================================================================
intrinsic ElementToWord(g::GrpElt : Method:="RWSGroup", NoInverse:=false) -> SeqEnum
/*
    Intrinsic: ElementToWord

    A group elememnt as a word in the generators.

    Declaration:
        :intrinsic ElementToWord(g::GrpElt : Method:="RWSGroup") -> SeqEnum

    Parameters:
       G - A group.

    Options:
        Method - The method to be used for computing a word presentation. Possible choices are "RWSGroup" and "FPGroup". Default is "RWSGroup".

    Description:
        The element +g+ of a group +G+ as a word in the generators of +G+.
*/
{}

    G := Parent(g);

    if Method eq "RWSGroup" then
        RWSGroup(~G);
        w := Eltseq(G`ToRWSGroup(g));
    elif Method eq "FPGroup" then
        FPGroup(~G);
        w := Eltseq(G`ToFPGroup(g));
    end if;

    if NoInverse then
        neww := [];
        for i:=1 to #w do
            if w[i] lt 0 then
                o := Order(G.Abs(w[i]));
                neww cat:= [Abs(w[i]) : j in [1..o-1] ];
            else
                neww cat:=[w[i]];
            end if;
        end for;
        w := neww;
    end if;

    return w;

end intrinsic;

//==============================================================================
intrinsic WordEmbedding(~G::Grp : Method:="RWSGroup", NoInverse:=false)
{}

    if assigned G`WordEmbedding then
        return;
    end if;

    Set(~G);
    NumberingMap(~G);

    words := [ ElementToWord(G`InverseNumberingMap(i) : NoInverse:=NoInverse) : i in [1..#G`Set]] ;

    G`WordEmbedding := map<G`Set -> words | [<G`InverseNumberingMap(i),words[i]> : i in [1..#G`Set]]>;

end intrinsic;

//==============================================================================
intrinsic WordToElement(G::Grp, w::SeqEnum) -> GrpElt
/*
    Intrinsic: WordToElement

    Construct group element from word in the generators.

    Declaration:
        :intrinsic WordToElement(G::Grp, w::SeqEnum) -> GrpElt

    Parameters:
       G - A group.
       w - A sequence describing a word in the generators of +G+.

    Description:
        If +w = [w1,w2,...,wn]+ is an integer sequence, return +G.w1*G.w2*...*G.wn+. If +wi+ is negative, then +G.wi+ is +(G.|wi|)^-1+.
*/
{}

    if Type(G) eq GrpPerm or Type(G) eq GrpMat or Type(G) eq GrpPermCox then
        return ArrayProduct([ G.(Abs(w[i]))^(Sign(w[i])) : i in [1..#w] ] : OneElement:=Identity(G) );
    elif Type(G) eq GrpFP or Type(G) eq GrpRWS then
        return G!w;
    end if;

end intrinsic;
