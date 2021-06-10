/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple intrinsics around sets (and maps of sets), sequences, lists.
*/


//==============================================================================
intrinsic ArrayProduct(X::SeqEnum : OneElement:=One(Integers())) -> .
/*
    Intrinsic: ArrayProduct

    The product of the elements in a sequence.

    Declaration:
        :intrinsic ArrayProduct(X::SeqEnum : OneElement:=One(Integers())) -> .

    Parameters:
       X - A sequence.

    Options:
        OneElement - The result for empty products. By default, it is +One(Integers())+.

    Description:
        The product of all the elements in the sequence +X+ in correct order. If +X+ is empty, then +OneElement+ is returned. This fixes Magma's &*[].

    History:
        * Friday, June 28, 2013 18:19:34: Changed empty return value to One(Integers()).
*/
{}

    if IsEmpty(X) then
        return OneElement;
    else
        return &*X;
    end if;

end intrinsic;


//==============================================================================
intrinsic ArraySum(X::SeqEnum : ZeroElement:=Zero(Integers())) -> .
/*
    Intrinsic: ArraySum

    The sum of the elements in a sequence.

    Declaration:
        :intrinsic ArraySum(X::SeqEnum : ZeroElement:=Zero(Integers())) -> .

    Parameters:
       X - A sequence.

    Options:
        ZeroElement - The result for empty sums. By default, it is +One(Integers())+.

    Description:
        The sum of all the elements in the sequence +X+. If +X+ is empty, then +ZeroElement+ is returned. This fixes Magma's &+[].

    History:
        * Friday, June 28, 2013 18:19:15: Changed empty return value to Zero(Integers()).
*/
{}

    if IsEmpty(X) then
        return ZeroElement;
    else
        return &+X;
    end if;

end intrinsic;

//=============================================================================
intrinsic RemoveDuplicates(X::SeqEnum) -> SeqEnum
/*
    Intrinsic: RemoveDuplicates

    Removes duplicates from a sequence in order

    Declaration:
        :intrinsic RemoveDuplicates(X::SeqEnum) -> SeqEnum

    Description:
        Removes duplicates from +X+ and preserves order.
*/
{Removes duplicates from list in order.}

	Y := [];
	for x in X do
		if Position(Y,x) eq 0 then
			Append(~Y,x);
		end if;
	end for;

	return Y;

end intrinsic;

//==============================================================================
intrinsic FlatFixed(X::SetIndx : Depth:=-1) -> SetIndx
/*
    Intrinsic: FlatFixed

    Flatten an indexed set, set, or sequence.

    Declaration:
        :intrinsic FlatFixed(X::SetIndx : Depth:=0) -> SetIndx
        :intrinsic FlatFixed(X::SeqEnum : Depth:=0) -> SeqEnum
        :intrinsic FlatFixed(X::SetEnum : Depth:=0) -> SetEnum

    Parameters:
        X - An indexed set.

    Options:
        Depth - Recursion depth. Default is 0.

    Description:
        Flattens +X+ with recursion depth +Depth+. If +Depth+ is set to 0, then depth is determined automatically. This also fixes a stupid Magma behavior, namely the following:
        :> Flat([[]]);
        :
        :Flat(
        :S: [ [] ]
        :)
        :In file "/Applications/Magma/package/Aggregate/SeqEnum/eseq_misc.m", line 127, column 15:
        :>>   while Type(L[1]) eq SeqEnum and d lt Depth do
        :         ^
        :Runtime error in '[]': Sequence element 1 not defined

        With this fix we get
        :> FlatFixed([[]]);
        :[]

        Our version of +Flat+ is even *faster* than Magma's version:
        :> X:=[ [ [Random({-10^8..10^8}) : i in [1..200] ] : k in [1..200]] : j in [1..200] ];
        :> time Y:=Flat(X);
        :Time: 7.580
        :> time Z:=FlatFixed(X);
        :Time: 0.600
        :> Y eq Z;
        :true

    Note:
    	Magma's bad performance has been fixed in version 2.20-6 after my inquiry.

    History:
        * Friday, May 30, 2014 at 11:29:57: Enabled Flat if Version > 2.20-6.
        * Wednesday, April 2, 2014 at 10:50:40: Renamed from Flatten to FlatFixed.
        * Wednesday, January 15, 2014 at 17:13:44: Added automatic depth detection.
        * Saturday, September 21, 2013 20:58:55: Initial; ported from old project.

*/
{}

    if IsEmpty(X) or Depth eq 0 then
        return X;
    end if;

    //automatic depth detection
    if Depth eq -1 then
        Y := X;
        while Type(Y) eq SetIndx do
            Depth +:= 1;
            if IsEmpty(Y) then
                break;
            end if;
            Y := Y[1];
        end while;
    end if;

    if Depth eq 0 then
        return X;
    end if;

	flat := {@@};
	for i:=1 to #X do
        flat join:=FlatFixed(X[i] : Depth:=Depth-1);
	end for;

	return flat;

end intrinsic;

//==============================================================================
intrinsic FlatFixed(X::SeqEnum : Depth:=-1) -> SeqEnum
{}

    V1,V2,V3 := GetVersion();
    if V1 ge 2 and V2 ge 20 and V3 ge 6 then
        return Flat(X);
    end if;

    if IsEmpty(X) or Depth eq 0 then
        return X;
    end if;

    //automatic depth detection
    if Depth eq -1 then
        Y := X;
        while Type(Y) eq SeqEnum do
            Depth +:= 1;
            if IsEmpty(Y) then
                break;
            end if;
            Y := Y[1];
        end while;
    end if;

    if Depth eq 0 then
        return X;
    end if;

	flat := [];
	for i:=1 to #X do
		flat cat:=FlatFixed(X[i] : Depth:=Depth-1);
	end for;

	return flat;

end intrinsic;

//==============================================================================
intrinsic FlatFixed(X::SetEnum : Depth:=-1) -> SetEnum
/*
    History:
        * Monday, May 19, 2014 at 12:59:28: Fixed.
*/
{}

    if IsEmpty(X) or Depth eq 0 then
        return X;
    end if;

    //automatic depth detection
    if Depth eq -1 then
        Y := X;
        while Type(Y) eq SetEnum do
            Depth +:= 1;
            if IsEmpty(Y) then
                break;
            end if;
            Y := Random(Y);
        end while;
    end if;

    if Depth eq 0 then
        return X;
    end if;

	flat := {};
	for x in X do
		flat join:=FlatFixed(x : Depth:=Depth-1);
	end for;

	return flat;

end intrinsic;



//=============================================================================
intrinsic SequenceToIndexedSet(X::SeqEnum) -> SetIndx
/*
    Intrinsic: SequenceToIndexedSet

    The indexed set defined by a sequence.

    Declaration:
        :intrinsic SequenceToIndexedSet(X::SeqEnum) -> SetIndx

    Parameters:
        X - A sequence.

    Description:
        The indexed set defined by the sequence +X+ preserving order.

    History:
        * Friday, June 28, 2013 18:48:03: Initial.

*/
{}

	return {@ X[i] : i in [1..#X] @};

end intrinsic;

//==============================================================================
intrinsic Diff(~X::SetIndx, Y::Setq : Method:="Rewrite")
/*
    Intrinsic: Diff

    Removes all elements from a sequence or indexed set contained in a sequence, indexed set or set, preserving the order.

    Declaration:
        :intrinsic Diff(~X::SetIndx, Y::Setq : Method:="Rewrite")
        :intrinsic Diff(X::SetIndx, Y::Setq : Method:="Rewrite") -> SetIndx
        :intrinsic Diff(~X::SeqEnum, Y::Setq : Method:="Rewrite")
        :intrinsic Diff(X::SeqEnum, Y::Setq : Method:="Rewrite") -> SeqEnum

    Parameters:
        X - A sequence or an indexed set.
        x - A sequence or an indexed set or a set.

    Options:
        Method - Possible choices are "Rewrite" and "Loop".

    Description:
        Removes all elements lying in +Y+ from +X+, and *preserves* the order of the elements in +X+. Magma's internal +diff+ for indexed sets *changes* the order for examle:

        : > {@2,1,3@} diff {@3@};
        :   {@ 1, 2 @}

        How stupid is this? That's why we need +Diff+.

         If +Method+ is "Rewrite" then +X+ is replaced by a new sequence/indexed set not containing elements from +Y+. If +Method+ is "Loop", then Magma's +Exlude+ function is looped until no element from +Y+ does occurr in +X+ any more. The Method "Rewrite" seems to be *much faster* although it might use *more memory* as a new sequence is constructed (which later replaces +X+). With a random sequence of 500,000 integers between -10 and 10 the running time factor between the two methods is around 140 (0.16 seconds compared to 22.45). It's idiotic that this seems to be not implemented in Magma.

    History:
        * Wednesday, January 15, 2014 at 16:34:13: Renamed to "Diff".
        * Friday, June 28, 2013 18:49:52: Initial.

*/
{}
    if Method eq "Rewrite" then
        X := {@ x : x in X | x notin Y @};
    elif Method eq "Loop" then
        for y in Y do
            while Position(X,y) ne 0 do
                Exclude(~X,y);
            end while;
        end for;
    end if;

end intrinsic;

//============================================================================
intrinsic Diff(X::SetIndx, Y::Setq : Method:="Rewrite") -> SetIndx
{}

    Diff(~X, Y : Method:=Method);
    return X;

end intrinsic;


//=============================================================================
intrinsic Diff(~X::SeqEnum, Y::Setq : Method:="Rewrite")
{}

    if Method eq "Rewrite" then
        X := [ x : x in X | x notin Y ];
    elif Method eq "Loop" then
        for y in Y do
            while Position(X,y) ne 0 do
                Exclude(~X,y);
            end while;
        end for;
    end if;

end intrinsic;

//=============================================================================
intrinsic Diff(X::SeqEnum, Y::Setq : Method:="Rewrite") -> SeqEnum
{}

    Diff(~X, Y : Method:=Method);
    return X;

end intrinsic;


//==============================================================================
intrinsic IsSubsequenceExtended(x::SeqEnum, X::SeqEnum) -> BoolElt, RngIntElt
/*
    Intrinsic: IsSubsequenceExtended

    Declaration:
        :intrinsic IsSubsequenceExtended(x::SeqEnum, X::SeqEnum) -> BoolElt, RngIntElt

    Parameters:
        x - a sequence
        X - a sequence

    Description:
        Checks if +x+ is a subsequence of +X+.
        If so, the first index of +X+ where +x+ becomes a subsequence is also returned. Magma's
        function +IsSubsequence+ does *not* return this index.
*/
{}

    if IsSubsequence(x,X) then
        n1 := #x-1;
        i := Position(X,x[1]);
        while i ne 0 do
            if X[i..i+n1] eq x then
                return true, i;
            end if;
            i := Position(X,x[1],i+1);
        end while;
    else
        return false,0;
    end if;

end intrinsic;

//==============================================================================
intrinsic NumberOfOccurrences(X::SeqEnum, x::.) -> RngIntElt
/*
    Intrinsic: NumberOfOccurrences

    Declaration:
        :intrinsic NumberOfOccurrences(X::SeqEnum, x::.) -> RngIntElt

    Parameters:
        X - a sequence
        x - an element of the universe of +X+

    Description:
        The number of times +x+ occurs in +X+.
*/
{}

    p := Position(X,x);
    count := 0;
    N := #X;
    while p ne 0 do
        count +:= 1;
        if p ge N then
            break;
        else
            p := Position(X,x,p+1);
        end if;
    end while;
    return count;

end intrinsic;


//=============================================================================
intrinsic Fibers(f::Map) -> SetEnum
/*
    Intrinsic: Fibers

    The non-empty fibers of a set map

    Declaration:
        :intrinsic Fibers(f::Map) -> SetEnum

    Parameters:
        f - A map between sets.

    Description:
        Returns the set of non-empty fibers of +f+.
*/
{The set of non-empty fibers of a map f.}

    return {{x : x in Domain(f) | f(x) eq y} : y in Codomain(f)} diff {{}};

end intrinsic;


//=============================================================================
intrinsic InverseMap(f::Map) -> Map
/*
    Intrinsic: InverseMap

    The inverse of a map.

    Declaration:
        :intrinsic InverseMap(f::Map) -> Map

    Parameters:
        f - A map of sets.

    Description:
        Sometimes Magma is too stupid to compute the inverse of a map. Although this intrinsic might not be optimal, it works at least. But note that *no* check is performed if +f+ is indeed invertible; this is up to the user.
*/
{}

    return map<Codomain(f)->Domain(f) | [ <f(x),x> : x in Domain(f) ]>;

end intrinsic;



//==============================================================================
intrinsic Last(S::SeqEnum) -> .
/*
    Intrinsic: Last

    The last element of a sequence.

    Declaration:
        :intrinsic Last(S::SeqEnum) -> .
*/
{The last element of the sequence.}

    return S[#S];

end intrinsic;

//==============================================================================
intrinsic PrimeFactors(X::SetEnum[RngIntElt]) -> SetEnum
/*
    Intrinsic: PrimeFactors

    The set of prime factors of a sequence of integers.

    Declaration:
        :intrinsic PrimeFactors(X::SetEnum[RngIntElt]) -> SetEnum
*/
{}

    return {d : d in PrimeFactors(x), x in X};

end intrinsic;

//=============================================================================
/*
    Intrinsic: IsEqual

    Checks if two records are equal.

    Declaration:
        :intrinsic IsEqual(X::Rec, Y::Rec) -> BoolElt

    Parameters:
        X - A record.
        Y - A record.

    Description:
        Checks if the two records +X+ and +Y+ are equal, i.e., if they have the same Names and all Names have the same value.

    History:
        * Wednesday, October 09, 2013 15:52:40: Initial.

*/
intrinsic IsEqual(X::Rec, Y::Rec) -> BoolElt
{}

    if Names(X) ne Names(Y) then
        return false;
    end if;

    for name in Names(X) do
        if not (Type(X``name) eq Type(Y``name) and X``name eq Y``name) then
            return false;
        end if;
    end for;

    return true;

end intrinsic;
