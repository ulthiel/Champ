/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

intrinsic DaggerPartition(lambda::SeqEnum[RngIntElt], a::RngIntElt, b::RngIntElt) -> SeqEnum
{
If lambda is a partition with lambda_1 <= a and of length <= b, then lambda is contained in the box partition a^b. The dagger of lambda is the north-east reflection of the complement of lambda in this box.

Needed in Bellamy-Thiel "Cuspidal..."
}

	assert (IsEmpty(lambda) or lambda[1] le a) and #lambda le b;

	//padding
	while #lambda lt b do
		Append(~lambda, 0);
	end while;

	mu := [ #{j : j in [1..b] | a-lambda[b+1-j] ge i} : i in [1..b] ];

	//remove zeros
	while not IsEmpty(mu) and Last(mu) eq 0 do
		Remove(~mu, #mu);
	end while;

	return mu;

end intrinsic;


//===========================================================================
intrinsic Multipartitions(n::RngIntElt,m::RngIntElt) -> SeqEnum
{The m-multipartitions of n.}

    X := [];
    if m eq 2 then
        for i:=0 to n do
            X cat:= SetToSequence(Set(CartesianProduct(Partitions(i), Partitions(n-i))));
        end for;
    else
        error "Not implemented for m>2";
    end if;

    return [[x[1],x[2]] : x in X];

end intrinsic;

//============================================================================
intrinsic Bipartitions(n::RngIntElt) -> SeqEnum
{The bipartitions of n.}

	return Multipartitions(n,2);

end intrinsic;
