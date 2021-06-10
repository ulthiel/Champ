/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


/*
    Simple extensions for records.
*/




//==============================================================================
intrinsic CopyRecord(~R::Rec, S::Rec : Exclude:={})
/*
    Intrinsic: CopyRecord

    Copies data from one record to another record

    Declaration:
        :intrinsic CopyRecord(~R::Rec, S::Rec : Exclude:={})

    Description:
        Copies all fields from +S+ to +R+ which also exist in +R+.

*/
{Copies all data from record S to record R if the corresponding names exist in R.}

    for name in SequenceToSet(Names(S)) diff Exclude do
        if name in Names(R) and assigned S``name then
            R``name := S``name;
        end if;
    end for;

end intrinsic;


//=============================================================================
intrinsic 'eq'(X::Rec, Y::Rec) -> BoolElt
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
*/
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
