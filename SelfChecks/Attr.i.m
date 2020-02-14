/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/

/*
    File: SelfCheck.i.m
    
    Internal CHAMP self check intrinsics (see also <SelfCheck.m>).
*/

/*
    Intrinsic: CheckAttributeAssignment
    
    Checks if non-procedural intrinsics modify attributes. This should be the case.
    
    Declaration:
        :intrinsic CheckAttributeAssignment(X::.)
*/
intrinsic CheckAttributeAssignment(X::.)
{}
    X`MyAttribute := 0;
end intrinsic;