/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


/*
    Intrinsics around characters of groups, in particular dealing with importing such data from the data base.
*/


//==============================================================================
/*
    Verbosity (Chtr already exists!)
*/
declare verbose Characters, 5;


//============================================================================
intrinsic CharacterData(~G::Grp, Data::Tup : Check:=false)
/*
    Intrinsic: CharacterData

    Sets character data of a group.

    Declaration:
        :intrinsic CharacterData(~G::Grp, Data::Tup : Check:=false)

    Parameters:
       G - A group.
       Data - 3-Tuple with character data.

    Options:
       Check - If set to +true+, then characters are checked. Check is +false+ by default.

    Description:
        With this intrinsic the character data (conjugacy classes, irreducible complex characters, and names of these characters) can be set at once to the attributes <Grp.Classes>, <Grp.CharacterTable>, and <Grp.CharacterNames> of +G+. It is used for importing character tables from CHEVIE for example to preserve CHEVIE's labeling. The point is that although we can use Magma's +AssignClasses+ function to manually set the conjugacy classes internally, the resulting list will usually be *resorted* so that we cannot afterwards import our character table! That's why we implemented this function which does all this at the same time by computing the column permutation caused by Magma's internal resort and then rearranging our given character table.

        +G+ can be any group and +Data+ is a 3-tuple giving the character data in the following form
            * The first component is a system of representatives of the conjugacy classes of +G+ (see below).
            * The second component is the character table with respect to the given conjugacy classes.
            * The third component are names for the characters.

        The conjugacy classes (first component of +Data+) can be given in a variety of forms
            * As a sequence of elements of +G+.
            * As a sequence of integer sequences descibing presentations of elements of +G+ in its generators.
            * As a sequence of 2-tuples consisting of an element of +G+ as above and the class length of this element.
            * As a sequence of 2-tuples consisting of an integer sequence as above and the class length of the corresponding element.

        <Grp.ClassMap> is automatically assigned. If the representatives were given as words (as for example in case of the CHEVIE data), then these words are also stored in <Grp.ClassWords>. If the option +Check+ is set to +true+, then it is checked whether the given characters are indeed characters and irreducible. This option is +false+ by default.



    History:
        * Wednesday, September 25, 2013 20:25:35: Initial.

*/
{}

    //we don't allow overwrite
    if assigned G`Classes and assigned G`CharacterTable and assigned G`CharacterNames then
        return;
    end if;

    //the components of data
    C := Data[1]; //the classes
    X := Data[2]; //the character table
    N := Data[3]; //the character names

    if assigned G`Classes and not assigned G`CharacterTable and IsEmpty(X) then
        return; //no new information
    end if;

    if assigned G`Classes and assigned G`CharacterTable and not assigned G`CharacterNames and IsEmpty(N) then
        return; //no new information
    end if;

    //first, we analyze the classes
    Clist := [];
    Creps := [];
    Cwords := [];

    vprint Characters, 5: "(ClassWork). Doing class work.";

    //analyze type of C and set Clist to the interpreted list which can be used by Magma
    if ISA(ExtendedType(C), SeqEnum[GrpElt]) then
        //this is a list of group elements.
        Clist := C;
        Creps := C;
    elif ExtendedType(C) eq SeqEnum[SeqEnum[RngIntElt]] then
        //this is a list of integer sequences, corresponding to representatives of class representatives in words.
        Clist := [ WordToElement(G, w) : w in C ];
        Creps := Clist;
        Cwords := C;
    elif ExtendedType(C) eq SeqEnum[Tup] then
        if ISA(Type(Component(Universe(C),1)), Grp) then
            //this is a list of tuples of group elements and class length
            Clist := C;
            words := C;
        elif Type(Component(Universe(C),1)) eq PowSeqEnum[RngIntElt] then
            //this is a list of tuples of integer sequences and class length
            Creps := [ WordToElement(G,C[i][1]) : i in [1..#C] ];
            Clist := [ <Creps[i], C[i][2]> : i in [1..#C] ];
            Cwords := [ C[i][1] : i in [1..#C] ];
        end if;
    end if;

    if not assigned G`Classes then
        //now, we will assign the classes and, if necessary, compute the permutation sigma cause by Magma's attribute assertion.
        if Type(G) eq GrpPerm or Type(G) eq GrpMat or Type(G) eq GrpPC then

            vprint Characters, 5: "(Classes). Using attribute assertion for classes.";

            //here it comes
            AssertAttribute(G, "Classes", Clist);

            //attribute assertion might change the order so we have to check this! construct permuation sigma.
            C1 := [ClassRepresentative(G, i) : i in [1..NumberOfClasses(G)]];
            sigma := [];
            remaining := {1..#C1};
            for i:=1 to #C1 do
                j := Position(Creps, C1[i]); //perhaps its really just a permutation on elements, not on classes
                if j ne 0 then
                    Append(~sigma, j);
                    remaining diff:={j};
                    continue;
                else
                    class := Class(G, C1[i]);
                    found := false;
                    for j in remaining do
                        if Creps[j] in class then
                            Append(~sigma, j);
                            remaining diff:={j};
                            found := true;
                            break;
                        end if;
                    end for;
                    if not found then
                        error "There is no bijection.";
                    end if;
                end if;
            end for;

            if sigma ne [1..#sigma] then
                vprint Characters, 5: "(Classes). WARNING: Order of the classes in Magma is different from given list. The permutation from the internal representatives to the given ones is "*Sprint(sigma)*".";
            end if;

            if Creps ne [] then
                G`ClassWords := [ Cwords[sigma[i]] : i in [1..#sigma] ];
            end if;

        else
            vprint Characters, 5: "(Classes). Setting classes.";
            G`Classes := [ <Order(Clist[i][1]),Clist[i][1],Clist[i][2]> : i in [1..#C] ];

            if Creps ne [] then
                G`ClassWords := Creps;
            end if;

            sigma := [1..#G`Classes];
        end if;
    end if;

    //now, the character table
    if not assigned G`CharacterTable and not IsEmpty(X) then
        CharacterRing(~G);
        G`CharacterTable := [];
        for i:=1 to #X do
            chi := G`CharacterRing![ X[i][sigma[j]] : j in [1..#X[i]] ]; //we use sigma here!
            Append(~G`CharacterTable, chi);
        end for;

        if Check then
            vprint Characters, 5: "(CharacterTable). Checking character table.";
            for chi in G`CharacterTable do
                if (not IsCharacter(chi) or not IsIrreducible(chi)) then
                    error "(CharacterTable). There is a reducible character.";
                end if;
            end for;

            if #G`CharacterTable ne #G`Classes or #SequenceToSet(G`CharacterTable) ne #G`Classes then
                error "(CharacterTable). Wrong number of characters.";
            end if;
        end if;
    end if;

    //now, the character names
    if not assigned G`CharacterNames and not IsEmpty(N) then
        G`CharacterNames := N; //order was preserved
    end if;

end intrinsic;


//============================================================================
intrinsic Classes(~G::Grp : Method:="PermutationGroup", UseDB:=true, Check:=false)
/*
    Intrinsic: Classes

    Sets conjugacy classes of a group.

    Declaration:
        :intrinsic Classes(~G::Grp : Method:="PermutationGroup", UseDB:=true, Check:=false)
        and
        :intrinsic Classes(G::Grp : Method:="PermutationGroup", UseDB:=true, Check:=false)

    Parameters:
       G - A group.

    Options:
       Method - Method to be used for group types not supporting conjugacy class computation.
       UseDB - If set to +true+, then load character data from data base if available. Default is +true+.
       Check - If set to +true+, then check if the representatives are indeed a complete system of representatives of the conjugacy classes.

    Description:
        Sets <Grp.Classes> to the conjugacy classes of +G+. The group +G+ can be of any type; if there is no algorithm in Magma for the type of +G+, then the classes are computed by representing +G+ in a different type using +Method+ and then pulling the data back. If +UseDB+ is set to +true+ (default), then the character data is loaded from the data base if available and <CharacterData> is invoked so that also the character table and character names may be set already.

    History:
        * Sunday, June 30, 2013 12:19:42: Initial.

*/
{}

    if assigned G`Classes then
        vprint Characters, 5: "(Classes). Classes already set. Quitting.";
        return;
    end if;

    vprint Characters, 5: "(Classes). Computing classes.";

    //first, check the DB if asked to
    foundindb := false;
    if UseDB and assigned G`DBName then
        vprint Characters, 5: "(Classes). Querying database.";
        if CHAMP_ExistsInDB("CharacterTables", G`DBName) then
            Data := CHAMP_GetFromDB("CharacterTables", G`DBName);
            vprint Characters, 5: "(Classes). Loaded classes from database.";
            CharacterData(~G, Data:Check:=Check); //this might automatically load character table and character names!
            G`ClassMap := ClassMap(G);
            foundindb := true;
        end if;
    end if;
    if foundindb then
        return;
    end if;

    //otherwise, compute the classes
    if Type(G) eq GrpPerm or Type(G) eq GrpMat or Type(G) eq GrpPC or Type(G) eq GrpAb then
        vprint Characters, 5: "(Classes). Computing classes using Magma.";
        C := Classes(G);
        G`ClassMap := ClassMap(G);
    elif Type(G) eq GrpFP then
        if Method eq "PermutationGroup" then
            vprint Characters, 5: "(Classes). Computing classes using the classes of a permutation representation.";
            PermutationGroup(~G);
            Classes(~G`PermutationGroup);
            PushforwardAttribute(~G`FromPermutationGroup, "Classes");
            G`ClassMap := G`PermutationGroup`ClassMap*G`ToPermutationGroup;
        end if;
    end if;

    CharacterRing(~G);

end intrinsic;

//============================================================================
intrinsic Classes(G::Grp : Method:="PermutationGroup",UseDB:=true, Check:=false) -> SeqEnum
/*
    History:
        Sunday, June 30, 2013 12:24:18: Initial.
*/
{}

    Classes(~G:Method:=Method,UseDB:=UseDB, Check:=Check);
    return G`Classes;

end intrinsic;


//============================================================================
intrinsic CharacterTable(~G::Grp : UseDB:=true, Check:=false)
/*
    Intrinsic: CharacterTable

    Sets the character table of a group.

    Declaration:
        :CharacterTable(~G::Grp : UseDB:=true, Check:=false)

    Parameters:
       G - A group.

    Options:
       UseDB - If set to +true+, then load character data from data base if available. Default is +true+.
       Check - If set to +true+, then check if the characters give indeed a complete system of the irreducible complex characters.

    Description:
        Sets <CharacterTable> to the table of irreducible complex characters of +G+. The group +G+ can be of any type supporting character rings. If +UseDB+ is set to +true+ (default), then the character data is loaded from the data base if available and <CharacterData> is invoked so that also the character names may be set already.

    History:
        * Monday, August 12, 2013 11:57:42: Initial.

*/
{}
    Classes(~G : UseDB:=UseDB, Check:=Check); //this might automatically load the character table since CharacterData is called!

    if assigned G`CharacterTable then
        vprint Characters, 5: "(CharacterTable). Character table already set. Quitting.";
        return;
    end if;

    vprint Characters, 5: "(CharacterTable). Computing character table.";

    //first check if representations are assigned.
    if assigned G`Representations and IsDefined(G`Representations,0) then
        vprint Characters, 5: "(CharacterTable). Using characters from given representations.";
        CharacterRing(~G);
        G`CharacterTable := [ Character(G`Representations[0][i]) : i in [1..#G`Representations[0]] ];
        return;
    elif assigned G`Modules and IsDefined(G`Modules,0) then
        vprint Characters, 5: "(CharacterTable). Using characters from given modules.";
        CharacterRing(~G);
        G`CharacterTable := [ Character(G`Modules[0][i]) : i in [1..#G`Modules[0]] ];
        return;
    end if;

    vprint Characters, 5: "(CharacterTable). Computing character table using Magma.";
    G`CharacterTable := CharacterTable(G);
    CharacterRing(~G);

end intrinsic;

//============================================================================
intrinsic CharacterNames(~G::Grp)
/*
    Intrinsic: CharacterNames

    Sets the character names of a group.

    Declaration:
        :intrinsic CharacterNames(~G::Grp)

    Parameters:
       G - A group.

    Description:
        Invokes <Classes> with option +UseDB+ set to +true+ so that the character data is imported from the data base (if available).

    History:
        * Sunday, September 15, 2013 00:47:23: Initial.

*/
{}
    Classes(~G : UseDB:=true); //this might automatically load the character names since CharacterData is called!

    if assigned G`CharacterNames then
        vprint Characters, 5: "(CharacterNames). Character names already set. Quitting.";
        return;
    end if;

    //otherwise, no information is available (where should it come from?).

end intrinsic;


//============================================================================
intrinsic CharacterRing(~G::Grp)
/*
    Intrinsic: CharacterRing

    Sets the character ring of a group.

    Declaration:
        :intrinsic CharacterRing(~G::Grp)

    Parameters:
       G - A group.

    Description:
        Sets the attribute <Grp.CharacterRing> to the character ring of +G+ which can be of any type supporting character rings.

    History:
        * Sunday, June 30, 2013 13:26:02: Initial.

*/
{}

    if assigned G`CharacterRing then
        vprint Characters, 5: "(CharacterRing). Character ring already set. Quitting.";
        return;
    end if;

    vprint Characters, 5: "(CharacterRing). Setting up character ring.";

    Classes(~G); //perhaps classes are stored!
    G`CharacterRing := CharacterRing(G);

end intrinsic;

//============================================================================
intrinsic ClassWords(~G::Grp : Method:="FPGroup")
/*
    Intrinsic: ClassWords

    Computes word presentations for the class representatives of a group.

    Declaration:
        :intrinsic ClassWords(~G::Grp : Method:="FPGroup")

    Parameters:
       G - A group.

    Options:
        Method - Method to be used for computing words. Possible choices are "FPGroup" and "RWSGroup". Default is "FPGroup".

    Description:
        Sets the attribute <Grp.ClassWords> to a sequence containing word presentations of the class representatives of +G+.

    History:
        * Wednesday, September 25, 2013 21:16:34: Deleted accidentally.

*/
{}

    Classes(~G);

    if assigned G`ClassWords then
        return;
    end if;

    if Method eq "FPGroup" then
        FPGroup(~G);
    elif Method eq "RWSGroup" then
        RWSGroup(~G);
    end if;

    G`ClassWords := [ ElementToWord(ClassRepresentative(G,i) : Method:=Method) : i in [1..#G`Classes] ];

end intrinsic;

//============================================================================
intrinsic ShortClassWords(~G::Grp : Method:="RWSGroup")
/*
    Intrinsic: ShortClassWords

    Find representatives of the conjugacy classes of a group with short word presentation.

    Declaration:
        :intrinsic ShortClassWords(~G::Grp : Method:="RWSGroup")

    Parameters:
       G - A group.

    Options:
        Method - Method to be used for computing words. Possible choices are "FPGroup" and "RWSGroup". Default is "RWSGroup".

    Description:
        Find representatives of the conjugacy classes of +G+ with short word presentation.

    History:
        * Wednesday, September 25, 2013 21:16:34: Deleted accidentally.

*/
{}

    vprint Characters, 5: "(ShortClassWords). Trying to find short class words.";

    //ClassWords(~G:Method:=Method);
    for i:=1 to #G`Classes do
        class := Class(G,G`Classes[i][3]);
        classwords := [ ElementToWord(g : Method:=Method) : g in class ];
        minlength := Minimum({#w : w in classwords});
        shortrep := Minimum([ w : w in classwords | #w eq minlength ]);
        G`ClassWords[i] := shortrep;
    end for;

end intrinsic;
