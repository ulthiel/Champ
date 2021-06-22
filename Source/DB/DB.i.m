/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


/*
    Database functions.

    The CHAMP database consists of files inside the directory +DB+ of the base directory of CHAMP. Each file defines some object (for example a group) and its content can be evaluated in Magma with the +eval+ command. If such a file is loaded with the intrinsic +CHAMP_LoadFromDB+, then the attributes +DBDir+ and +DBFile+ of this object are set accordingly. The philosophy is that inside the database directory +DBDir+ of this object there can be stored additional data about this object (for example character data or irreducible representations) which can later be loaded by the appropriate function by checking if this data exists in +DBDir+. This prevents heavy calculations (like the invariant degrees for the reflection group E8 or G31 which are known anyways) and allows to use data consistently (like labelings of the irreducible characters).

    So far, the database only works properly with matrix groups and associated data (this also being the main motivation for implementing the database). Realizations of specific matrix groups (like exceptional complex reflection groups) are stored inside specific directories inside the directory +GrpMat+ of the CHAMP database directory. In each such subdirectory the group is defined in the file +GrpMat.m+. Character data (see <GrpRep.i.m>) is stored inside +CharacterData.m+, explicit realizations of the irreducible representations in characteristic zero are stored in +Representations_0.m+. The invariant degrees are stored in +Degrees.m+. The following data is stored so far

    * +Gi_CHEVIE+: Exceptional complex reflection group Gi (for 4 <= i <= 37) realized as in <CHEVIE> (the matrices are always transposed due to Magma's right action). The character data (including the labeling of the characters) and the irreducible characteristic zero representations (with transposed matrices as usual) are always the *same* as in CHEVIE. Character data and invariant degrees are stored for all i, whereas the irreducible representations are only stored for i <= 31.

    * +Gi_CHEVIE_Dual+: The dual of Gi_CHEVIE. This is mainly used to store the invariant degrees of the dual groups (which are the same as for the original group) so that these do not have to be computed.

    * +Gi_Magma+: Exceptional complex reflection group Gi realized as in Magma 2.19 (for 4 <= i <= 37).

    * +Gi_Magma_Dual+: Dual of Gi_Magma.

    * +B2_Cedric+: Weyl group of type +B2+ and irreducible characteristic zero representations as realized in <BR13>.

*/


//============================================================================
intrinsic CHAMP_GetDBDir() -> MonStgElt
{Returns the directory of the database of CHAMP (full path).}

    if CHAMP_GetOS() eq "Win" then
    	return CHAMP_GetDir()*"Champ-DB\\";
    else
    	return CHAMP_GetDir()*"Champ-DB/";
    end if;

end intrinsic;

//============================================================================
intrinsic CHAMP_ExistsInDB(dir::MonStgElt, name::MonStgElt) -> BoolElt
{Checks of an object exists in the CHAMP database.}

    if CHAMP_GetOS() eq "Win" then
    	if dir[#dir] ne "\\" then
			dir *:= "\\";
		end if;
    else
		if dir[#dir] ne "/" then
			dir *:= "/";
		end if;
	end if;

    file := CHAMP_GetDBDir()*dir*name;

    return FileExists(file*".o.m") or FileExists(file*".o.m.bz2");

end intrinsic;

//============================================================================
intrinsic CHAMP_GetFromDB(dir::MonStgElt, name::MonStgElt) -> .
{Returns object defined in name.m in directory dir inside the CHAMP database. If the category of the created object has the attributes DBDir and DBFile, then these are set to dir and to the path of the object file inside the CHAMP database directory, respectively.}

    if CHAMP_GetOS() eq "Win" then
    	if dir[#dir] ne "\\" then
			dir *:= "\\";
		end if;
    else
		if dir[#dir] ne "/" then
			dir *:= "/";
		end if;
	end if;

    file := CHAMP_GetDBDir()*dir*name;

    if FileExists(file*".o.m") then
    	X := eval Read(file*".o.m");
    elif FileExists(file*".o.m.bz2") then
    	X := eval ReadCompressed(file*".o.m.bz2");
    else
    	error "File does not exist in database.";
    end if;

    if "DBDir" in GetAttributes(Category(X)) then
        X`DBDir := dir;
    end if;
    if "DBFile" in GetAttributes(Category(X)) then
        X`DBFile := dir*name*".m";
    end if;

    return X;

end intrinsic;



//==============================================================================
/*
    Intrinsic: CHAMP_FindInDB

    Tries to find an object in CHAMP's database.

    Declaration:
        :intrinsic CHAMP_FindInDB(X::.) -> MonStgElt, MonStgElt

    Parameters:
        X - Any Magma object.

    Description:
        Constructs each object in the CHAMP database and checks if it is equal to +X+. If so, its file name and directory inside the CHAMP database are returned.

    Comments:
        So far, this intrinsic works only for matrix groups as these are the only objects where we need this function at the moment.

    History:
        * Saturday, September 14, 2013 16:10:40: Initial.

*/
intrinsic CHAMP_FindInDB(X::.) -> MonStgElt, MonStgElt
{Tries to find object X in the CHAMP database and returns the directory and the name if it is found. The object can the be obtained via CHAMP_GetFromDB(dir, name).}

    if assigned X`DBFile then
        dir := Directory(X`DBFile);
        name := BaseName(X`DBFile);
        return X`DBFile, dir;
    end if;

    objtype := "";
    if Type(X) eq GrpMat then
        objtype := "GrpMat";
    else
        return "", "";
    end if;

    res := Split(Pipe("cd "*CHAMP_GetDBDir()*" && find * -name "*objtype*"*.m", ""), "\n");
    if #res eq 0 then
        return "","";
    end if;
    for f in res do
        Y := CHAMP_GetFromDB(f);
        if objtype eq "GrpMat" and Type(BaseRing(X)) eq Type(BaseRing(Y)) and BaseRing(X) eq BaseRing(Y) and Dimension(X) eq Dimension(Y) and Y eq X then
           return f, Directory(f), BaseName(f);
        end if;
    end for;

    return "","";

end intrinsic;

//==============================================================================
intrinsic CHAMP_FindInDB(~X::. : Overwrite:=false)
/*
    History:
        Saturday, September 14, 2013 17:26:39: Initial.
*/
{Sets X`DBDir to its DBDir if it is in the database.}

    if assigned X`DBDir then
        return;
    end if;

    f, dir := CHAMP_FindInDB(X);

    if f ne "" then
        X`DBDir := dir;
    end if;

end intrinsic;

//============================================================================
intrinsic CHAMP_SaveToDB(X::MonStgElt, dir::MonStgElt, f::MonStgElt)
{Save object X to object name f in the database.}

    if CHAMP_GetOS() eq "Win" then
    	Pipe(CHAMP_GetDir()*"Win\\mkdir.exe -p "*CHAMP_GetDBDir()*dir, "");
    	file := CHAMP_GetDBDir()*dir*"\\"*f*".o.m";
    else
    	Pipe("mkdir -p "*CHAMP_GetDBDir()*dir, "");
    	file := CHAMP_GetDBDir()*dir*"/"*f*".o.m";
    end if;

    if #X gt 10*1024 then
    	WriteCompressed(file, X);
	else
    	Write(file, X : Overwrite:=true);
    end if;

end intrinsic;


//=============================================================================
/*
    Namespace: Grp

    Additions to the category +Grp+.
*/

//============================================================================
declare attributes Grp:
    DBDir,
    DBFile
    ;
/*
    Attribute: DBDir

    Data base directory of the group.
*/

/*
    Attribute: DBFile

    Data base file of the group.
*/
