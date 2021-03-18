//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Database functions.
//
//==============================================================================


//============================================================================
intrinsic CHAMP_GetDBDir() -> MonStgElt
{Returns the directory of the database of CHAMP (full path).}

	if CHAMP_GetOS() eq "Win" then
		return CHAMP_GetDir()*"DB\\";
	else
		return CHAMP_GetDir()*"DB/";
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
{Returns object with the given name in the directory dir of the CHAMP database. If the category of the created object has the attributes DBName, DBDir and DBFile, then these are set to the name, dir and to the path of the object file inside the CHAMP database directory, respectively.}

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

		ext := "";
	if FileExists(file*".o.m") then
		X := eval Read(file*".o.m");
		ext := ".o.m";
	elif FileExists(file*".o.m.bz2") then
		X := eval ReadCompressed(file*".o.m.bz2");
		ext := ".o.m.bz2";
	else
		error "File does not exist in database.";
	end if;

	if "DBDir" in GetAttributes(Category(X)) then
		X`DBDir := dir;
	end if;
	if "DBFile" in GetAttributes(Category(X)) then
		X`DBFile := dir*name*ext;
	end if;
	if "DBName" in GetAttributes(Category(X)) then
		X`DBName := name;
	end if;

	return X;

end intrinsic;


//============================================================================
intrinsic CHAMP_SaveToDB(X::MonStgElt, dir::MonStgElt, f::MonStgElt)
{Save object X with name f in the direcotory dir of the database.}

	//Create directory dir (if it not exists)
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
