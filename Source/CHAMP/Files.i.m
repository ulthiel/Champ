//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// File handling
//
//==============================================================================


//##############################################################################
intrinsic CHAMP_GetUnixTool(name::MonStgElt) -> MonStgElt
{Returns the path of the (CHAMP) Unix tool that can be used under Windows.}

  if GetOSType() eq "Unix" then
    return name;
  else
    return MakePath([CHAMP_GetDir(), "Tools", "UnixTools", name*".exe"]);
  end if;

end intrinsic;

//==============================================================================
intrinsic DirectorySeparator() -> MonStgElt
{Returns / if the operating system is Unix and returns \ if it is Win.}

	if GetOSType() eq "Unix" then
		return "/";
	else
		return "\\";
	end if;
end intrinsic;

//==============================================================================
intrinsic MakePath(X::SeqEnum) -> MonStgElt
{Concetenates the elements of X to a path using the OS specific directory separator.}

	dir := "";
	for i:=1 to #X do
		dir *:= X[i];
		if i lt #X then
			dir *:= DirectorySeparator();
		end if;
	end for;

	return dir;

end intrinsic;

//==============================================================================
intrinsic FileExists(file::MonStgElt) -> BoolElt
{Returns true if file exists, false otherwise.}

	if GetOSType() eq "Win" then
		ret := System(CHAMP_GetDir()*"Win\\test.exe -f "*file)/256;
	else
		ret := System("test -f "*file)/256;
	end if;
	if ret eq 0 then
		return true;
	else
		return false;
	end if;

end intrinsic;


//==============================================================================
intrinsic DirectoryExists(dir::MonStgElt) -> BoolElt
{Returns true if directory exists, false otherwise.}

	if GetOSType() eq "Win" then
		ret := System(CHAMP_GetDir()*"Win\\test.exe -d "*dir)/256;
	else
		ret := System("test -d "*dir)/256;
	end if;
	if ret eq 0 then
		return true;
	else
		return false;
	end if;

end intrinsic;

intrinsic MakeDirectory(dir::MonStgElt)
{Creates directory with name dir.}

	if GetOSType() eq "Win" then
		Pipe(CHAMP_GetDir()*"Win\\mkdir.exe -p "*CHAMP_GetDBDir()*dir, "");
		file := CHAMP_GetDBDir()*dir*"\\"*f*".o.m";
	else
		Pipe("mkdir -p "*CHAMP_GetDBDir()*dir, "");
		file := CHAMP_GetDBDir()*dir*"/"*f*".o.m";
	end if;

end intrinsic;

//==============================================================================
intrinsic BaseName(f::MonStgElt) -> MonStgElt
{Returns the base name of a file name, i.e., split off the directory}

	if GetOSType() eq "Win" then
		pos := Position(Reverse(f), "\\"); //mathces last slash in f
	else
		pos := Position(Reverse(f), "/"); //mathces last slash in f
	end if;
	if pos eq 0 then
		name := f;
	else
		name := f[#f-pos+2..#f];
	end if;

	return name;

end intrinsic;

//==============================================================================
intrinsic Directory(f::MonStgElt) -> MonStgElt
{Returns the base name of a file name, i.e., split off the directory}

	if GetOSType() eq "Win" then
		pos := Position(Reverse(f), "\\"); //matches last slash in f
	else
		pos := Position(Reverse(f), "/"); //matches last slash in f
	end if;
	if pos eq 0 then
		dir := ".";
	else
		dir := f[1..#f-pos];
	end if;

	return dir;

end intrinsic;

//==============================================================================
intrinsic GetFileSize(F::MonStgElt) -> RngIntElt
{Returns the size of a file}

	return StringToInteger(Pipe("wc -c "*F*" | awk '{print $1}'", ""));

end intrinsic;


//==============================================================================
intrinsic WriteCompressed(F::MonStgElt, X::.)
{Saves X to file F and compresses F using bzip, producing the file F.gz.}

	Write(F,X : Overwrite:=true);
	if GetOSType() eq "Win" then
		ret := System(CHAMP_GetDir()*"Win\\bzip2.exe --best -f \""*F*"\"");
	else
		ret := System("bzip2 --best -f \""*F*"\"");
	end if;

end intrinsic;

//==============================================================================
intrinsic ReadCompressed(F::MonStgElt) -> MonStgElt
{Decompresses the bzip compressed file F and reads the data.}

	if GetOSType() eq "Win" then
		return Pipe(CHAMP_GetDir()*"Win\\bzip2.exe -c -d \""*F*"\"", "");
	else
		return Pipe("bzip2 -c -d \""*F*"\"", "");
	end if;

end intrinsic;
