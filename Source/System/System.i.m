/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple extensions for I/O, system, and Magma operations.
*/



//============================================================================
intrinsic StartProfile()
{Starts profiling.}

    ProfileReset();
    SetProfile(true);

end intrinsic;

//============================================================================
intrinsic StopProfile() -> GrphDir
{Stops profiling and returns the profile graph.}

    SetProfile(false);
    return ProfileGraph();

end intrinsic;

//============================================================================
intrinsic HTMLProfile(P::GrphDir)
{}

	date := ShortDate();
	dir := CHAMP_GetDir()*"Profiles/"*date;
	ret := System("mkdir -p "*dir);
	ProfileHTMLOutput(P, "profile" : Dir:=dir);
	ret := System("cp "*CHAMP_GetDir()*"Profiles/sorttable.js "*dir);
	ret := System("cp "*CHAMP_GetDir()*"Profiles/leaf.gif "*dir);

end intrinsic;

//============================================================================
intrinsic Date( : Unix:=false) -> MonStgElt
{The current date and time.}

	if CHAMP_GetOS() eq "Win" then
		if Unix eq false then
			pipe := POpen(CHAMP_GetDir()*"Win\\date.exe", "r");
		else
			pipe := POpen(CHAMP_GetDir()*"Win\\date.exe +%s", "r");
		end if;
	else
		if Unix eq false then
			pipe := POpen("date", "r");
		else
			pipe := POpen("date +%s", "r");
		end if;
	end if;
	return Gets(pipe);

end intrinsic;

//============================================================================
intrinsic ShortDate() -> MonStgElt
{}

	if CHAMP_GetOS() eq "Win" then
		pipe := POpen(CHAMP_GetDir()*"Win\\date.exe +%Y%m%d%H%M%S", "r");
	else
		pipe := POpen("date +%Y%m%d%H%M%S", "r");
	end if;
	return Gets(pipe);

end intrinsic;

//============================================================================
intrinsic HumanReadableTime(t::RngIntElt) -> MonStgElt
{Prints a time in seconds in human readable format.}

	return HumanReadableTime(t*1.0);

end intrinsic;

//============================================================================
intrinsic HumanReadableTime(t::FldReElt) -> MonStgElt
{Prints a time in seconds in human readable format.}
	F := RealField(2);
	t := F!t;
	if t lt 60 then
		return Sprint(t)*"s";
	elif t lt 3600 then
		return Sprint(t/60)*"m";
	elif t lt 86400 then
		return Sprint(t/3600)*"h";
	elif t lt 604800 then
		return Sprint(t/86400)*"d";
	elif t lt 2629800 then
		return Sprint(t/604800)*"w";
	elif t lt 31557600 then
		return Sprint(t/2629800)*"M";
	else
		return Sprint(t)*"s";
	end if;

end intrinsic;

//===========================================================================
intrinsic PrintAndDelete(str::MonStgElt)
{}

	delstr := "";
	for i:=1 to #str do
		delstr *:="\b";
	end for;

	printf str;
    printf delstr;

end intrinsic;

//===========================================================================
intrinsic GetPercentage(P::RngIntElt, G::RngIntElt) -> MonStgElt
{P/G*100 with four decimal digits as string.}

    p := ChangePrecision(RealField()!(P/G*100), 4);
    return Sprint(p)*"%";

end intrinsic;

//===========================================================================
intrinsic PrintPercentage(P::RngIntElt, G::RngIntElt)
{Prints P/G*100 with four decimal digits and puts back cursor so that percentage can be overwritten.}

    p := ChangePrecision(RealField()!(P/G*100), 4);
    pstr := Sprint(p);
    pstr := pstr[1..Minimum({5,#pstr})];
    printf "%6.2o%%", pstr;
    printf "\b\b\b\b\b\b\b";
end intrinsic;

//============================================================================
intrinsic GetVersionString() -> MonStgElt
{Returns the Magma version as a string.}

    a,b,c := GetVersion();

    return Sprint(a)*"."*Sprint(b)*"-"*Sprint(c);

end intrinsic;


//============================================================================
intrinsic GetMemoryUsageInMB() -> RngIntElt
{Returns the memory used by Magma in mega bytes.}

    return Round( GetMemoryUsage()/1024^2 );

end intrinsic;



//============================================================================
intrinsic PrintAndDelete(s::MonStgElt)
{Print s and set back cursor to the beginning of s so that s can be overwritten.}
    printf s;
    delstring := "";
    for i:=1 to #s do
        delstring *:="\b";
    end for;
    printf delstring;

end intrinsic;


//============================================================================
intrinsic Sleep(n::RngIntElt)
{Sleeps for n seconds.}

    System("sleep "*Sprint(n));

end intrinsic;

//==============================================================================
intrinsic GetOS() -> MonStgElt
{Returns the operating system}

    str := Pipe("uname -a", "");
    return Split(str, " ")[1];

end intrinsic;


//============================================================================
intrinsic FileExists(file::MonStgElt) -> BoolElt
{Returns true if file exists, false otherwise.}

	if CHAMP_GetOS() eq "Win" then
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


//============================================================================
intrinsic DirectoryExists(dir::MonStgElt) -> BoolElt
{Returns true if directory exists, false otherwise.}

	if CHAMP_GetOS() eq "Win" then
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


//============================================================================
intrinsic BaseName(f::MonStgElt) -> MonStgElt
{Returns the base name of a file name, i.e., split off the directory}

	if CHAMP_GetOS() eq "Win" then
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

//============================================================================
intrinsic Directory(f::MonStgElt) -> MonStgElt
{Returns the base name of a file name, i.e., split off the directory}

	if CHAMP_GetOS() eq "Win" then
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

//============================================================================
intrinsic GetFileSize(F::MonStgElt) -> RngIntElt
{Returns the size of a file}

	return StringToInteger(Pipe("wc -c "*F*" | awk '{print $1}'", ""));

end intrinsic;


//============================================================================
intrinsic WriteCompressed(F::MonStgElt, X::.)
{Saves X to file F and compresses F using bzip, producing the file F.gz.}

    Write(F,X : Overwrite:=true);
    if CHAMP_GetOS() eq "Win" then
    	ret := System(CHAMP_GetDir()*"Win\\bzip2.exe --best -f \""*F*"\"");
    else
    	ret := System("bzip2 --best -f \""*F*"\"");
	end if;

end intrinsic;

//============================================================================
intrinsic ReadCompressed(F::MonStgElt) -> MonStgElt
{Decompresses the bzip compressed file F and reads the data.}

	if CHAMP_GetOS() eq "Win" then
		return Pipe(CHAMP_GetDir()*"Win\\bzip2.exe -c -d \""*F*"\"", "");
	else
		return Pipe("bzip2 -c -d \""*F*"\"", "");
	end if;

end intrinsic;

//============================================================================
intrinsic Hostname() -> MonStgElt
{The name of the host Magma is running on.}

	ret := Pipe("hostname", "");
	return ret[1..#ret-1];

end intrinsic;

//===========================================================================
intrinsic GetCPUType() -> MonStgElt
{}
	//currently for mac only
	cpu := Pipe("sysctl -n machdep.cpu.brand_string", "");
	return cpu[1..#cpu-1];

end intrinsic;

//============================================================================
intrinsic GetOSVersion() -> MonStgElt
{}
	//currently for mac only
	os := Pipe("sw_vers -productVersion", "");
	return "OS X "*os[1..#os-1];

end intrinsic;
