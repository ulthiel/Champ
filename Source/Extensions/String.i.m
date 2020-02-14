/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Simple extensions for string operations
*/

//============================================================================
intrinsic Replace(str::MonStgElt, reg::MonStgElt, rep::MonStgElt) -> MonStgElt
{Replaces all occurences of the regular expression reg by rep in the string str.}

	//Buggy: Replace("q^11 + q^1", "\\^1$", "^{1}"); does not work.
	/*
	newstr := "";
	while true do
		t,m := Regexp(reg,str);
		if t eq false then
			return newstr*str;
		else
			N := Position(str,m);
			newstr *:= str[1..N-1]*rep;
			str := str[N+#m..#str];
		end if;
	end while;
	*/

	str := Pipe("sed -e \"s/"*reg*"/"*rep*"/g\"", str);
	return str[1..#str-1];

end intrinsic;
