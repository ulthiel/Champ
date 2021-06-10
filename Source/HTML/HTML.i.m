/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple HTML file creation.
*/

//============================================================================
intrinsic WriteHTMLStart(file::MonStgElt : Title:="")
{Writes header of a HTML file.}

	str := "<html>\n<head>\n";
	str *:= "<title>"*Title*"</title>\n";
	str *:= "<link href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
	str *:= "</head>";
	str *:= "<body>\n";

	Write(file, str : Overwrite:=true);

end intrinsic;

//============================================================================
intrinsic WriteHTMLEnd(file::MonStgElt)
{Writes end of a HTML file.}

	str := "</body>\n</html>";

	Write(file, str : Overwrite:=false);

end intrinsic;

//============================================================================
intrinsic ShowHTML(file::MonStgElt)
{}

	System("open "*file);

end intrinsic;

//============================================================================
intrinsic HTMLTable(X::SeqEnum : Header:=[], WithStartEnd:=true) -> MonStgElt
{}

	if WithStartEnd then
		str := "<table>\n";
	else
		str := "";
	end if;
	if not IsEmpty(Header) then
		str *:= "    <tr>\n";
		for j:=1 to #Header do
			str *:= "        <th>"*Sprint(Header[j])*"</th>\n";
		end for;
		str *:= "    </tr>\n";
	end if;
	for i:=1 to #X do
		str *:= "    <tr>\n";
		for j:=1 to #X[i] do
			str *:= "        <td>"*Sprint(X[i][j])*"</td>\n";
		end for;
		str *:= "    </tr>\n";
	end for;
	if WithStartEnd then
		str *:= "</table>\n";
	end if;

	return str;

end intrinsic;

//============================================================================
intrinsic HTMLMatrix(X::SeqEnum : WithStartEnd:=true) -> MonStgElt
{}

	if WithStartEnd then
		str := "<table class=\"matrix\">\n";
	else
		str := "";
	end if;
	for i:=1 to #X do
		str *:= "    <tr>\n";
		for j:=1 to #X[i] do
			str *:= "        <td>"*Sprint(X[i][j])*"</td>\n";
		end for;
		str *:= "    </tr>\n";
	end for;
	if WithStartEnd then
		str *:= "</table>\n";
	end if;

	return str;

end intrinsic;




//============================================================================
intrinsic HTML(f::RngSerLaurElt) -> MonStgElt
/*
	Prints an element of a Laurent series ring.
*/
{}

     if f eq 0 then
            return "0";
    end if;

    range := Exponents(f);
    varname := Sprint(Parent(f).1);

    str := "";
    for p in range do
        coeff := Coefficient(f,p);
        if coeff eq 0 then
            continue;
        end if;
        if coeff eq 1 then
            if p ne Minimum(range) then
                str *:= " + ";
            end if;
            if p ne 0 then
                if p eq 1 then
                    str *:= varname;
                else
                    str *:= varname*"<sup>"*Sprint(p)*"</sup>";
                end if;
            else
                str *:= "1";
            end if;
            continue;
        end if;
        coeffstr := Sprint(coeff);
        if coeffstr[1] eq "-" then
            if p ne Minimum(range) then
                if coeff ne -1 then
                    coeffstr := " - "*Latex(-coeff);
                else
                    coeffstr := " - ";
                end if;
            end if;
        else
            if p ne Minimum(range) then
                coeffstr := " + "*Latex(coeff);
            end if;
        end if;

            if p eq 0 then
                str *:= coeffstr;
            elif p eq 1 then
                str *:= coeffstr*varname;
            else
                str *:= coeffstr*varname*"<sup>"*Sprint(p)*"</sup>";
            end if;
    end for;

    return str;

end intrinsic;

intrinsic HTML(f::FldFunRatUElt) -> MonStgElt
{}

	R := Parent(f);
	S := PolynomialRing(BaseRing(R));
	AssignNames(~S, Names(R));
	if IsOne(Denominator(f)) then
		return HTML(S!Numerator(f));
	else
		return "("*Sprint(S!Numerator(f))*")/("*Sprint(S!Denominator(f))*")";
	end if;

end intrinsic;

//============================================================================
intrinsic HTML(f::RngUPolElt) -> MonStgElt
{}

    /*K := BaseRing(Parent(f));

    if Type(K) eq FldCyc then
        oldname := [ Sprint(K.1) ];
        AssignNames(~K, ["\\zeta_{"*Sprint(CyclotomicOrder(K))*"}"]);
    end if;

    str := Sprint(f);

    if Type(K) eq FldCyc then
        AssignNames(~K, oldname);
    end if;

    str := Replace(str, "\\*", "");

    //replace powers
    for p in Exponents(f) do
        str := Replace(str, "\\^"*Sprint(p), "<sup>"*Sprint(p)*"<\\/sup>");
    end for;

    return str;*/

		str := Sprint(f, "Latex");
		for p in Exponents(f) do
			str := Replace(str, "\\^"*Sprint(p), "<sup>"*Sprint(p)*"<\\/sup>");
			str := Replace(str, "\\^{"*Sprint(p)*"}", "<sup>"*Sprint(p)*"<\\/sup>");
		end for;

		//name := Names(Parent(f))[1];
		//str := Replace(str, name, "<i>"*name*"<\\/i>");

		return str;


end intrinsic;
