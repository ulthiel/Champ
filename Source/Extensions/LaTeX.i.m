/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Basic support for LaTeX output
*/

//==============================================================================
/*
    Intrinsic: Latex

    Prints LaTeX code for several types of objects.

    Declaration:
        :intrinsic Latex(n::RngIntElt) -> MonStgElt
        :intrinsic Latex(x::FldFinElt) -> MonStgElt
        :intrinsic Latex(A::AlgMatElt) -> MonStgElt
        :intrinsic Latex(A::ModMatRngElt) -> MonStgElt
        :intrinsic Latex(A::GrpMatElt) -> MonStgElt
        :intrinsic Latex(z::FldCycElt) ->MonStgElt
        :intrinsic Latex(f::RngSerLaurElt : IncludeDollars:=false) -> MonStgElt
        :intrinsic Latex(f::RngUPolElt : IncludeDollars:=false) -> MonStgElt
        :intrinsic Latex(T::Tup) -> MonStgElt
        :intrinsic Latex(S::MonStgElt) -> MonStgElt
        :intrinsic Latex(X::SeqEnum) -> MonStgElt
*/
intrinsic Latex(n::RngIntElt) -> MonStgElt
{}

    return Sprint(n);

end intrinsic;

//==============================================================================
intrinsic Latex(x::FldFinElt) -> MonStgElt
{}

    return Sprint(x);

end intrinsic;

//==============================================================================
intrinsic Latex(A::AlgMatElt) -> MonStgElt
{}

    str := "";
    n := Nrows(A);
    m := Ncols(A);
    str := "\\left( \\begin{array}{";
    for j:=1 to m do
        str *:= "c";
    end for;
    str *:= "}\n";
    for i:=1 to n do
        for j:=1 to m do
            str *:= Latex(A[i][j]);
            if j lt m then
                str *:= " & ";
            end if;
            if j eq m and i lt n then
                str *:= " \\\\\n";
            end if;
        end for;
    end for;
    str *:="\n\\end{array} \\right)";

    return str;

end intrinsic;

//==============================================================================
intrinsic Latex(A::ModMatRngElt) -> MonStgElt
{}

    str := "";
    n := Nrows(A);
    m := Ncols(A);
    str := "\\left( \\begin{array}{";
    for j:=1 to m do
        str *:= "c";
    end for;
    str *:= "}\n";
    for i:=1 to n do
        for j:=1 to m do
            str *:= Latex(A[i][j]);
            if j lt m then
                str *:= " & ";
            end if;
            if j eq m and i lt n then
                str *:= " \\\\\n";
            end if;
        end for;
    end for;
    str *:="\n\\end{array} \\right)";

    return str;

end intrinsic;

//==============================================================================
intrinsic Latex(A::GrpMatElt) -> MonStgElt
{}

    return Latex(Matrix(A));

end intrinsic;

//==============================================================================
intrinsic Latex(z::FldCycElt) ->MonStgElt
{}

    K := Parent(z);
    oldname := [ Sprint(K.1) ];
    AssignNames(~K, ["\\zeta_{"*Sprint(CyclotomicOrder(K))*"}"]);
    str := Sprint(z);
    AssignNames(~K, oldname);
    return str;

end intrinsic;

//==============================================================================
intrinsic Latex(f::RngSerLaurElt : IncludeDollars:=false) -> MonStgElt
/*
	Prints an element of a Laurent series ring.
*/
{}

     if f eq 0 then
        if IncludeDollars then
            return "$0$";
        else
            return "0";
        end if;
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
                if IncludeDollars then
                    str *:= " $";
                end if;
                str *:= " + ";
                if IncludeDollars then
                    str *:= " $";
                end if;
            end if;
            if IncludeDollars then
                str *:= " $";
            end if;
            if p ne 0 then
                if p eq 1 then
                    str *:= varname;
                else
                    str *:= varname*"^{"*Sprint(p)*"}";
                end if;
            else
                str *:= "1";
            end if;
            if IncludeDollars then
                str *:= " $";
            end if;
            continue;
        end if;
        coeffstr := Sprint(coeff);
        if coeffstr[1] eq "-" then
            if p ne Minimum(range) then
                if IncludeDollars then
                    coeffstr := "$-$ $"*Latex(-coeff);
                else
                    if coeff ne -1 then
                        coeffstr := " - "*Latex(-coeff);
                    else
                        coeffstr := " - ";
                    end if;
                end if;
            end if;
        else
            if p ne Minimum(range) then
                if IncludeDollars then
                    coeffstr := "$+$ $"*Latex(coeff);
                else
                    coeffstr := " + "*Latex(coeff);
                end if;
            else
                if IncludeDollars then
                    coeffstr := "$"*coeffstr;
                end if;
            end if;
        end if;
        if IncludeDollars then
            str *:= coeffstr*varname*"^{"*Sprint(p)*"}$ ";
        else
            if p eq 0 then
                str *:= coeffstr;
            elif p eq 1 then
                str *:= coeffstr*varname;
            else
                str *:= coeffstr*varname*"^{"*Sprint(p)*"}";
            end if;
        end if;
    end for;

    return str;

end intrinsic;

//============================================================================
intrinsic Latex(f::RngUPolElt : IncludeDollars:=false) -> MonStgElt
/*
    History:
        Monday, October 21, 2013 18:32:28: Initial.
*/
{}

    K := BaseRing(Parent(f));

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
    	if Coefficient(f,p) eq 0 then
    		continue;
    	end if;
        str := Replace(str, "\\^"*Sprint(p)*"$", "^{"*Sprint(p)*"}"); //exponent at line end
        str := Replace(str, "\\^"*Sprint(p)*" ", "^{"*Sprint(p)*"} "); //exponent followed by white space
    end for;

    return str;

end intrinsic;

//============================================================================
intrinsic Latex(f::FldFunRatUElt : IncludeDollars:=false, VarName:="") -> MonStgElt
{}

	P := Parent(Numerator(f));
	AssignNames(~P, Names(Parent(f)));

	if Denominator(f) eq 1 then
		return Latex(Numerator(f));
	else
		str := "\\frac{"*Latex(Denominator(f))*"}{"*Latex(Numerator(f))*"}";
	end if;

	if IncludeDollars then
		return "$"*str*"$";
	else
		return str;
	end if;

end intrinsic;



//============================================================================
intrinsic Latex(T::Tup) -> MonStgElt
{}

    str := "(";
    for i:=1 to #T do
        str *:= Latex(T[i]);
        if i lt #T then
            str *:= ", ";
        end if;
    end for;
    str *:= ")";
    return str;

end intrinsic;

//==============================================================================
intrinsic Latex(S::MonStgElt) -> MonStgElt
{}

    return S;

end intrinsic;

//==============================================================================
intrinsic Latex(X::SeqEnum) -> MonStgElt
{}

    str := "";
    for i:=1 to #X do
        str *:= Latex(X[i]);
        if i lt #X then
            str *:= ", ";
        end if;
    end for;

    return str;


end intrinsic;

//==============================================================================
intrinsic LatexTable(Data::SeqEnum, Rowlabels::SeqEnum, Columnlabels::SeqEnum : Caption:="", Colwidth:=[]) -> MonStgElt
/*
    Intrinsic: LatexTable

    Prints a LaTeX table

    Declaration:
        :intrinsic LatexTable(Data::SeqEnum, Rowlabels::SeqEnum, Columnlabels::SeqEnum : Caption:="", Colwidth:=[]) -> MonStgElt

    Parameters:
        Data - this is a sequence of sequences, where each sequence is one row of the table
        Rowlabels - a sequence encoding the row labels
        Columnlabels - a sequence encoding the column labels

    Options:
    	Caption - caption for the table
    	Colwidth - a sequence encoding the columns widths

*/
{}

    str := "\\begin{center}\n\\begin{longtable}{|";
    ncols := #Data[1];
    nrows := #Data;
    if not IsEmpty(Rowlabels) then
        str *:= "c||";
    end if;
    for i:=1 to ncols do
        if IsDefined(Colwidth, i) then
            str *:= "p{"*Sprint(Colwidth[i])*"cm}|";
        else
            str *:= "c|";
        end if;
    end for;
    str *:= "}\n\\hline\n";
    if not IsEmpty(Columnlabels) then
        if not IsEmpty(Rowlabels) then
            str *:= Sprint(Columnlabels[1])*"&";
        end if;
        for i:=2 to ncols+1 do
            str *:= Sprint(Columnlabels[i]);
            if i lt ncols+1 then
                str *:= " & ";
            end if;
        end for;
        str *:= "\\\\ \\hline \\hline \\endhead \n";
    end if;
    for i:=1 to nrows do
        if not IsEmpty(Rowlabels) then
            str *:= Sprint(Rowlabels[i])*" & ";
        end if;
        for j:=1 to ncols do
            str *:= Latex(Data[i][j]);
            if j lt ncols then
                str *:= " & ";
            end if;
        end for;
        str *:= "\\\\ \\hline \n";
    end for;
    if Caption ne "" then
        str *:= "\\caption{"*Sprint(Caption)*"}\n";
    end if;
    str *:= "\\end{longtable}\n\\end{center}";

    return str;

end intrinsic;
