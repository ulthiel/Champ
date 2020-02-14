groups := [13];



for i in groups do
    print i;

    G := ExceptionalComplexReflectionGroup(i);
    c := CherednikParameter(G : Rational:=false);
    R := Codomain(c);
    CharacterTable(~G);
    mariavars := eval Read("GAP3_Data/var"*Sprint(i));

    mariaring := PolynomialRing(Rationals(), #mariavars);

    if i eq 4 then
        phi := hom<mariaring -> R | [0, R.1, R.2] >;
    elif i eq 5 then
        phi := hom<mariaring -> R | [0, R.1, R.2, 0, R.3, R.4] >;
    elif i eq 6 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3] >;
    elif i eq 7 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, 0, R.4, R.5] >;
    elif i eq 8 then
        phi := hom<mariaring -> R | [0, R.1, R.2, R.3] >;
    elif i eq 9 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, R.4] >;
    elif i eq 10 then
        phi := hom<mariaring -> R | [0, R.1, R.2, 0, R.3, R.4, R.5] >;
    elif i eq 11 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, 0, R.4, R.5, R.6] >;
    elif i eq 12 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 13 then
        phi := hom<mariaring -> R | [0, R.2, 0, R.1] >; //wrong way probably!!
    elif i eq 14 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3] >;
    elif i eq 15 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, 0, R.4] >;
    elif i eq 16 then
        phi := hom<mariaring -> R | [0, R.1, R.2, R.3, R.4] >;
    elif i eq 17 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, R.4, R.5] >;
    elif i eq 18 then
        phi := hom<mariaring -> R | [0, R.1, R.2, 0, R.3, R.4, R.5, R.6] >;
    elif i eq 19 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3, 0, R.4, R.5, R.6, R.7] >;
    elif i eq 20 then
        phi := hom<mariaring -> R | [0, R.1, R.2] >;
    elif i eq 21 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3] >;
    elif i eq 22 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 23 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 24 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 25 then
        phi := hom<mariaring -> R | [0, R.1, R.2] >;
    elif i eq 26 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2, R.3] >;
    elif i eq 27 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 28 then
        phi := hom<mariaring -> R | [0, R.1, 0, R.2] >;
    elif i eq 29 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 30 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 31 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 32 then
        phi := hom<mariaring -> R | [0, R.1, R.2] >;
    elif i eq 33 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 34 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 35 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 36 then
        phi := hom<mariaring -> R | [0, R.1] >;
    elif i eq 37 then
        phi := hom<mariaring -> R | [0, R.1] >;
    end if;


    mariares := AssociativeArray({""});
    myres := AssociativeArray(R);

    FP := Open("GAP3_Data/Rou"*Sprint(i), "r");
    curfams := "";
    curplane := "";
    while true do
        s := Gets(FP);
        if IsEof(s) then
            stopnow := true;
        end if;
        if IsEof(s) or s eq "No essential hyperplane" or s[1] in {"a","b","c","d","e","f","g","2","3","4","5","6","7","8","9"} then //beginning of new plane
            newplane := true;
            if curfams ne "" then
                mariares[curplane] := { {Position(G`CharacterNames, Replace(y, "phi{", "\\phi_{")) : y in SequenceToSet(x)} : x in eval curfams};
            end if;
            if IsEof(s) then
                break;
            end if;
            curplane := s;
            curfams := "";
        else
            newplane := false;
            curfams *:= s;
        end if;
    end while;



    //print Keys(mariares);

    for H in Keys(mariares) do
        if H eq "No essential hyperplane" then
            Hpol := One(R);
        else
            Hstr := Replace(H, "=0", "");
            for j:=1 to #mariavars do
                Hstr := Replace(Hstr, mariavars[j], "mariaring."*Sprint(j));
            end for;
            Hstr := Replace(Hstr, "2maria", "2*maria");
            Hstr := Replace(Hstr, "3maria", "3*maria");
            Hstr := Replace(Hstr, "4maria", "4*maria");
            Hstr := Replace(Hstr, "5maria", "5*maria");
            Hpol := NormalizeRationalHyperplaneEquation(phi(eval Hstr));
        end if;

        myres[Hpol] := mariares[H];
    end for;

    myresstr := "/*\n";
    myresstr *:= "Rouquier families for exceptional complex reflection group G"*Sprint(i)*" imported from CHEVIE-V4.\n";
    myresstr *:= "Date: "*Sprint(Date())*"\n";
    myresstr *:= "Version: "*Sprint(CHAMP_GetVersion())*"\n";
    myresstr *:= "*/\n";
    myresstr *:= "P := PolynomialRing("*Sprint(BaseRing(G), "Magma")*", "*Sprint(Ngens(R))*");\n";
    myresstr *:= "AssignNames(~P,[";
    for j:=1 to Ngens(R) do
        myresstr *:= "\""*Sprint(Name(R,j))*"\"";
        if j lt Ngens(R) then
            myresstr *:= ",";
        end if;
    end for;
    myresstr *:= "]);\n";
    myresstr *:= "res := AssociativeArray(P);\n";
    AssignNames(~R, ["P."*Sprint(j) : j in [1..Ngens(R)]]);
    for H in Keys(myres) do
        myresstr *:= "res["*Sprint(H)*"] := "*Sprint(myres[H])*";\n";
    end for;
    myresstr *:= "return res";

    //res := eval myresstr;

    Write(CHAMP_GetDBDir()*"/"*G`DBDir*"/RouquierFamilies.m", myresstr : Overwrite:=true);

end for;


quit;