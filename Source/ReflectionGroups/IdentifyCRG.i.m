/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

function IdentifyIrreducibleCRG(G)

	if Order(G) eq 1 then
		return <>;
	end if;

	//first, check if isomorphic to some G(m,p,n)
	//see 15.36 in my thesis
	candidates := {<Order(G),1,1>};
	maxn := 2;
	while Factorial(maxn) le Order(G) do
		maxn +:= 1;
	end while;
	for n:=2 to maxn do
		minm := Maximum(Floor(Root(Order(G)/Factorial(n), n)), 1);
		maxm := Ceiling(Root(Order(G)/Factorial(n), n-1));
		if maxm eq 0 then
			continue;
		end if;
		for m:=minm to maxm do
			for p in Divisors(m) do
				if (m^n*Factorial(n))/p eq Order(G) then
					candidates join:={<m,p,n>};
				end if;
			end for;
		end for;
	end for;

	for cand in candidates do
		if IsIsomorphic(G, ShephardTodd(cand[1],cand[2],cand[3])) then
			return cand;
		end if;
	end for;

	//exceptionals
	order := Order(G);

	if order eq 24 then
        if IsIsomorphic(G, ShephardTodd(4)) then
        	return <4>;
        end if;
    elif order eq 48 then
        if IsIsomorphic(G, ShephardTodd(6)) then
        	return <6>;
    	elif IsIsomorphic(G, ShephardTodd(12)) then
    		return <12>;
    	end if;
    elif order eq 72 then
        if IsIsomorphic(G, ShephardTodd(5)) then
        	return <5>;
        end if;
    elif order eq 96 then
         if IsIsomorphic(G, ShephardTodd(8)) then
         	return <8>;
         elif IsIsomorphic(G, ShephardTodd(13)) then
         	return <13>;
         end if;
    elif order eq 120 then
         if IsIsomorphic(G, ShephardTodd(23)) then
         	return <23>;
         end if;
    elif order eq 144 then
        if IsIsomorphic(G, ShephardTodd(7)) then
        	return <7>;
        elif IsIsomorphic(G, ShephardTodd(14)) then
        	return <14>;
        end if;
    elif order eq 192 then
        if IsIsomorphic(G, ShephardTodd(9)) then
        	return <9>;
        end if;
    elif order eq 240 then
        if IsIsomorphic(G, ShephardTodd(22)) then
        	return <22>;
        end if;
    elif order eq 288 then
        if IsIsomorphic(G, ShephardTodd(10)) then
        	return <10>;
        elif IsIsomorphic(G, ShephardTodd(15)) then
        	return <15>;
        end if;
    elif order eq 336 then
        if IsIsomorphic(G, ShephardTodd(24)) then
        	return <24>;
        end if;
    elif order eq 360 then
        if IsIsomorphic(G, ShephardTodd(20)) then
        	return <20>;
        end if;
    elif order eq 576 then
        if IsIsomorphic(G, ShephardTodd(11)) then
        	return <11>;
        end if;
    elif order eq 600 then
        if IsIsomorphic(G, ShephardTodd(16)) then
        	return <16>;
        end if;
    elif order eq 648 then
        if IsIsomorphic(G, ShephardTodd(25)) then
        	return <25>;
        end if;
    elif order eq 720 then
        if IsIsomorphic(G, ShephardTodd(21)) then
        	return <21>;
        end if;
    elif order eq 1152 then
        if IsIsomorphic(G, ShephardTodd(28)) then
        	return <28>;
        end if;
    elif order eq 1200 then
        if IsIsomorphic(G, ShephardTodd(17)) then
        	return <17>;
        end if;
    elif order eq 1296 then
        if IsIsomorphic(G, ShephardTodd(26)) then
        	return <26>;
        end if;
    elif order eq 1800 then
        if IsIsomorphic(G, ShephardTodd(18)) then
        	return <18>;
        end if;
    elif order eq 2160 then
        if IsIsomorphic(G, ShephardTodd(27)) then
        	return <27>;
        end if;
    elif order eq 3600 then
        if IsIsomorphic(G, ShephardTodd(19)) then
        	return <19>;
        end if;
    elif order eq 7680 then
        if IsIsomorphic(G, ShephardTodd(29)) then
        	return <29>;
        end if;
    elif order eq 14400 then
        if IsIsomorphic(G, ShephardTodd(30)) then
        	return <30>;
        end if;
    elif order eq 46080 then
        if IsIsomorphic(G, ShephardTodd(31)) then
        	return <31>;
        end if;
    elif order eq 51840 then
        if IsIsomorphic(G, ShephardTodd(33)) then
        	return <33>;
        elif IsIsomorphic(G, ShephardTodd(35)) then
        	return <35>;
        end if;
    elif order eq 155520 then
        if IsIsomorphic(G, ShephardTodd(32)) then
        	return <32>;
        end if;
    elif order eq 2903040 then
        if IsIsomorphic(G, ShephardTodd(36)) then
        	return <36>;
        end if;
    elif order eq 39191040 then
        if IsIsomorphic(G, ShephardTodd(34)) then
        	return <34>;
        end if;
    elif order eq 696729600 then
        if IsIsomorphic(G, ShephardTodd(37)) then
        	return <37>;
        end if;
    end if;

    error "Could not identify group";

end function;

//============================================================================
intrinsic IdentifyCRG(G::GrpMat) -> SeqEnum
{}
	if Characteristic(BaseRing(G)) ne 0 then
		error "Only works in characteristic zero.";
	end if;
	if not IsReflectionGroup(G) then
		error "Only works for reflection groups.";
	end if;
	rho := NaturalRepresentation(G);
	if IsIrreducible(rho) then
		return [* IdentifyIrreducibleCRG(G) *];
	else
		M:=GModule(rho);
		Msummands:=DirectSumDecomposition(M);
		Gsummands := [ MatrixGroup<Dimension(N), BaseRing(G) | ActionGenerators(N)> : N in Msummands ];
		return [* IdentifyIrreducibleCRG(H) : H in Gsummands *];
	end if;

end intrinsic;

//============================================================================
intrinsic ParabolicSubgroupsOfCRG(G::GrpMat) -> SeqEnum
{}

	ParabolicSubgroups(~G);
	return [* IdentifyCRG(G`Subgroups[i]`subgroup) : i in G`ParabolicSubgroups *];

end intrinsic;
