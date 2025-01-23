import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;
// Following [FH] Section 2

// Start with an implmentation based on modular symbols
function CheckModularNonALInvolutionModSym(X)
    assert X`g ne 0 and X`N mod 8 eq 0;
    MDN := ModularSymbols(X`D*X`N, 2, 0);
    SDN := CuspidalSubspace(MDN);
    // We want the D-new subspace
    ps := PrimeDivisors(X`D);
    SDN_new := SDN;
    for p in ps do
        SDN_new := NewSubspace(SDN_new, p);
    end for;
    S2_MDN := ActionOnModularSymbolsBasis([2,1,0,2], MDN);
    W2_MDN := AtkinLehnerOperator(MDN, 2^Valuation(X`N, 2));
    V2_MDN := S2_MDN*W2_MDN*S2_MDN;
    SDN_new_basis := Matrix([Representation(v) : v in Basis(SDN_new)]);
    V2_SN := Solution(SDN_new_basis, SDN_new_basis * V2_MDN);
    als := [Q : Q in Divisors(X`D*X`N) | GCD(Q, X`D*X`N div Q) eq 1]; 
    ws := X`W diff {1};
    other_ws := [w : w in als | w notin ws];
    for other_w in other_ws do
        other_W_MDN := AtkinLehnerOperator(MDN, other_w);
        other_W_SN := Solution(SDN_new_basis, SDN_new_basis * other_W_MDN);
        fixed_subspace := Kernel(Matrix(V2_SN * other_W_SN) - 1);
        for w in ws do
            W_MDN := AtkinLehnerOperator(MDN, w);
            W_SN := Solution(SDN_new_basis, SDN_new_basis * W_MDN);
            al_sign := (X`D mod w eq 0) select -1 else 1;
            kerw := Kernel(Matrix(W_SN) - al_sign);
            fixed_subspace meet:= kerw;
        end for;
        d := Dimension(fixed_subspace);
        // print "d = ", d;
        assert IsEven(d);
        g := d div 2;
        // print "g = ", g;
        if (g eq 0) then
            return true, other_w;
        end if;
    end for;    
    return false;
end function;