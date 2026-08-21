// The FFT S-action must agree with the dense matrix, exactly.
AttachSpec("ShimuraQuotients.spec");
CC := ComplexField(40);
for DN in [[6,1],[10,1],[15,2]] do
    D := DN[1]; N := DN[2];
    t0 := Cputime();
    Ld := ShimuraCurveLattice(D,N);
    S, Tdiag, elts, i0 := WeilRepresentationComplex(Ld, CC : Dual := true);
    tdense := Cputime(t0);
    n := #elts;
    t1 := Cputime();
    data := VVWeilFFT(Ld, CC : Dual := true);
    tfft := Cputime(t1);
    assert data[7] eq elts and data[8] eq i0;
    assert Maximum([Abs(data[6][i] - Tdiag[i]) : i in [1..n]]) lt 1e-25;
    // compare on e_0 and on a couple of deterministic vectors
    err := 0.0;
    for trial in [1..3] do
        v := [CC | ((i*trial) mod 7) - 3 + CC.1*(((i+trial) mod 5) - 2) : i in [1..n]];
        dense := S * Matrix(CC, n, 1, v);
        fast := VVApplyS(data, v);
        err := Maximum(err, Maximum([Abs(dense[i][1] - fast[i]) : i in [1..n]]));
    end for;
    printf "%o_%o: |G| = %-6o  max|dense - fft| = %-12o  build dense %o s, fft %o s\n",
           D, N, n, ChangePrecision(err,6), tdense, tfft;
    assert err lt 1e-25;
    // and the full coset vectors must agree
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    reps := VVCosetReps(M);
    words := [VVSTWord(g) : g in reps];
    e2 := 0.0;
    t2 := Cputime();
    for w in words[1..Minimum(8,#words)] do
        a := VVRhoInvE0(S, Tdiag, w, i0);
        b := VVRhoInvE0FFT(data, w);
        e2 := Maximum(e2, Maximum([Abs(a[i][1] - b[i]) : i in [1..n]]));
    end for;
    printf "        coset rho-vectors agree to %o\n", ChangePrecision(e2,6);
    assert e2 lt 1e-22;
end for;
printf "FFT S-action validated.\n";
quit;
