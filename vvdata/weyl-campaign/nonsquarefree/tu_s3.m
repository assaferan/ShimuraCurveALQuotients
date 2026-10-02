// Tu's t_4 on X_0^15(4)/<w_3,w_5> (Pacific J. Math. 269 (2014), Lemma 13): values
//   +-1/sqrt(-3) at discriminant -12, +-sqrt(-15)/5 at -15, (+-1 +- sqrt(-15))/8 at -60.
// The S_3 Galois group of this curve over X_0^15(1)^* acts on the t_4-line by Mobius maps; the two
// -12 values are the fixed points of the 3-cycle (the (3,3) fibre over the -3 star point), and the
// six values at -15 and -60 should be one orbit (the unramified fibre over the -15 star point).
// Test: is there an order-3 Mobius map rho fixing +-1/sqrt(-3) that permutes the six values?
K<s3> := QuadraticField(-3); L<s15> := QuadraticField(-15);
F := Compositum(K, L); a := F!s3; b := F!s15;
f1 := 1/a; f2 := -1/a;                                  // fixed points of rho
six := [b/5, -b/5, (1+b)/8, (1-b)/8, (-1+b)/8, (-1-b)/8];
// an order-3 Mobius map with fixed points f1, f2:  (rho(t) - f1)/(rho(t) - f2) = w (t - f1)/(t - f2), w^3 = 1
w := (-1 + a)/2;                                        // a primitive cube root of unity
rho := func<t | (f1 - w*f2*(t - f1)/(t - f2)) / (1 - w*(t - f1)/(t - f2))>;
img := [rho(t) : t in six];
printf "rho permutes the six values: %o\n", Set(img) eq Set(six);
printf "orbits of rho: %o\n", [[Position(six, t), Position(six, rho(t)), Position(six, rho(rho(t)))] : t in six[1..2]];
// and an involution swapping f1, f2 and normalising rho: sigma(t) = -t swaps +-1/sqrt(-3)
sig := func<t | -t>;
printf "t -> -t permutes the six values: %o\n", Set([sig(t) : t in six]) eq Set(six);
printf "sigma rho sigma = rho^-1: %o\n", forall{t : t in six | sig(rho(sig(t))) eq rho(rho(t))};
