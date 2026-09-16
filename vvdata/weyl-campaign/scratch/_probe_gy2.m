// Every base this repo treats as a Guo-Yang target, and whether we have a model for it.
GY := ["10_11","10_13","10_19","10_23","111_1","134_1","146_1","14_3","14_5","15_1","15_2","15_4",
       "194_1","206_1","21_2","22_3","22_5","26_1","26_3","35_1","38_1","39_1","39_2","51_1",
       "55_1","57_1","58_1","62_1","6_11","6_17","6_19","6_29","6_31","6_37","74_1","82_1",
       "86_1","87_1","93_1","94_1"];
printf "Guo-Yang targets considered: %o\n\n", #GY;
nfull := 0; issues := [];
for g in GY do
    fn := "data/models/models_" cat g cat ".m";
    if not FileExists(fn) then Append(~issues, g cat " : NO MODEL FILE"); continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    ks := Keys(models);
    pop := [k : k in ks | #models[k] gt 0];
    k1 := [Integers()|1];
    full := IsDefined(models, k1) and #models[k1] gt 0;
    if #pop eq #ks and full then
        nfull +:= 1;
    else
        Append(~issues, Sprintf("%o : %o/%o keys populated%o", g, #pop, #ks,
                                full select "" else "  <== W=[1] EMPTY"));
    end if;
end for;
printf "fully populated: %o of %o\n\n", nfull, #GY;
printf "NOT fully populated (%o):\n", #issues;
for i in issues do printf "   %o\n", i; end for;
