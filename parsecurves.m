data := Read("nonmodularinvolutiondata.out");
data := Split(data, "\n");

for curve in data do
    curvedata := Split(curve,":");

    if (eval curvedata[2]) eq 1 then
        index := eval curvedata[1];
        curves[index]`IsHyp := true;
        curves[index]`IsSubhyp := true;
        inv := curvedata[3];
        curves[index]`TestInWhichProved := "ModularNonALInvolution " cat inv;
    elif (eval curvedata[2]) eq 0 then
        inv := curvedata[3];
        num := curvedata[4];
        index := eval curvedata[1];
        curves[index]`IsHyp := false;
        curves[index]`IsSubhyp := false;
        curves[index]`TestInWhichProved := "ModularNonALInvolution " cat inv cat " " cat num;
    end if;
end for;
