declare type SchoferTable;

declare attributes SchoferTable: Values, Keys_fs, sIndex, sTildeIndex, Discs, Xstar, BorcherdsForms, RowScales, FldsOfDefn, K_idxs, Curves;

intrinsic Print(x::SchoferTable, L::MonStgElt)
    {Print X at level L}
    if L eq "Default" then
        print x`Values;
    elif L eq "Magma" then
        printf "CreateSchoferTable(%m, %m, %m)", x`Values, x`Keys_Ffs, x`Discs;
    end if;
end intrinsic;

intrinsic CreateSchoferTable(vals::SeqEnum, keys_fs::SeqEnum, ds::SeqEnum, curves::SeqEnum[ShimuraQuot], Xstar::ShimuraQuot) -> SchoferTable
{}
    x := New(SchoferTable);
    x`Values := vals;
    x`Keys_fs := keys_fs;
    x`Discs := ds;
    x`sIndex := Index(keys_fs, -1);
    x`sTildeIndex := Index(keys_fs, -2);
    x`K_idxs := [i : i->k in keys_fs | k gt 0];
    x`Curves := curves;
    x`Xstar := Xstar;
    curveid := Xstar`CurveID;

    flds := AssociativeArray();
    for k in x`K_idxs do
        flds[keys_fs[k]] :=  AssociativeArray();
        for d in x`Discs do
            flds[keys_fs[k]][d] :=  FieldsOfDefinitionOfCMPoint(curves[keys_fs[k]], d);
        end for;
    end for;
    flds[curveid] := AssociativeArray();
    for d in x`Discs do
        flds[curveid][d] :=  FieldsOfDefinitionOfCMPoint(curves[curveid], d);
    end for;
    x`FldsOfDefn := flds;

    return x;
end intrinsic;

intrinsic CreateSchoferTable(vals::List, keys_fs::SeqEnum, ds::SeqEnum, curves::SeqEnum[ShimuraQuot], Xstar ::ShimuraQuot) -> SchoferTable
{}
    x := New(SchoferTable);
    x`Values := vals;
    x`Keys_fs := keys_fs;
    x`Discs := ds;
    x`sIndex := Index(keys_fs, -1);
    x`sTildeIndex := Index(keys_fs, -2);
    x`K_idxs := [i : i->k in keys_fs | k gt 0];
    x`Curves := curves;
    x`Xstar := Xstar;
    curveid := Xstar`CurveID;

    flds := AssociativeArray();
    for k in x`K_idxs do
        flds[keys_fs[k]] :=  AssociativeArray();
        for d in x`Discs do
            flds[keys_fs[k]][d] :=  FieldsOfDefinitionOfCMPoint(curves[keys_fs[k]], d);
        end for;
    end for;
    flds[curveid] := AssociativeArray();
    for d in x`Discs do
        flds[curveid][d] :=  FieldsOfDefinitionOfCMPoint(curves[curveid], d);
    end for;

    x`FldsOfDefn := flds;
    return x;
end intrinsic;

intrinsic UpdateFieldsOfDefn(x::SchoferTable, d::RngIntElt)
    {Update fields of definition to include d}
    curves := x`Curves;
    flds := x`FldsOfDefn;
    k_idxs := x`K_idxs;
    keys_fs := x`Keys_fs;
    curveid := (x`Xstar)`CurveID;
    for k in k_idxs do
        flds[keys_fs[k]][d] :=  FieldsOfDefinitionOfCMPoint(curves[keys_fs[k]], d);
    end for;
    flds[curveid][d] :=  FieldsOfDefinitionOfCMPoint(curves[curveid], d);
    x`FldsOfDefn := flds;
end intrinsic;
