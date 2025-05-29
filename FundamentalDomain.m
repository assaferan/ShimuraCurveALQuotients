import !"Geometry/GrpPSL2/Graphics/graphics.m" : AutoScale, endPsfile, firstNotSeq, labeling, showpsfile, startPsfile;

// if P is a hyperbolic polygon which is a fundamental domain for the Fuchsian group G
// inside a hyperbolic unit disc D, and z is a point in D
// Returns delta in G such that delta z is inside P
function Reduce(z, G, P)
   D := Parent(z);
   O := QuaternionOrder(G);
   deltared := ShimuraReduceUnit(G!O!1,Generators(G),D : z0 := z);
   return deltared[1];
end function;

// if P is a hyperbolic polygon which is a fundamental domain for the Fuchsian group G
// inside a hyperbolic unit disc D, and z is a point in D
// Returns whether z is inside P 
function PointInDomain(z, G, P)
   gamma := Reduce(z, G, P);
   is_in_domain := IsScalar(Quaternion(gamma));
   return is_in_domain;
end function;

function rho(P)
   return Max([Abs(x) : x in P]);
end function;

function j(g,z)
   mat_g := Matrix(g);
   return mat_g[2,1]*ComplexValue(z) + mat_g[2,2];
end function;

function Ka(n,w,G,P,k)
   g := Reduce(w, G, P);
   w_prime := g*w;
   H := UpperHalfPlane();
   z := DiscToPlane(H, w);
   w := ComplexValue(w);
   w_prime := ComplexValue(w_prime);
   ret := j(g,z)^k*(1-w)^k*w^n - (1-w_prime)^k*w_prime^n;
   return ret;
end function;

function SolveForExpansion(G, p, k, eps)
   CC<i> := Parent(p);
   D := UnitDisc(: Center := p);
   P := FundamentalDomain(G, D);
   N := 2*Ceiling(Log(eps)/Log(rho(P)));
   Q := N + 1;
   repeat
      Q +:= 1;
      ws := [D | rho(P)*Exp(2*Pi(CC)*i*m/Q) : m in [0..Q-1]];
      coeffs := [[Ka(n,w,G,P,k) : n in [0..N-1]] : w in ws];
      mat_coeffs := Matrix(coeffs);
      ker := Kernel(Transpose(mat_coeffs));
   until Dimension(ker) le 1;
   assert Dimension(ker) eq 1;
   b := Basis(ker)[1];
   return b;
end function;

procedure TestExampleVW2( : eps := 10^(-5))
   CC<i> := ComplexField();
   Omega := 0.321211772390;
   Theta := -4*Pi(CC)*Omega^2;
   a := Sqrt(-3)*(2 - Sqrt(3))*Omega^4;
   b := [a*x : x in [1/12, 0, 5/24*Theta^2, -45/(8*Factorial(4))*Theta^4, 0, 555/(4*Factorial(6))*Theta^6, 0, 57165/(8*Factorial(8))*Theta^8]];
   B := QuaternionAlgebra(6);
   O := MaximalOrder(B);
   G := FuchsianGroup(O);
   p := 1/2*(Sqrt(6) - Sqrt(2))*i;
   v := SolveForExpansion(G,p,4,eps);
   v := Vector([v[i] : i in [1..#b]]);
   assert Re(Norm(v[1]/b[1]*Vector(b) - v)) lt eps;
end procedure;

function Kc(n, r, w, g, z, w_prime, k)
   Q := #w;
   return 1/Q*&+[j(g[m],z[m])^(-k) * w_prime[m]^r*(1-w_prime[m])^k / (w[m]^n*(1-w[m])^k) : m in [1..Q]];
end function;

function J(n, m, w,g,z,k)
   return j(g[m],z[m])^(-k) / (w[m]^n*(1-w[m])^k);
end function;

function Wprime(m,r, w_prime,k)
   return w_prime[m]^r*(1-w_prime[m])^k;
end function;

function SolveForExpansion2(G,p,k,eps)
   CC<i> := Parent(p);
   H := UpperHalfPlane();
   D := UnitDisc(: Center := p);
   P := FundamentalDomain(G, D);
   N := 2*Ceiling(Log(eps)/Log(rho(P)));
   Q := N + 1;
   repeat
      Q +:= 1;
      w := [D | rho(P)*Exp(2*Pi(CC)*i*m/Q) : m in [0..Q-1]];
      g := [Reduce(ww, G, P) : ww in w];
      w_prime := [ComplexValue(g[j]*w[j]) : j in [1..#w]];
      z := [DiscToPlane(H,ww) : ww in w];
      w := [ComplexValue(ww) : ww in w];
      coeffs := [[Kc(n,r,w,g,z,w_prime,k) : r in [0..N]] : n in [0..N]];
      mat_coeffs := Matrix(coeffs);
      Jmat := Matrix([[J(n,m,w,g,z,k) : m in [1..Q]] : n in [0..N]]);
      Wprimemat := Matrix([[Wprime(m,r,w_prime,k) : r in [0..N]] : m in [1..Q]]);
      assert Re(Norm(Vector(Eltseq(Q*mat_coeffs - Jmat*Wprimemat)))) lt eps^2;
      ker := Kernel(Transpose(mat_coeffs) - 1);
   until Dimension(ker) le 1;
   assert Dimension(ker) eq 1;
   b := Basis(ker)[1];
   return b;
end function;

function AutoHeight(P)
    max := 0;
    // the following loop takes each edge e_i:=PP[i],PP[i+1] of each
    // polygon PP in the list P of polygons
    // for e_i we compute the point on e_i with largest imaginary part.
    // If the slope of the euclidean straight line between the end points
    // if e_i is greater than 1, then one of the end points has greatest
    // imaginary value.  Otherwise the greatest value is at the mid point
    // of the semi circle which extends e_i
    if #P eq 0 then
       return 0;
    end if;
    Pnontriv := [p : p in P | #p gt 0];
    H := Parent(Pnontriv[1][1]);
    if Type(H) cmpeq SpcHyd then
        return 1; // check that this is what we really want
    end if;
    inf := H!Infinity();
    for PP in Pnontriv do
       pLast:=PP[#PP];
       realLast := Real(pLast);
       imLast := Imaginary(pLast);
       for p in PP do
	  im := Imaginary(p);
	  real := Real(p);
	  // for determining width of picture:
	  if p eq inf or pLast eq inf then
	     contains_infinity := true;
	  elif real ne realLast and
	     (Abs((im-imLast)/(real - realLast)) lt 1) then
	     // find the height of the mid point on an arc
	     // between two points:
	     endpoints := ExtendGeodesic([Parent(p)|p,pLast],Parent(p));
	     // it is assumed that if (real - realLast) gt 0.01 then
	     // neither of the end points of the extended geodesic between
	     // p and pLast can be infinity.
	     im := Abs((Real(endpoints[1]) - Real(endpoints[2]))/2);
	  end if;
	  if im gt max and p ne inf then
	     max := im;
	  end if;
	  pLast := p;
	  realLast := real;
	  imLast := im;
       end for;   
    end for;
    // to get some clearance from top point:
    max := max*1.05;
    return max;    
end function;

function AutoWidth(P)
    if #P eq 0 then
       return 0;
    end if;
    Pnontriv := [p : p in P | #p gt 0];
    H := Parent(Pnontriv[1][1]);
    if Type(H) cmpeq SpcHyd then
        return -1, 1;
    end if; 
    H := UpperHalfPlaneWithCusps();
    i,j,s := firstNotSeq(P,H!Infinity());
    if i*j eq 0 then return 10;
    end if;
    max := Real(P[i][j]);
    min := Real(P[i][j]);
    for PP in P do
	for p in PP do
	    if p ne H!Infinity() then
		re := Real(p);
		if re gt max then
		    max := re;
		elif re lt min then
		    min := re;
		end if;
	    end if;
	end for;
    end for;
    return min,max;        
end function;

function DrawEdgeBetween(a,b,scale,translation,topOfScreen)
    // this function will take elements a,b of UpperHalfplaneunioncusps
    // and give the postscript for the line between them.
    // note that start and end functions for a whole path will
    // be needed to draw something.    
    MAX_ALLOWED := 10000;
    too_big := false;
    //    scale:=150;
    oldR := GetDefaultRealField();
    SetDefaultRealField(RealField(10));
    if Type(a) cmpeq SpcHydElt then
        H := Parent(a);
        assert H eq Parent(b);
        inf := H!0;
    else
        H := UpperHalfPlaneWithCusps();
        inf := H!Infinity();
    end if;
    R := RealField();
    if a eq b then
	    return "";
    end if;
    
    if a eq inf then
	    XVAL := Real(b)*scale-translation;
	    if XVAL gt MAX_ALLOWED then return ""; end if;
	    outputstring :=
	        Sprintf("%o %o lineto %o %o lineto\n",
	                XVAL,topOfScreen,
	                XVAL,Imaginary(b)*scale);
	    return outputstring;
    end if;
    if b eq inf then
	    XVAL := Real(a)*scale-translation;
	    if XVAL gt MAX_ALLOWED then return ""; end if;
	    outputstring :=
	        Sprintf("%o %o lineto %o %o lineto\n",
	                XVAL,Imaginary(a)*scale,
	                XVAL,topOfScreen);
	    return outputstring;
    end if;
    if (not (Type(a) cmpeq SpcHydElt)) and (Abs((Real(b) - Real(a))*scale) lt 0.001) then
	    XVAL := Real(a)*scale-translation;
	    if XVAL gt MAX_ALLOWED then return ""; end if;
	        outputstring :=
	            Sprintf("%o %o lineto %o %o lineto\n",
	                    XVAL,Imaginary(a)*scale,
	                    XVAL,Imaginary(b)*scale);
	    return outputstring;
    end if;
    if IsCusp(a) and IsCusp(b) then	
        A := R!Real(a)*scale;
        B := R!Real(b)*scale;
        if Real(b) gt Real(a) then
            XVAL := (A+B)/2-translation;
            if XVAL gt MAX_ALLOWED then return ""; end if;
            outputstring:=Sprintf("%o 0 %o 180 0 arcn\n",
                                    XVAL,(B-A)/2);
        else
            XVAL := (A+B)/2-translation;
            if XVAL gt MAX_ALLOWED then return ""; end if;
            outputstring:=Sprintf("%o 0 %o 180 0 arc\n",
                                    XVAL,(B-A)/2);
        end if;
        return outputstring;
    end if;
    x := Real(a);
    y := Imaginary(a);
    u := Real(b);
    v := Imaginary(b);
    if Type(a) cmpeq SpcHydElt then
        cntr, r := Geodesic(a,b);
        c := Re(cntr);
        c_im := Im(cntr);
    else
        c := (v^2 + u^2 - x^2 - y^2)/(u-x)/2;
        r := SquareRoot((x-c)^2 + y^2);
        c_im := 0;
    end if;
    dist1 := x - c;
    dist2 := u - c;  
    PI:=Pi(R);
    cos1 := dist1/r;  cos2 := dist2/r;
    // make sure these values are in the domain of Arccos
    if Abs(cos1) gt 1 then cos1 := Round(cos1); end if;
    if Abs(cos2) gt 1 then cos2 := Round(cos2); end if;
    ang1 := Real((Arccos(cos1))*180/PI); // mod 180;
    ang2 := Real((Arccos(cos2))*180/PI); // mod 180;
    if (c_im gt y) then
        ang1 := 360 - ang1;
    end if;
    if (c_im gt v) then
        ang2 := 360 - ang2;
    end if;   
    if (Type(a) cmpeq SpcHydElt) or (ang1 gt ang2) then
	    XVAL := c*scale-translation;
        YVAL := c_im*scale-translation;
	    if XVAL gt MAX_ALLOWED then return ""; end if;
	    outputstring:=Sprintf("%o %o %o %o %o	
	                            arcn\n",XVAL, YVAL, r*scale,ang1,ang2);
    else
	    XVAL := c*scale-translation;
        YVAL := c_im*scale-translation;
	    if XVAL gt MAX_ALLOWED then return ""; end if;
	    outputstring:=Sprintf("%o %o %o %o %o
	                            arc\n",XVAL, YVAL, r*scale,ang1,ang2);
    end if;
    SetDefaultRealField(oldR);
    return outputstring;
end function;

intrinsic IsCusp(x::SpcHydElt) -> BoolElt
{Returns true is and only if the element x of a unit disc is a cusp.}
    D := Parent(x);
    infty := D!0;
    on_boundary := Abs(Abs(x - infty) - 1) lt D`eps;
    is_infty := Abs(x - infty) lt D`eps;
    return on_boundary or is_infty;
end intrinsic;

function DrawPolygon(S,colour,pencolour,
   scale,outline,fill,translation,topOfScreen,
   Radius)
    // S is a sequence of elements of Upperhalfplaneunioncusps or a hyperbolic unit disc
    // colour is a triple of numbers between 0 and 1.
    // try scale around 100 - it's the number of pixels per unit
    // try top of screen about 400 - or as large as you want.
    
    //   topOfScreen := 400;
    //    scale:=150;
    MAX_ALLOWED := 10000;
    if Type(S) ne SeqEnum or  #S eq 0 then
	    return "";
    end if;
    c1,c2,c3 := Explode(colour);
    pc1,pc2,pc3 := Explode(pencolour);
    //    printf "%o %o\n",Parent(S),Parent([S[1]]);
    elt_of_S := func<i | i le #S select S[i] else S[1]>;
    S1 := [Parent(S[1])|elt_of_S(i) : i in [1..(#S+1)]];   
    //       S1 := S cat [S[1]];
    if #S eq 1 and (not (IsCusp(S[1]) and
       (Type(ExactValue(S[1])) eq SetCspElt and ExactValue(S[1])
       eq Cusps()!Infinity()))) then
       x := scale*Real(S[1]) - translation;
       y := scale*Imaginary(S[1]);
       startString := Sprintf("%o %o %o setrgbcolor\n",c1,c2,c3);
       polystring :=
       Sprintf("newpath  \n %o %o %o 0 360 arc fill\n\n",x,y,Radius);
       return startString cat polystring;
    end if;
       
    
    startStringA := Sprintf("%o %o %o setrgbcolor\n",c1,c2,c3);
    startStringB := Sprintf("%o %o %o setrgbcolor\n",pc1,pc2,pc3);
    if not
       (IsCusp(S[1]) and Type(ExactValue(S[1])) eq SetCspElt and
       ExactValue(S[1]) eq Cusps()!Infinity())
       then
       x := scale*Real(S[1]) - translation;
       y := scale*Imaginary(S[1]);
    else
       y := topOfScreen;
       x:=0 - translation;
    end if;    
    polystring := Sprintf("newpath\n %o %o moveto\n",x,y);
    if x gt MAX_ALLOWED then return ""; end if;
    for i in [2..#S1] do
	    newedge := DrawEdgeBetween(S1[i-1],S1[i],scale,translation,topOfScreen);
        // give up on the polygon if it's not in drawable area.
        // this gives a bug ! need to fix!
        //	if newedge eq "" then return ""; end if;
        polystring := polystring cat newedge;
    end for;
    endstringA := "fill\n\n";
    endstringB := "stroke\n\n";
    if outline and fill then
       outputstring := startStringA cat polystring cat endstringA
       cat startStringB cat polystring cat endstringB;
    elif fill then
       outputstring := startStringA cat polystring cat endstringA;
    else
       outputstring := startStringB cat polystring cat endstringB;
    end if;
    return outputstring;
end function;

intrinsic DisplayPolygonsUnitDisc(P::SeqEnum,filename::MonStgElt:
    Colours := [1,1,0],
    PenColours := [0,0,0],
    Outline := [true : p in P],
    Fill := [true : p in P],
    Fontsize := 2,
    Labels := [Rationals()|0,1],
    Radius := 0.5,
    Pixels := 300,
    Size := [],
    Show := false,
    Overwrite := true) -> SeqEnum
    {Create/display a postscript drawing of the polygon with given vertices}

    // P is a sequence of sequences, each representing a polygon
    // colour is a list of triples of numbers between 0 and 1.
    // filename is some string where the file should be written
    // Show is a boolean; if it's true pop up the ps file with
    // a system command.
    //
    require Type(Labels) cmpeq SeqEnum:
    "Labels must be a sequence of rationals or cusps";
    require (#Labels eq 0 or Parent(Labels[1]) cmpeq Rationals()
    or Parent(Labels[1]) cmpeq Integers()
    or Parent(Labels[1]) cmpeq Cusps()
    or Type(Labels[1]) cmpeq SpcHypElt
    or Type(Labels[1]) cmpeq SpcHydElt):
    "Labels must be a sequence of rationals or cusps";    
    if  #Labels gt 0 and Parent(Labels[1]) cmpeq Integers() then
       Labels := [Rationals()|x : x in Labels];
    end if;

    // for hyperbolic space elements, only add labels if they are
    // cusps
    if  #Labels gt 0 and Type(Labels[1]) cmpeq SpcHypElt then
       Labels := [Cusps()|ExactValue(x) : x in Labels | IsCusp(x)];
    end if;

    require Parent(Pixels) cmpeq Integers(): "Pixels should be an integer";
    if Pixels lt 10 then
       Pixels := 10;
    end if;
    
    if not ((Type(Size) eq SeqEnum and #Size eq 4 and Type(Universe(Size)) in {RngInt, FldRe, FldRat}))     then
       Autoscale := true;
    else
       if (Type(Universe(Size)) ne FldRe) then
	  Size := [RealField() | x : x in Size];
       end if;
       Autoscale := false;
    end if;
       
    require #P gt 0: "Please give a list of polygons";    
    require Type(P[1]) in {SeqEnum, SpcHypElt, SpcHydElt, SetCspElt}:
    "Please give a list of polygons";
    require Type(P[1]) in {SpcHypElt, SpcHydElt, SetCspElt}
    or Type(P[1][1]) in {SpcHypElt, SpcHydElt, SetCspElt}:
    "Please give a list of polygons";
    if Type(P[1]) in {SpcHypElt, SpcHydElt, SetCspElt} then
       P := [P];
    end if;    
    if Type(P[1][1]) eq SetCspElt then
       H := UpperHalfPlaneWithCusps();
       P := [[H|x : x in p] : p in P];
    end if;
    
    // note requirements for Colours must come after requirements for P
    require Type(Colours) eq SeqEnum and #Colours gt 0:
    "Colours must be a sequence of sequences of numbers";
    require Type(Colours[1]) in {SeqEnum,RngIntElt,FldRatElt,FldReElt}:
    "Colours must be a sequence of sequences of numbers";
    require (Type(Colours[1]) in {RngIntElt,FldRatElt,FldReElt}) or
    (Type(Colours[1][1]) in {RngIntElt,FldRatElt,FldReElt}):
    "Colours must be a sequence of sequences of 3 numbers";
    require ((Type(Colours[1]) in {RngIntElt,FldRatElt,FldReElt})
    and #Colours ge 3)
    or (Min([#Colours[i] : i in [1..#Colours]]) ge 3):
    "Colours must be a sequence of sequences of 3 integers";
    if Type(Colours[1]) in {RngIntElt,FldRatElt,FldReElt} then
       Colours := [[Colours[i] : i in [1..3]]];
    end if;
    if #Colours lt #P then       
       Colours := Colours cat [Colours[1] : i in [(#Colours+1)..#P]];
    end if;
    if Type(Colours[1][1]) eq FldRatElt then
       Colours := [[RealField() | 1.0*i : i in c] : c in Colours];
    end if;


    // note requirements for PenColours must come after requirements for P
    require Type(PenColours) eq SeqEnum and #PenColours gt 0:
    "PenColours must be a sequence of sequences of numbers";
    require Type(PenColours[1]) in {SeqEnum,RngIntElt,FldRatElt,FldReElt}:
    "PenColours must be a sequence of sequences of numbers";
    require (Type(PenColours[1]) in {RngIntElt,FldRatElt,FldReElt}) or
    (Type(PenColours[1][1]) in {RngIntElt,FldRatElt,FldReElt}):
    "PenColours must be a sequence of sequences of 3 numbers";
    require ((Type(PenColours[1]) in {RngIntElt,FldRatElt,FldReElt})
    and #PenColours ge 3)
    or (Min([#PenColours[i] : i in [1..#PenColours]]) ge 3):
    "PenColours must be a sequence of sequences of 3 integers";
    if Type(PenColours[1]) in {RngIntElt,FldRatElt,FldReElt} then
       PenColours := [[PenColours[i] : i in [1..3]]];
    end if;
    if #PenColours lt #P then       
       PenColours := PenColours cat [PenColours[1] : i in [(#PenColours+1)..#P]];
    end if;
    if Type(PenColours[1][1]) eq FldRatElt then
       PenColours := [[RealField() | 1.0*i : i in c] : c in PenColours];
    end if;

    
       
    require (Type(Outline) eq BoolElt)
    or (Type(Outline) eq SeqEnum and #Outline gt 0):
    "Outline must be a boolean or sequence of booleans";
    require (Type(Outline) eq BoolElt) or
    (#Outline gt 0 and Type(Outline[1]) eq BoolElt):
    "Outline must be a boolean or sequence of booleans";
    if Type(Outline) eq BoolElt then
       Outline := [Outline];
    end if;
    if #Outline lt #P then
       Outline := Outline cat [Outline[1] : i in [(#Outline+1)..#P]];
    end if;

    
    require Type(Fill) in {SeqEnum,BoolElt} and
    (Type(Fill) eq BoolElt or #Fill gt 0):
    "Fill must be a boolean or sequence of booleans";
    require (Type(Fill) eq BoolElt) or
    (#Fill gt 0 and Type(Fill[1]) eq BoolElt):
    "Fill must be a boolean or sequence of booleans";
    if Type(Fill) eq BoolElt then
       Fill := [Fill];
    end if;
    if #Fill lt #P then       
       Fill := Fill cat [Fill[1] : i in [(#Fill+1)..#P]];
    end if;

    require Parent(Fontsize) cmpeq Integers(): "Fontsize should be an integer";
    require Type(Show) cmpeq BoolElt: "Show must be a boolean";
    require Type(Autoscale) cmpeq BoolElt: "Autoscale must be a boolean";

    TOO_BIG := false;
    // note scale in ps file, currently 4 4 scale, above, thus the following:
    scale4 := 4;
    if Type(P[1][1]) cmpeq SpcHydElt then
        scale4 := 8;
    end if;    

    if Autoscale then
       h := AutoHeight(P);    
       w1,w2 := AutoWidth(P);
       S := AutoScale(Pixels,h,w1,w2);           
       // following needs improving!
       // the points is that currently the height is computed
       // so that all vertices of a polygon are included in a
       // displayed area, but if the vertices are all cusps in
       // the real line, then they all have 0 hight, so one
       // should compute the maximum points on the lines between
       // the instead, but so far, this is not done, and we
       // just take a hight 0.5.
       if h eq 0 then h:=0.5; end if; 
    else
       w1 := Size[1];
       w2 := Size[2];
       h := Size[3];
       S := Size[4];
    end if;

    translation := Integers()!Ceiling(S*w1/scale4);
    w := w2 - w1;
    if Abs(w*S) lt 0.001
       then width:=0.5;
    elif Abs(w*S) gt 10000 then
       width:=10000;      
    else
       width := Min(Integers()!Ceiling(w*S),10000);
    end if;
    height := Min(Integers()!Ceiling(h*S),10000);
    topOfScreen := height;

    writable := startPsfile(filename,translation,width,height,
                            Overwrite,Fontsize);
    if writable then

        oldR := GetDefaultRealField();
        SetDefaultRealField(RealField(10));
        for i in [1..#P] do
	        polystring := DrawPolygon(P[i],Colours[i],PenColours[i],1.0*(S/scale4),
	                                    Outline[i],Fill[i],translation,topOfScreen,Radius);
	        if #P[i] gt 0 and polystring eq "" then
	            TOO_BIG := true;
	        end if;
	        Write(filename,polystring);
        end for;
        if #Labels gt 0 then
	        Write(filename,labeling(Labels,1.0*(S/4),Fontsize,translation));
        end if;
        endPsfile(filename);
        if Show then
	        showpsfile(filename);
        end if;
        if TOO_BIG then
	        printf "Some polygons have not been drawn, due to size
	                limitations";
        end if;
        SetDefaultRealField(oldR);
    else
        printf "No file created.\n";
    end if;
    return [w1,w2,h,S];
end intrinsic;