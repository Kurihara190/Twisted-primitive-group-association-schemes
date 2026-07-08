# GAP verification for the section:
#   "A twisted fission of the group association scheme of A8"
#
# This script avoids constructing a 20160 by 20160 relation matrix.  It checks
# the assertions used in the section directly from permutations:
#   class sizes, S-identities, the affine subgroup H, H-orbits on Sylow-7
#   cells, the construction of D_1 and D_2, T-identities, local graph
#   connectedness, and the triangle-count witness for non-Schurity.

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
if LoadPackage("grape") = fail then
    Error("The GAP package GRAPE is required for local graph calculations.");
fi;
LogTo(Concatenation(OutDir, "/result_A8.txt"));

Print("\n");
Print("============================================================\n");
Print("Twisted fission checks for A8\n");
Print("============================================================\n\n");

G := AlternatingGroup(8);
elts := Elements(G);
e := One(G);

posG := NewDictionary(e, true);
for i in [1..Length(elts)] do
    AddDictionary(posG, elts[i], i);
od;

ClassOf := x -> AsList(ConjugacyClass(G, x));

# Conjugacy classes in the order used in the text.
C0 := [e];
C3 := ClassOf((1,2,3));
C22 := ClassOf((1,2)(3,4));
C2222 := ClassOf((1,2)(3,4)(5,6)(7,8));
C3221 := ClassOf((1,2,3)(4,5)(6,7));
C3311 := ClassOf((1,2,3)(4,5,6));
C44 := ClassOf((1,2,3,4)(5,6,7,8));
C4211 := ClassOf((1,2,3,4)(5,6));
C5 := ClassOf((1,2,3,4,5));
C531 := ClassOf((1,2,3,4,5)(6,7,8));
C532 := ClassOf((1,3,5,2,4)(6,7,8));
C62 := ClassOf((1,2,3,4,5,6)(7,8));
x71 := (1,2,3,4,5,6,7);
C711 := ClassOf(x71);
C712 := ClassOf(x71^3);

classNames := ["0","(3)","(2,2)","(2,2,2,2)","(3,2,2,1)",
               "(3,3,1,1)","(4,4)","(4,2,1,1)","(5)",
               "(5,3)_1","(5,3)_2","(6,2)","(7,1)_1","(7,1)_2"];
classes := [C0,C3,C22,C2222,C3221,C3311,C44,C4211,C5,C531,C532,C62,C711,C712];

ValueSetsOnParts := function(vec, parts)
    return List(parts, C -> Set(List(C, g -> vec[LookupDictionary(posG, g)])));
end;

Print("1. Conjugacy class sizes\n");
Print("   Class sizes in order ", classNames, " = ",
      List(classes, Length), "\n\n");
Print("   Label representative (1,2,3,4,5,6,7) in C_(7,1)_1? ",
      x71 in C711, "\n");
Print("   Label representative (1,2,3,4,5)(6,7,8) in C_(5,3)_1? ",
      (1,2,3,4,5)(6,7,8) in C531, "\n\n");

U := Concatenation(C711, C712);
signS := List(elts, g -> 0);
for g in C711 do signS[LookupDictionary(posG, g)] := 1; od;
for g in C712 do signS[LookupDictionary(posG, g)] := -1; od;

Print("2. Original S identities for split (7,1) classes\n");
fixedClassData := [
    [C3,"C_(3)"], [C3221,"C_(3,2,2,1)"], [C3311,"C_(3,3,1,1)"],
    [C5,"C_(5)"], [C531,"C_(5,3)_1"], [C532,"C_(5,3)_2"],
    [C62,"C_(6,2)"], [C22,"C_(2,2)"], [C2222,"C_(2,2,2,2)"],
    [C44,"C_(4,4)"], [C4211,"C_(4,2,1,1)"]
];
eigData := [];
sLinearDefects := [];

for item in fixedClassData do
    coeff := List(elts, g -> 0);
    for a in item[1] do
        for b in U do
            p := LookupDictionary(posG, a*b);
            coeff[p] := coeff[p] + signS[LookupDictionary(posG, b)];
        od;
    od;
    lambda := coeff[LookupDictionary(posG, C711[1])];
    defects := 0;
    for idx in [1..Length(elts)] do
        if coeff[idx] <> lambda * signS[idx] then
            defects := defects + 1;
        fi;
    od;
    Add(eigData, [item[1], lambda, item[2]]);
    Add(sLinearDefects, [item[2], lambda, defects]);
    Print("   ", item[2], " * S value sets on classes ",
          classNames, " = ", ValueSetsOnParts(coeff, classes), "\n");
od;
Print("   Computed eigenvalue data [class, lambda, defects] = ",
      sLinearDefects, "\n\n");

coeff := List(elts, g -> 0);
for a in U do
    ia := LookupDictionary(posG, a);
    for b in U do
        coeff[LookupDictionary(posG, a*b)] :=
            coeff[LookupDictionary(posG, a*b)] + signS[ia] * signS[LookupDictionary(posG, b)];
    od;
od;
Print("   S^2 value sets on classes ", classNames, " = ",
      ValueSetsOnParts(coeff, classes), "\n\n");
sSquareCoeff := ShallowCopy(coeff);

# Build H from the affine formula on F_2^2 x F_2.
Print("3. The affine subgroup H\n");
labelToVec := [[1,0,0],[0,1,0],[1,1,0],[0,0,1],
               [1,0,1],[0,1,1],[1,1,1],[0,0,0]];
gl22 := List(Elements(GL(2,2)),
    A -> List([1,2], i -> List([1,2], j -> Int(A[i][j]))));

Hgens := [];
for A in gl22 do
    for b1 in [0,1] do for b2 in [0,1] do
    for c1 in [0,1] do for c2 in [0,1] do
    for eta in [0,1] do
        img := [];
        for lab in [1..8] do
            u1 := labelToVec[lab][1]; u2 := labelToVec[lab][2]; z := labelToVec[lab][3];
            v := [ (A[1][1]*u1 + A[1][2]*u2 + b1 + z*c1) mod 2,
                   (A[2][1]*u1 + A[2][2]*u2 + b2 + z*c2) mod 2,
                   (z + eta) mod 2 ];
            img[lab] := Position(labelToVec, v);
        od;
        Add(Hgens, PermList(img));
    od; od; od; od; od;
od;
H := Group(Set(Hgens));
eltsH := Elements(H);
Print("   |H| = ", Size(H), "\n");
Print("   H <= A8? ", IsSubgroup(G, H), "\n");
Print("   StructureDescription(H) = ", StructureDescription(H), "\n");

reps := [
    (1,2,7,3,5,6,4), (1,2,7,6,3,5,4), (1,2,7,5,6,3,4),
    (1,2,3,7,6,5,4), (1,2,3,7,5,6,4), (1,2,6,3,5,7,4),
    (1,2,5,3,7,6,4)
];
cellOrbitSizes := [];
seenCells := [];
for r in reps do
    orbit := [];
    for h in eltsH do
        AddSet(orbit, Set(Elements(Group(r^h))));
    od;
    Add(cellOrbitSizes, Length(orbit));
    UniteSet(seenCells, orbit);
od;
Print("   H-orbit sizes of the seven listed cells = ", cellOrbitSizes, "\n");
Print("   They cover all 960 Sylow-7 cells? ",
      Sum(cellOrbitSizes)=960 and Length(seenCells)=960, "\n\n");

Print("4. Construct D_1 and D_2\n");
Q := [1,2,4];
Nres := [3,5,6];
D1 := [];
D2 := [];
for r in reps do
    for h in eltsH do
        for q in Q do AddSet(D1, (r^q)^h); od;
        for q in Nres do AddSet(D2, (r^q)^h); od;
    od;
od;
twClasses := Concatenation(classes{[1..12]}, [D1, D2]);
Print("   |D_1| = ", Length(D1), ", |D_2| = ", Length(D2), "\n");
Print("   D_1 and D_2 partition all (7,1)-elements? ",
      Intersection(D1,D2)=[] and Set(Concatenation(D1,D2))=Set(U), "\n");
Print("   D_1^-1 = D_2 ? ", Set(List(D1, x -> x^-1)) = Set(D2), "\n\n");

Print("5. Twisted T identities\n");
signT := List(elts, g -> 0);
for g in D1 do signT[LookupDictionary(posG, g)] := 1; od;
for g in D2 do signT[LookupDictionary(posG, g)] := -1; od;

for item in eigData do
    coeff := List(elts, g -> 0);
    for a in item[1] do
        for b in U do
            p := LookupDictionary(posG, a*b);
            coeff[p] := coeff[p] + signT[LookupDictionary(posG, b)];
        od;
    od;
    defects := 0;
    for idx in [1..Length(elts)] do
        if coeff[idx] <> item[2] * signT[idx] then
            defects := defects + 1;
        fi;
    od;
    Print("   ", item[3], " * T value sets on [D_1,D_2] = ",
          ValueSetsOnParts(coeff, [D1,D2]), "\n");
    Print("      defects against computed lambda ", item[2],
          " = ", defects, "\n");
    Print("      values on unchanged classes ",
          classNames{[1..12]}, " = ",
          ValueSetsOnParts(coeff, classes{[1..12]}), "\n");
od;

coeff := List(elts, g -> 0);
for a in U do
    ia := LookupDictionary(posG, a);
    for b in U do
        coeff[LookupDictionary(posG, a*b)] :=
            coeff[LookupDictionary(posG, a*b)] + signT[ia] * signT[LookupDictionary(posG, b)];
    od;
od;
Print("   T^2 value sets on classes ", classNames, " = ",
      ValueSetsOnParts(coeff, classes), "\n");
squareDefects := 0;
for idx in [1..Length(elts)] do
    if coeff[idx] <> sSquareCoeff[idx] then
        squareDefects := squareDefects + 1;
    fi;
od;
Print("   T^2 coefficient defects against computed S^2 = ",
      squareDefects, "\n\n");

Print("6. Local graph connectedness\n");

for pair in [[C711,C3,"Gamma_C(3)(C_(7,1)_1)"],
             [C712,C3,"Gamma_C(3)(C_(7,1)_2)"],
             [D1,C3,"Gamma_C(3)(D_1)"],
             [D2,C3,"Gamma_C(3)(D_2)"]] do
    graphLocal := Graph(Group(()), pair[1],
        function(x, g) return x; end,
        function(x, y) return x <> y and y*x^-1 in pair[2]; end,
        true);
    comps := List(ConnectedComponents(graphLocal), Length);
    Sort(comps);
    Print("   ", pair[3], ": degree set = ",
          VertexDegrees(graphLocal), ", component sizes = ", comps, "\n");
od;

Print("\n7. Triangle counts in Gamma_C(3)(D_1)\n");

graphD1C3 := Graph(Group(()), D1,
    function(x, g) return x; end,
    function(x, y) return x <> y and y*x^-1 in C3; end,
    true);
adj := graphD1C3.adjacencies;
neighborSets := List(adj, Set);

triangleCounts := [];
for i in [1..graphD1C3.order] do
    tau := 0;
    for jj in [1..Length(adj[i])] do
        for kk in [jj+1..Length(adj[i])] do
            j := adj[i][jj];
            k := adj[i][kk];
            if k in neighborSets[j] then
                tau := tau + 1;
            fi;
        od;
    od;
    Add(triangleCounts, tau);
od;

Print("   Set of triangle counts per vertex = ", Set(triangleCounts), "\n");
Print("   The graph is not vertex-transitive? ",
      Set(triangleCounts) <> [triangleCounts[1]], "\n");
Print("   Triangle-count distribution = ",
      Collected(triangleCounts), "\n");

PrintTo(Concatenation(OutDir, "/A8_twisted_Schur_partition.txt"),
        "A8TwistedSchurPartition := ", twClasses, ";\n");

Print("\nDone.\n");

LogTo();
