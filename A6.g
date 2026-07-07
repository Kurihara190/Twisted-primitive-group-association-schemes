# Verification script for the A6 group association scheme and its twist.
#
# This file follows the section
#   "The group association scheme X(A6) and an example of a twist"
# in the order in which the assertions appear.

# Write printed verification output to out/result_A6.txt relative to the
# directory from which this script is run.
OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
if LoadPackage("grape") = fail then
    Error("The GAP package GRAPE is required for local graph calculations.");
fi;
LogTo(Concatenation(OutDir, "/result_A6.txt"));

Print("\n");
Print("============================================================\n");
Print("A6 group association scheme and the twisted partition\n");
Print("============================================================\n\n");

G := AlternatingGroup(6);
elts := Elements(G);
e := One(G);
x0 := (1,2,3,4,5);

posG := NewDictionary(e, true);
for i in [1..Length(elts)] do
    AddDictionary(posG, elts[i], i);
od;

ValueSetsOnParts := function(vec, parts)
    return List(parts, C -> Set(List(C, g -> vec[LookupDictionary(posG, g)])));
end;

Print("1. Basic group data\n");
Print("   |A6| = ", Size(G), "\n");
Print("   x0 = ", x0, "\n\n");

# ---------------------------------------------------------------------------
# Conjugacy classes in the order used in the text:
#   0, (2,2), (3), (3,3), (4,2), (5)_1, (5)_2.
# The class (5)_1 is chosen to be the A6-conjugacy class of x0=(1,2,3,4,5).
# ---------------------------------------------------------------------------

C0 := [e];
C22 := [];
C3 := [];
C33 := [];
C42 := [];
All5 := [];

for g in elts do
    if Order(g) = 2 then
        Add(C22, g);
    elif Order(g) = 3 then
        moved := 0;
        for i in [1..6] do
            if i^g <> i then
                moved := moved + 1;
            fi;
        od;
        if moved = 3 then
            Add(C3, g);
        else
            Add(C33, g);
        fi;
    elif Order(g) = 4 then
        Add(C42, g);
    elif Order(g) = 5 then
        Add(All5, g);
    fi;
od;

C51 := AsList(ConjugacyClass(G, x0));
C52 := Difference(All5, C51);

classes := [C0, C22, C3, C33, C42, C51, C52];
classNames := ["0", "(2,2)", "(3)", "(3,3)", "(4,2)", "(5)_1", "(5)_2"];

Print("2. Conjugacy class sizes in the paper's order\n");
for i in [1..Length(classes)] do
    Print("   C_", classNames[i], " has size ", Length(classes[i]), "\n");
od;
Print("   Total size = ", Sum(List(classes, Length)), "\n\n");

# ---------------------------------------------------------------------------
# Intersection matrices B_i=(p_{ij}^k)_{j,k}.
# Here p_{ij}^k is the coefficient of an element of C_k in
#   underline(C_i) * underline(C_j).
# For every pair (i,j), we count all products a*b and distribute them among
# the seven basic sets.
# ---------------------------------------------------------------------------

Print("3. Intersection matrices B_i\n");

B := [];
for i in [1..7] do
    B[i] := [];
    for j in [1..7] do
        row := List([1..7], t -> 0);
        for a in classes[i] do
            for b in classes[j] do
                prod := a*b;
                for k in [1..7] do
                    if prod in classes[k] then
                        row[k] := row[k] + 1;
                        break;
                    fi;
                od;
            od;
        od;
        # Divide by |C_k| to get the coefficient for a fixed element of C_k.
        for k in [1..7] do
            row[k] := row[k] / Length(classes[k]);
        od;
        B[i][j] := row;
    od;
od;


for i in [1..7] do
    Print("   B_", classNames[i], " = ",
          B[i], "\n");
od;
Print("\n");


# ---------------------------------------------------------------------------
# The original S = underline(C_(5)_1) - underline(C_(5)_2) identities.
# We verify the coefficients on every group element.
# ---------------------------------------------------------------------------

Print("4. Original S identities\n");

signS := List(elts, g -> 0);
for g in C51 do
    signS[Position(elts, g)] := 1;
od;
for g in C52 do
    signS[Position(elts, g)] := -1;
od;

sActionOnFixedClasses := [[C0, 1]];

for pair in [[C22,"C_(2,2)"], [C42,"C_(4,2)"],
             [C3,"C_(3)"], [C33,"C_(3,3)"]] do
    coeff := List(elts, g -> 0);
    for a in pair[1] do
        for b in All5 do
            coeff[Position(elts, a*b)] :=
                coeff[Position(elts, a*b)] + signS[Position(elts, b)];
        od;
    od;
    Add(sActionOnFixedClasses, [pair[1], coeff[LookupDictionary(posG, C51[1])]]);
    Print("   ", pair[2], " * S value sets on classes ",
          classNames, " = ", ValueSetsOnParts(coeff, classes), "\n");
od;

coeff := List(elts, g -> 0);
for term in [[C51, C51, 1], [C52, C52, 1],
             [C51, C52, -1], [C52, C51, -1]] do
    for a in term[1] do
        for b in term[2] do
            coeff[Position(elts, a*b)] :=
                coeff[Position(elts, a*b)] + term[3];
        od;
    od;
od;
Print("   S^2 value sets on classes ", classNames, " = ",
      ValueSetsOnParts(coeff, classes), "\n\n");
sSquareCoeff := ShallowCopy(coeff);

# ---------------------------------------------------------------------------
# Define cr(x), split each 5-cycle class by cr=2 or cr=4, and define
#   D_1 = C_(5)_1^(2) union C_(5)_2^(4),
#   D_2 = C_(5)_1^(4) union C_(5)_2^(2).
# The loops below compute cr(x) directly from the undirected 5-cycle.
# ---------------------------------------------------------------------------

Print("5. The cut number cr(x) and the sets D_1, D_2\n");

A := [1,2,3];
BB := [4,5,6];
C51cr2 := [];
C51cr4 := [];
C52cr2 := [];
C52cr4 := [];

for x in All5 do
    movedPoints := [];
    for i in [1..6] do
        if i^x <> i then
            Add(movedPoints, i);
        fi;
    od;
    cyc := [];
    cur := Minimum(movedPoints);
    for i in [1..5] do
        Add(cyc, cur);
        cur := cur^x;
    od;
    cr := 0;
    for i in [1..5] do
        u := cyc[i];
        if i < 5 then
            v := cyc[i+1];
        else
            v := cyc[1];
        fi;
        if (u in A and v in BB) or (u in BB and v in A) then
            cr := cr + 1;
        fi;
    od;
    if x in C51 then
        if cr = 2 then
            Add(C51cr2, x);
        elif cr = 4 then
            Add(C51cr4, x);
        else
            Print("   Unexpected cr value for ", x, ": ", cr, "\n");
        fi;
    else
        if cr = 2 then
            Add(C52cr2, x);
        elif cr = 4 then
            Add(C52cr4, x);
        else
            Print("   Unexpected cr value for ", x, ": ", cr, "\n");
        fi;
    fi;
od;

D1 := Union(C51cr2, C52cr4);
D2 := Union(C51cr4, C52cr2);

Print("   |C_(5)_1 with cr=2| = ", Length(C51cr2), "\n");
Print("   |C_(5)_1 with cr=4| = ", Length(C51cr4), "\n");
Print("   |C_(5)_2 with cr=2| = ", Length(C52cr2), "\n");
Print("   |C_(5)_2 with cr=4| = ", Length(C52cr4), "\n");
Print("   |D_1| = ", Length(D1), ", |D_2| = ", Length(D2), "\n");
Print("   |D_1 cap D_2| = ", Length(Intersection(D1, D2)), "\n");
Print("   |D_1 union D_2| = ", Length(Union(D1, D2)), "\n");
Print("   Number of inverses of D_1 outside D_1 = ",
      Length(Filtered(D1, g -> not (g^-1 in D1))), "\n");
Print("   Number of inverses of D_2 outside D_2 = ",
      Length(Filtered(D2, g -> not (g^-1 in D2))), "\n\n");

# ---------------------------------------------------------------------------
# The Sylow-5 description of D_1.
# The program checks that the 18 listed generators generate 18 distinct
# Sylow 5-subgroups, and that their nonidentity elements are exactly D_1.
# ---------------------------------------------------------------------------

Print("6. Sylow-5 description of D_1\n");

GensPlus := [
    (1,2,3,4,5), (1,2,3,5,6), (1,2,3,6,4),
    (1,2,4,3,5), (1,2,4,5,3), (1,2,4,6,5),
    (1,2,5,3,6), (1,2,5,4,6), (1,2,5,6,3),
    (1,2,6,3,4), (1,2,6,4,3), (1,2,6,5,4),
    (1,3,4,5,6), (1,3,5,6,4), (1,3,6,4,5),
    (2,3,4,6,5), (2,3,5,4,6), (2,3,6,5,4)
];

sylowSets := [];
unionSylowNonId := [];
for g in GensPlus do
    H := Group(g);
    Add(sylowSets, Set(Elements(H)));
    for h in Elements(H) do
        if h <> e then
            AddSet(unionSylowNonId, h);
        fi;
    od;
od;

Print("   Number of listed generators in D_1 = ",
      Length(Filtered(GensPlus, g -> g in D1)), "\n");
Print("   Number of distinct generated subgroups = ",
      Length(Set(sylowSets)), "\n");
Print("   |D_1 minus Sylow union| = ",
      Length(Difference(D1, unionSylowNonId)), "\n");
Print("   |Sylow union minus D_1| = ",
      Length(Difference(unionSylowNonId, D1)), "\n");
Print("   Hence |D_1| from Sylow union = ", Length(unionSylowNonId), "\n\n");

# ---------------------------------------------------------------------------
# The twisted T = underline(D_1) - underline(D_2) identities.
# These are the concrete identities used to prove that the twisted partition
# is a Schur partition and that the two adjacency algebras are algebraically
# isomorphic.
# ---------------------------------------------------------------------------

Print("7. Twisted T identities\n");

twClasses := [C0, C22, C3, C33, C42, D1, D2];
twClassNames := ["0", "(2,2)", "(3)", "(3,3)", "(4,2)", "D_1", "D_2"];

signT := List(elts, g -> 0);
for g in D1 do
    signT[Position(elts, g)] := 1;
od;
for g in D2 do
    signT[Position(elts, g)] := -1;
od;

for pair in [[C22,"C_(2,2)"], [C42,"C_(4,2)"],
             [C3,"C_(3)"], [C33,"C_(3,3)"]] do
    coeff := List(elts, g -> 0);
    for a in pair[1] do
        for b in All5 do
            coeff[Position(elts, a*b)] :=
                coeff[Position(elts, a*b)] + signT[Position(elts, b)];
        od;
    od;
    Print("   ", pair[2], " * T value sets on classes ",
          twClassNames, " = ", ValueSetsOnParts(coeff, twClasses), "\n");
od;

coeff := List(elts, g -> 0);
for term in [[D1, D1, 1], [D2, D2, 1],
             [D1, D2, -1], [D2, D1, -1]] do
    for a in term[1] do
        for b in term[2] do
            coeff[Position(elts, a*b)] :=
                coeff[Position(elts, a*b)] + term[3];
        od;
    od;
od;
Print("   T^2 value sets on classes ", twClassNames, " = ",
      ValueSetsOnParts(coeff, twClasses), "\n\n");

# ---------------------------------------------------------------------------
# Local graph Gamma_Y(X):
#   vertices are elements x in X,
#   x is adjacent to y when y*x^-1 is in Y.
# All relevant basic sets are inverse-closed, so this is an undirected graph.
# ---------------------------------------------------------------------------

Print("8. Local graphs in the original and twisted partitions\n");

graphC51C3 := Graph(Group(()), C51,
    function(x, g) return x; end,
    function(x, y) return x <> y and y*x^-1 in C3; end,
    true);
componentSizes := List(ConnectedComponents(graphC51C3), Length);
Sort(componentSizes);
Print("   Original Gamma_C(3)(C_(5)_1): degree set = ",
      VertexDegrees(graphC51C3), ", component sizes = ", componentSizes, "\n");

IcoNeighbors := [
    [2,3,4,5,6],
    [1,3,6,7,8],
    [1,2,4,8,9],
    [1,3,5,9,10],
    [1,4,6,10,11],
    [1,2,5,7,11],
    [12,2,6,8,11],
    [12,2,3,7,9],
    [12,3,4,8,10],
    [12,4,5,9,11],
    [12,5,6,7,10],
    [7,8,9,10,11]
];
IcoGraph := Graph(Group(()), [1..12],
    function(x, g) return x; end,
    function(x, y) return y in IcoNeighbors[x]; end,
    true);
componentGraphs := List(ConnectedComponents(graphC51C3),
    comp -> InducedSubgraph(graphC51C3, comp));
icoChecks := List(componentGraphs, g -> IsIsomorphicGraph(g, IcoGraph));
Print("   Componentwise isomorphic to the icosahedral graph = ",
      icoChecks, "\n");

Print("   Neighbors of x0 in original Gamma_C(3)(C_(5)_1):\n");
for j in graphC51C3.adjacencies[Position(graphC51C3.names, x0)] do
    Print("      ", graphC51C3.names[j], "\n");
od;

candidatePairs := [[D1,C3,"Gamma_C(3)(D_1)"],
                   [D1,C33,"Gamma_C(3,3)(D_1)"],
                   [D2,C3,"Gamma_C(3)(D_2)"],
                   [D2,C33,"Gamma_C(3,3)(D_2)"]];

for pair in candidatePairs do
    graphLocal := Graph(Group(()), pair[1],
        function(x, g) return x; end,
        function(x, y) return x <> y and y*x^-1 in pair[2]; end,
        true);
    componentSizes := List(ConnectedComponents(graphLocal), Length);
    Sort(componentSizes);
    Print("   ", pair[3], ": degree set = ",
          VertexDegrees(graphLocal), ", component sizes = ", componentSizes, "\n");
od;
Print("   Thus no twisted candidate local graph has six 12-vertex components.\n\n");

# ---------------------------------------------------------------------------
# Non-Schurian witness in the text:
# For Gamma_C(3)(D_1), count the number of triangles containing each vertex.
# If the graph were vertex-transitive, this number would be constant.
# ---------------------------------------------------------------------------

Print("9. Triangle counts in Gamma_C(3)(D_1)\n");

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
Print("   Number of vertices with tau=1 = ",
      Length(Filtered(triangleCounts, t -> t = 1)), "\n");
Print("   Number of vertices with tau=3 = ",
      Length(Filtered(triangleCounts, t -> t = 3)), "\n");
Print("   Triangle-count distribution = ",
      Collected(triangleCounts), "\n\n");

# ---------------------------------------------------------------------------
# The remark about choosing any 3-subset A:
# A and its complement define the same cut, so there are 10 cuts.  We verify
# that each cut gives the same T-identities as above, and store the resulting
# partitions for the final exhaustive comparison.
# ---------------------------------------------------------------------------

Print("10. The ten 3+3 cuts\n");

all3sets := Combinations([1..6], 3);
canonicalCuts := [];
for Sset in all3sets do
    comp := Difference([1..6], Sset);
    if Minimum(Sset) < Minimum(comp) then
        Add(canonicalCuts, Sset);
    fi;
od;
Print("   Number of 3+3 cuts up to complement = ", Length(canonicalCuts), "\n");

cutsOKCount := 0;
cutConditionDefects := [];
cutPartitions := [];
for Sset in canonicalCuts do
    localDefects := 0;
    Tset := Difference([1..6], Sset);
    localC51cr2 := [];
    localC51cr4 := [];
    localC52cr2 := [];
    localC52cr4 := [];

    for x in All5 do
        movedPoints := [];
        for i in [1..6] do
            if i^x <> i then
                Add(movedPoints, i);
            fi;
        od;
        cyc := [];
        cur := Minimum(movedPoints);
        for i in [1..5] do
            Add(cyc, cur);
            cur := cur^x;
        od;
        cr := 0;
        for i in [1..5] do
            u := cyc[i];
            if i < 5 then
                v := cyc[i+1];
            else
                v := cyc[1];
            fi;
            if (u in Sset and v in Tset) or (u in Tset and v in Sset) then
                cr := cr + 1;
            fi;
        od;
        if x in C51 then
            if cr = 2 then
                Add(localC51cr2, x);
            else
                Add(localC51cr4, x);
            fi;
        else
            if cr = 2 then
                Add(localC52cr2, x);
            else
                Add(localC52cr4, x);
            fi;
        fi;
    od;

    localD1 := Union(localC51cr2, localC52cr4);
    localD2 := Union(localC51cr4, localC52cr2);
    Add(cutPartitions, [localD1, localD2]);
    localSignT := List(elts, g -> 0);
    for g in localD1 do
        localSignT[Position(elts, g)] := 1;
    od;
    for g in localD2 do
        localSignT[Position(elts, g)] := -1;
    od;

    for pair in [[C22,0], [C42,0], [C3,-5], [C33,-5]] do
        coeff := List(elts, g -> 0);
        for a in pair[1] do
            for b in All5 do
                coeff[Position(elts, a*b)] :=
                    coeff[Position(elts, a*b)] + localSignT[Position(elts, b)];
            od;
        od;
        for pos in [1..Length(elts)] do
            if coeff[pos] <> pair[2] * localSignT[pos] then
                localDefects := localDefects + 1;
            fi;
        od;
    od;

    coeff := List(elts, g -> 0);
    for a in All5 do
        for b in All5 do
            coeff[Position(elts, a*b)] :=
                coeff[Position(elts, a*b)] +
                localSignT[Position(elts, a)] * localSignT[Position(elts, b)];
        od;
    od;
    for pos in [1..Length(elts)] do
        if coeff[pos] <> sSquareCoeff[pos] then
            localDefects := localDefects + 1;
        fi;
    od;
    Add(cutConditionDefects, [Sset, localDefects]);
    if localDefects = 0 then
        cutsOKCount := cutsOKCount + 1;
    fi;
od;

Print("   Cut condition defect counts [3-subset, coefficient defects] = ",
      cutConditionDefects, "\n");
Print("   Number of 3+3 cuts satisfying all T-identities = ",
      cutsOKCount, "\n\n");

# ---------------------------------------------------------------------------
# Exhaustive check for the remark in the paper.
# We form the linear equations for a balanced +/-1 sign vector on All5: it
# must satisfy the same C_r*S=lambda_r*S equations as the original S.  We then
# enumerate the +/-1 vectors in this linear solution space, discard the
# original C_(5)_1/C_(5)_2 split, impose the same quadratic S^2 condition, and
# compare the resulting partitions with the ten 3+3 cuts and with the A6-orbit
# of one cut.
# ---------------------------------------------------------------------------

Print("11. Exhaustive check of nontrivial A6 partitions\n");

linearConditionMatrix := [];
Add(linearConditionMatrix, List(All5, b -> 1));
for actionData in sActionOnFixedClasses do
    activeClass := actionData[1];
    lambda := actionData[2];
    for g in elts do
        row := [];
        for b in All5 do
            entry := 0;
            if g * b^-1 in activeClass then
                entry := entry + 1;
            fi;
            if g = b then
                entry := entry - lambda;
            fi;
            Add(row, entry);
        od;
        if ForAny(row, x -> x <> 0) then
            Add(linearConditionMatrix, row);
        fi;
    od;
od;

basis := NullspaceMat(TransposedMat(linearConditionMatrix));
signVectors := [];
d := Length(basis);
if d > 0 then
    # Pivot coordinates of a row-space basis give independent coordinates;
    # choosing their signs determines a unique vector in the linear space.
    rowSpaceBasis := BaseMat(basis);
    pivots := List(rowSpaceBasis, PositionNonZero);

    rows := [];
    for basisVec in basis do
        Add(rows, basisVec{pivots});
    od;
    for signs in Tuples([1,-1], d) do
        coeffs := SolutionMat(rows, signs);
        vec := coeffs * basis;
        if vec[1] = 1 and ForAll(vec, x -> x = 1 or x = -1) and Sum(vec) = 0 then
            Add(signVectors, vec);
        fi;
    od;
fi;

foundPartitions := [];
for vec in signVectors do
    all5WithSigns := TransposedMat([All5, vec]);
    Dp := List(Filtered(all5WithSigns, pair -> pair[2] = 1), pair -> pair[1]);
    Dm := List(Filtered(all5WithSigns, pair -> pair[2] = -1), pair -> pair[1]);
    if Set([Set(Dp), Set(Dm)]) <> Set([Set(C51), Set(C52)]) then
        candidateSplitSign := List(elts, x -> 0);
        for a in Dp do
            candidateSplitSign[LookupDictionary(posG, a)] := 1;
        od;
        for a in Dm do
            candidateSplitSign[LookupDictionary(posG, a)] := -1;
        od;

        conditionOK := true;
        for actionData in sActionOnFixedClasses do
            activeClass := actionData[1];
            lambda := actionData[2];
            vNew := List(elts, x -> 0);
            for a in activeClass do
                for b in All5 do
                    vNew[LookupDictionary(posG, a*b)] :=
                        vNew[LookupDictionary(posG, a*b)]
                        + candidateSplitSign[LookupDictionary(posG, b)];
                od;
            od;

            for idx in [1..Length(elts)] do
                if vNew[idx] <> lambda * candidateSplitSign[idx] then
                    conditionOK := false;
                    break;
                fi;
            od;
            if not conditionOK then
                break;
            fi;
        od;

        if conditionOK then
            sqNew := List(elts, x -> 0);
            for a in All5 do
                for b in All5 do
                    sqNew[LookupDictionary(posG, a*b)] :=
                        sqNew[LookupDictionary(posG, a*b)]
                        + candidateSplitSign[LookupDictionary(posG, a)]
                        * candidateSplitSign[LookupDictionary(posG, b)];
                od;
            od;
            if sSquareCoeff = sqNew then
                Add(foundPartitions, [Dp, Dm]);
            fi;
        fi;
    fi;
od;

cutPartitionSets := Set(List(cutPartitions,
    pair -> Set([Set(pair[1]), Set(pair[2])])));
foundPartitionSets := Set(List(foundPartitions,
    pair -> Set([Set(pair[1]), Set(pair[2])])));

Print("   Linear nullity for the split 5-cycle pair = ", Length(basis), "\n");
Print("   Balanced sign vectors in the linear solution space = ",
      Length(signVectors), "\n");
Print("   Nontrivial partitions satisfying all conditions = ",
      Length(foundPartitions), "\n");
Print("   Number of cut partitions from the 3+3 cut search = ",
      Length(cutPartitionSets), "\n");
Print("   Number of found partitions that are 3+3 cut partitions = ",
      Length(Intersection(foundPartitionSets, cutPartitionSets)), "\n");

basePair := cutPartitions[1];
baseOrbit := [];
for sigma in elts do
    AddSet(baseOrbit, Set([
        Set(List(basePair[1], x -> x^sigma)),
        Set(List(basePair[2], x -> x^sigma))
    ]));
od;
Print("   A6-orbit size of the first cut partition = ",
      Length(baseOrbit), "\n");
Print("   Number of cut partitions in that A6-orbit = ",
      Length(Intersection(baseOrbit, cutPartitionSets)), "\n\n");

# Write the {1,2,3}|{4,5,6} twisted Schur partition as GAP input for
# later use by checktwistA6.g.
PrintTo(Concatenation(OutDir, "/A6_twisted_Schur_partition.txt"),
        "A6TwistedSchurPartition := ", twClasses, ";\n");

Print("Done.\n");

LogTo();
