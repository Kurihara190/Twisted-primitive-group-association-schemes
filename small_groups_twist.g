# GAP search for the two-class twisting criterion.
#
# The search bound is set to 200.
# For every group in the search range,
# this script finds the pairs of conjugacy classes satisfying (C1)--(C4).  For
# each pair it exhaustively enumerates unordered, non-original partitions
# {D_1,D_2} satisfying (D1)--(D4).  Conditions (D1) and (D3) are imposed as
# linear equations on the sign vector of T, condition (D2) is checked under
# inversion, and condition (D4), T^2=S^2, is verified exactly.
#
# Detailed results are written as one row per (C1)--(C4) class pair to
# result_small_groups_twist.csv.  The CSV records the linear-search statistics
# and the combinatorial-isomorphism outcome.  The text file
# result_small_groups_twist.txt contains only the final search summary.

OutDir := "out";
if IsExistingFile("code/small_groups_twist.g") then
    OutDir := "code/out";
fi;
Exec(Concatenation("mkdir -p ", OutDir));

if LoadPackage("SmallGrp") = fail then
    Error("The GAP package SmallGrp is required for the order-at-most-200 search.");
fi;
if LoadPackage("AssociationSchemes") = fail then
    Error("The GAP package AssociationSchemes is required for scheme isomorphism checks.");
fi;

SummaryFile := Concatenation(OutDir, "/result_small_groups_twist.txt");
CsvFile := Concatenation(OutDir, "/result_small_groups_twist.csv");

# Quote a CSV text field, doubling embedded quotation marks.  Text fields are
# always quoted; numeric fields are written directly so that spreadsheet and
# statistics software read them as numbers.
CsvQuote := function(value)
    local text;
    text := String(value);
    text := ReplacedString(text, "\r", " ");
    text := ReplacedString(text, "\n", " ");
    return Concatenation("\"", ReplacedString(text, "\"", "\"\""), "\"");
end;

CsvStream := OutputTextFile(CsvFile, false);
SetPrintFormattingStatus(CsvStream, false);
PrintTo(CsvStream,
    "group_order,small_group_id,group_structure,class_i,class_j,class_size,",
    "union_size,linear_nullity,unordered_sign_solutions,",
    "nontrivial_partitions,isomorphic_partitions,new_twists,isomorphism_flags\n");

ClassPairConditionsHold := function(classes, i, j)
    local U, inversePair, r, k, a, b, productCoeff, lambda,
          sSquareCoeff, sa, sb;
    U := Union(classes[i], classes[j]);

    # (C1): the two conjugacy classes have the same size.
    if Length(classes[i]) <> Length(classes[j]) then return false; fi;

    # (C2): inversion preserves the unordered pair of classes.
    inversePair := Set([
        Set(List(classes[i], x -> x^-1)),
        Set(List(classes[j], x -> x^-1))
    ]);
    if inversePair <> Set([Set(classes[i]), Set(classes[j])]) then
        return false;
    fi;

    # (C3): for S=underline(C_i)-underline(C_j), each unchanged class sum
    # satisfies underline(C_r) S = lambda_r S.
    for r in Difference([1..Length(classes)], [i,j]) do
        productCoeff := List([1..Length(classes)], t -> 0);
        for a in classes[r] do
            for b in classes[i] do
                for k in [1..Length(classes)] do
                    if a*b in classes[k] then
                        productCoeff[k] := productCoeff[k] + 1;
                        break;
                    fi;
                od;
            od;
            for b in classes[j] do
                for k in [1..Length(classes)] do
                    if a*b in classes[k] then
                        productCoeff[k] := productCoeff[k] - 1;
                        break;
                    fi;
                od;
            od;
        od;
        for k in [1..Length(classes)] do
            productCoeff[k] := productCoeff[k] / Length(classes[k]);
        od;
        lambda := productCoeff[i];
        if productCoeff[j] <> -lambda then return false; fi;
        for k in Difference([1..Length(classes)], [i,j]) do
            if productCoeff[k] <> 0 then return false; fi;
        od;
    od;

    # (C4): S^2 is central automatically.  It lies in the span of the
    # unchanged class sums and U=underline(C_i)+underline(C_j) exactly when
    # its coefficients on C_i and C_j agree.
    sSquareCoeff := List([1..Length(classes)], t -> 0);
    for a in U do
        for b in U do
            for k in [1..Length(classes)] do
                if a*b in classes[k] then
                    if a in classes[i] then
                        sa := 1;
                    else
                        sa := -1;
                    fi;
                    if b in classes[i] then
                        sb := 1;
                    else
                        sb := -1;
                    fi;
                    sSquareCoeff[k] := sSquareCoeff[k] + sa * sb;
                    break;
                fi;
            od;
        od;
    od;
    for k in [1..Length(classes)] do
        sSquareCoeff[k] := sSquareCoeff[k] / Length(classes[k]);
    od;
    return sSquareCoeff[i] = sSquareCoeff[j];
end;

# The relation indexed by a part X consists of the pairs (x,y) for which
# y*x^-1 lies in X, matching the convention used in the paper.
RelationMatrixFromPartition := function(G, part)
    local elts, M, i, j, k, prod;
    elts := Elements(G);
    M := [];
    for i in [1..Length(elts)] do
        M[i] := [];
        for j in [1..Length(elts)] do
            prod := elts[j] * elts[i]^-1;
            for k in [1..Length(part)] do
                if prod in part[k] then
                    M[i][j] := k - 1;
                    break;
                fi;
            od;
        od;
    od;
    return M;
end;

ReplacementConditionsHold := function(G, classes, i, j, D1, D2)
    local U, inversePair, signS, signT, elts, pos, r, idx, a, b,
          sAction, tAction, sSquareCoeff, tSquareCoeff, lambda;
    elts := Elements(G);
    pos := NewDictionary(One(G), true);
    for idx in [1..Length(elts)] do AddDictionary(pos, elts[idx], idx); od;
    U := Union(classes[i], classes[j]);

    # (D1), together with the requirement that D_1,D_2 partition U.
    if Length(D1) <> Length(D2)
       or Intersection(D1, D2) <> []
       or Set(Union(D1, D2)) <> Set(U) then
        return false;
    fi;

    # (D2): inversion preserves the unordered replacement pair.
    inversePair := Set([
        Set(List(D1, x -> x^-1)),
        Set(List(D2, x -> x^-1))
    ]);
    if inversePair <> Set([Set(D1), Set(D2)]) then return false; fi;

    signS := List(elts, x -> 0);
    signT := List(elts, x -> 0);
    for a in classes[i] do signS[LookupDictionary(pos,a)] := 1; od;
    for a in classes[j] do signS[LookupDictionary(pos,a)] := -1; od;
    for a in D1 do signT[LookupDictionary(pos,a)] := 1; od;
    for a in D2 do signT[LookupDictionary(pos,a)] := -1; od;

    # (D3): every unchanged class sum acts on T by the same scalar by which
    # it acts on S.  Recompute all coefficients here after the linear filter.
    for r in Difference([1..Length(classes)], [i,j]) do
        sAction := List(elts, x -> 0);
        tAction := List(elts, x -> 0);
        for a in classes[r] do
            for b in U do
                idx := LookupDictionary(pos, a*b);
                sAction[idx] := sAction[idx] + signS[LookupDictionary(pos,b)];
                tAction[idx] := tAction[idx] + signT[LookupDictionary(pos,b)];
            od;
        od;
        idx := LookupDictionary(pos, classes[i][1]);
        lambda := sAction[idx];
        for idx in [1..Length(elts)] do
            if tAction[idx] <> lambda * signT[idx] then return false; fi;
        od;
    od;

    # (D4): compare the coefficients of T^2 and S^2 on every group element.
    sSquareCoeff := List(elts, x -> 0);
    tSquareCoeff := List(elts, x -> 0);
    for a in U do
        for b in U do
            idx := LookupDictionary(pos, a*b);
            sSquareCoeff[idx] := sSquareCoeff[idx]
                + signS[LookupDictionary(pos,a)] * signS[LookupDictionary(pos,b)];
            tSquareCoeff[idx] := tSquareCoeff[idx]
                + signT[LookupDictionary(pos,a)] * signT[LookupDictionary(pos,b)];
        od;
    od;
    return sSquareCoeff = tSquareCoeff;
end;

LinearReplacementBasis := function(G, classes, i, j)
    local elts, pos, U, signOld, matrix, r, lambda, idx, g, row, col, a, b;

    elts := Elements(G);
    pos := NewDictionary(One(G), true);
    for idx in [1..Length(elts)] do AddDictionary(pos, elts[idx], idx); od;

    U := Union(classes[i], classes[j]);
    signOld := List(elts, x -> 0);
    for a in classes[i] do signOld[LookupDictionary(pos,a)] := 1; od;
    for a in classes[j] do signOld[LookupDictionary(pos,a)] := -1; od;

    matrix := [];

    # Condition (D1): |D_1|=|D_2|.
    Add(matrix, List([1..Length(U)], t -> 1));

    # For every unchanged conjugacy class C_r, impose (D3):
    # underline(C_r) T = lambda_r T, including coefficient equations on G\U.
    # The corresponding equation for U follows from (D1) and these equations.
    for r in Difference([1..Length(classes)], [i,j]) do
        lambda := 0;
        g := classes[i][1];
        for a in classes[r] do
            b := a^-1 * g;
            if b in U then
                lambda := lambda + signOld[LookupDictionary(pos,b)];
            fi;
        od;

        for g in elts do
            row := List([1..Length(U)], t -> 0);
            for col in [1..Length(U)] do
                b := U[col];
                if g * b^-1 in classes[r] then
                    row[col] := row[col] + 1;
                fi;
                if g = b then
                    row[col] := row[col] - lambda;
                fi;
            od;
            if ForAny(row, x -> x <> 0) then Add(matrix, row); fi;
        od;
    od;

    # NullspaceMat gives a left nullspace.  Transpose to get column solutions
    # of matrix * x = 0, represented as ordinary row vectors.
    return NullspaceMat(TransposedMat(matrix));
end;

UnorderedSignVectorsInSpan := function(basis)
    local d, m, pivots, col, test, rows, signs, mask, tmp, coeffs, v, out, i;

    d := Length(basis);
    if d = 0 then return []; fi;
    m := Length(basis[1]);

    # Pick coordinates that determine the coefficients in the nullspace basis.
    pivots := [];
    for col in [1..m] do
        test := Concatenation(pivots, [col]);
        rows := List(basis, v -> v{test});
        if RankMat(rows) = Length(test) then
            pivots := test;
            if Length(pivots) = d then break; fi;
        fi;
    od;

    out := [];
    for mask in [0..2^d-1] do
        tmp := mask;
        signs := [];
        for i in [1..d] do
            if tmp mod 2 = 0 then Add(signs, 1); else Add(signs, -1); fi;
            tmp := QuoInt(tmp, 2);
        od;

        rows := List(basis, v -> v{pivots});
        coeffs := SolutionMat(rows, signs);
        v := List([1..m], col -> Sum([1..d], i -> coeffs[i] * basis[i][col]));

        # The vectors v and -v give the same unordered partition {D_1,D_2}.
        # Count it once by keeping the orientation with first coordinate +1.
        if v[1] = 1 and ForAll(v, x -> x = 1 or x = -1) and Sum(v) = 0 then
            Add(out, v);
        fi;
    od;
    return out;
end;

classPairRecords := [];
nontrivialPartitionRecords := [];
groupsChecked := 0;

for n in [1..200] do
    numberOfGroups := NumberSmallGroups(n);
    if numberOfGroups = fail then
        Error("Small Groups Library data are unavailable for order ", n, ".");
    fi;
    for id in [1..numberOfGroups] do
        groupsChecked := groupsChecked + 1;
        G := SmallGroup(n,id);
        groupDescription := StructureDescription(G);
        classes := List(ConjugacyClasses(G), AsList);
        goodPairs := [];
        for i in [1..Length(classes)] do
            for j in [i+1..Length(classes)] do
                if ClassPairConditionsHold(classes, i, j) then
                    Add(goodPairs, [i,j,Length(classes[i])]);
                fi;
            od;
        od;
        if goodPairs <> [] then
            Add(classPairRecords, [n,id,groupDescription,goodPairs]);
        fi;

        for pair in goodPairs do
            i := pair[1]; j := pair[2];
            U := Union(classes[i], classes[j]);
            found := [];
            basis := LinearReplacementBasis(G, classes, i, j);
            signVectors := UnorderedSignVectorsInSpan(basis);
            for vec in signVectors do
                D1 := U{Filtered([1..Length(U)], t -> vec[t] = 1)};
                D2 := U{Filtered([1..Length(U)], t -> vec[t] = -1)};
                if Set([Set(D1),Set(D2)])
                   <> Set([Set(classes[i]),Set(classes[j])])
                   and ReplacementConditionsHold(G, classes, i, j, D1, D2) then
                    Add(found, [D1,D2]);
                fi;
            od;
            isoList := [];
            if found <> [] then
                originalPart := ShallowCopy(classes);
                originalScheme := HomogeneousCoherentConfiguration(
                    RelationMatrixFromPartition(G, originalPart));
                for entry in found do
                    replacementPart := ShallowCopy(classes);
                    replacementPart[i] := entry[1];
                    replacementPart[j] := entry[2];
                    replacementScheme := HomogeneousCoherentConfiguration(
                        RelationMatrixFromPartition(G, replacementPart));
                    Add(isoList, AreIsomorphicHomogeneousCoherentConfigurations(
                        originalScheme, replacementScheme));
                od;
                Add(nontrivialPartitionRecords,
                    [n,id,[i,j],Length(found),isoList]);
            fi;
            isomorphicCount := Length(Filtered(isoList, flag -> flag = true));
            newTwistCount := Length(Filtered(isoList, flag -> flag = false));
            PrintTo(CsvStream,
                n, ",", id, ",", CsvQuote(groupDescription), ",",
                i, ",", j, ",", pair[3], ",", Length(U), ",",
                Length(basis), ",", Length(signVectors), ",", Length(found), ",",
                isomorphicCount, ",", newTwistCount, ",",
                CsvQuote(JoinStringsWithSeparator(List(isoList, String), ";")),
                "\n");
        od;
    od;
od;

CloseStream(CsvStream);

LogTo(SummaryFile);
Print("Small-group search summary\n");
Print("   Maximum group order examined = 200\n");
Print("   Small groups examined = ", groupsChecked, "\n");
Print("   Groups with at least one (C1)--(C4) class pair = ",
      Length(classPairRecords), "\n");
Print("   Nontrivial (non-original) partition records\n");
Print("   [order, id, pair, number, combinatorially-isomorphic flags] =\n   ",
      nontrivialPartitionRecords, "\n");
Print("   Total nontrivial (D1)--(D4) partitions = ",
      Sum(List(nontrivialPartitionRecords, entry -> entry[4])), "\n");
Print("   Number not combinatorially isomorphic to the original (new twists) = ",
      Sum(List(nontrivialPartitionRecords,
               entry -> Length(Filtered(entry[5], flag -> flag = false)))), "\n");

LogTo();
