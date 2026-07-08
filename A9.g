# GAP verification for the section:
#   "Rigidity in A9"
#
# This script uses the character table and class multiplication coefficients.
# It verifies the class sizes, the class-algebra conditions for the split
# classes, the listed same-size non-example, and the eigenvalue action on the
# original difference vector.
#
# At the end it also checks the centralizer-orbit quotient of the two local
# graphs used in the rigidity argument.  This is much smaller than the full
# local graph and is a useful exact consistency check over GF(101).

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
LogTo(Concatenation(OutDir, "/result_A9.txt"));

LoadPackage("cvec");

Print("\n");
Print("============================================================\n");
Print("Class-algebra checks for A9\n");
Print("============================================================\n\n");

G := AlternatingGroup(9);
e := One(G);
x9 := (1,2,3,4,5,6,7,8,9);
x9alt := (1,2,3,4,5,6,7,9,8);
y53 := (1,2,3,4,5)(6,7,8);

C0 := [e];
C22 := AsList(ConjugacyClass(G, (1,2)(3,4)));
C2222 := AsList(ConjugacyClass(G, (1,2)(3,4)(5,6)(7,8)));
C3 := AsList(ConjugacyClass(G, (1,2,3)));
C322 := AsList(ConjugacyClass(G, (1,2,3)(4,5)(6,7)));
C33 := AsList(ConjugacyClass(G, (1,2,3)(4,5,6)));
C333 := AsList(ConjugacyClass(G, (1,2,3)(4,5,6)(7,8,9)));
C42 := AsList(ConjugacyClass(G, (1,2,3,4)(5,6)));
C432 := AsList(ConjugacyClass(G, (1,2,3,4)(5,6,7)(8,9)));
C44 := AsList(ConjugacyClass(G, (1,2,3,4)(5,6,7,8)));
C5 := AsList(ConjugacyClass(G, (1,2,3,4,5)));
C522 := AsList(ConjugacyClass(G, (1,2,3,4,5)(6,7)(8,9)));
C531 := AsList(ConjugacyClass(G, y53));
C532 := AsList(ConjugacyClass(G, y53^7));
C62 := AsList(ConjugacyClass(G, (1,2,3,4,5,6)(7,8)));
C7 := AsList(ConjugacyClass(G, (1,2,3,4,5,6,7)));
C91 := AsList(ConjugacyClass(G, x9));
C92 := AsList(ConjugacyClass(G, x9alt));

classes := [C0, C22, C2222, C3, C322, C33, C333, C42, C432,
            C44, C5, C522, C531, C532, C62, C7, C91, C92];
classReps := [e, (1,2)(3,4), (1,2)(3,4)(5,6)(7,8), (1,2,3),
              (1,2,3)(4,5)(6,7), (1,2,3)(4,5,6),
              (1,2,3)(4,5,6)(7,8,9), (1,2,3,4)(5,6),
              (1,2,3,4)(5,6,7)(8,9), (1,2,3,4)(5,6,7,8),
              (1,2,3,4,5), (1,2,3,4,5)(6,7)(8,9),
              y53, y53^7, (1,2,3,4,5,6)(7,8),
              (1,2,3,4,5,6,7), x9, x9alt];
classLabels := ["C0", "C22", "C2222", "C3", "C322", "C33", "C333",
                "C42", "C432", "C44", "C5", "C522", "C531", "C532",
                "C62", "C7", "C91", "C92"];
ClassNo := function(C)
    return PositionProperty(classes, D -> D = C);
end;

ct := CharacterTable(G);
orders := List(classReps, Order);
sizes := List(classes, Length);
classCount := Length(classes);

if orders <> OrdersClassRepresentatives(ct)
   or sizes <> SizesConjugacyClasses(ct) then
    Error("The concrete A9 class order does not match the character table.");
fi;

Print("1. Character table class data\n");
Print("   Concrete classes = ", classLabels, "\n");
Print("   Orders      = ", orders, "\n");
Print("   Sizes       = ", sizes, "\n\n");

Print("2. Same-size non-example C_(3,2,2) and C_(4,2)\n");
i := ClassNo(C322);
j := ClassNo(C42);
r := ClassNo(C22);
v := List([1..classCount], k ->
    ClassMultiplicationCoefficient(ct, r, i, k)
    - ClassMultiplicationCoefficient(ct, r, j, k));
Print("   Coefficient vector for C_(2,2)*(C_(3,2,2)-C_(4,2)) = ",
      v, "\n");
sq := List([1..classCount], k ->
    ClassMultiplicationCoefficient(ct, i, i, k)
    - ClassMultiplicationCoefficient(ct, i, j, k)
    - ClassMultiplicationCoefficient(ct, j, i, k)
    + ClassMultiplicationCoefficient(ct, j, j, k));
Print("   S^2 coefficient on C_(3,2,2) = ", sq[i], "\n");
Print("   S^2 coefficient on C_(4,2)   = ", sq[j], "\n\n");

Print("3. Split 9-cycle classes\n");
p := ClassNo(C91);
q := ClassNo(C92);
rel := ClassNo(C22);
eig := ClassMultiplicationCoefficient(ct, rel, p, p)
       - ClassMultiplicationCoefficient(ct, rel, p, q);
Print("   Eigenvalue from C_(2,2) on chi_S = ", eig, "\n");

csData := [];
for rel2 in Difference([1..classCount], [p,q]) do
    v := List([1..classCount], k ->
        ClassMultiplicationCoefficient(ct, rel2, p, k)
        - ClassMultiplicationCoefficient(ct, rel2, q, k));
    lambda := v[p];
    Add(csData, [classLabels[rel2], lambda, v]);
od;
Print("   C*S coefficient data [class, lambda, coefficient vector] = ",
      csData, "\n");
sq := List([1..classCount], k ->
    ClassMultiplicationCoefficient(ct, p, p, k)
    - ClassMultiplicationCoefficient(ct, p, q, k)
    - ClassMultiplicationCoefficient(ct, q, p, k)
    + ClassMultiplicationCoefficient(ct, q, q, k));
Print("   S^2 coefficient vector = ", sq, "\n\n");

Print("4. Split (5,3)-classes\n");
p := ClassNo(C531);
q := ClassNo(C532);
rel := ClassNo(C3);
eig := ClassMultiplicationCoefficient(ct, rel, p, p)
       - ClassMultiplicationCoefficient(ct, rel, p, q);
Print("   Eigenvalue from C_(3) on chi_S = ", eig, "\n");

csData := [];
for rel2 in Difference([1..classCount], [p,q]) do
    v := List([1..classCount], k ->
        ClassMultiplicationCoefficient(ct, rel2, p, k)
        - ClassMultiplicationCoefficient(ct, rel2, q, k));
    lambda := v[p];
    Add(csData, [classLabels[rel2], lambda, v]);
od;
Print("   C*S coefficient data [class, lambda, coefficient vector] = ",
      csData, "\n");
sq := List([1..classCount], k ->
    ClassMultiplicationCoefficient(ct, p, p, k)
    - ClassMultiplicationCoefficient(ct, p, q, k)
    - ClassMultiplicationCoefficient(ct, q, p, k)
    + ClassMultiplicationCoefficient(ct, q, q, k));
Print("   S^2 coefficient vector = ", sq, "\n\n");

Print("5. Centralizer-orbit quotient checks over GF(101)\n");
Print("   These quotient matrices record the action on functions fixed by the\n");
Print("   centralizer of one vertex.  The character checks below verify that\n");
Print("   every irreducible constituent of the conjugation permutation module\n");
Print("   has a nonzero vector fixed by this centralizer, so a larger full\n");
Print("   eigenspace would also be visible in the quotient.\n\n");

prime := 101;


Print("   Split 9-cycle quotient for A_(2,2)+54I\n");
rep := x9;
H := Centralizer(G, rep);
OmegaSet := Concatenation(C91, C92);
rels := C22;
orbs := OrbitsDomain(H, OmegaSet, OnPoints);
Print("      |Omega| = ", Length(OmegaSet), ", |C_(2,2)| = ", Length(rels),
      ", |C_G(x)| = ", Size(H), ", quotient size = ", Length(orbs), "\n");
idx := NewDictionary(rep, true);
for a in [1..Length(orbs)] do
    for y in orbs[a] do AddDictionary(idx, y, a); od;
od;
rows := [];
for a in [1..Length(orbs)] do
    row := List([1..Length(orbs)], t -> 0);
    x := orbs[a][1];
    for r in rels do
        y := r^-1*x;
        b := LookupDictionary(idx, y);
        if b <> fail then row[b] := RemInt(row[b] + 1, prime); fi;
    od;
    row[a] := RemInt(row[a] + 54, prime);
    Add(rows, CVec(row, prime, 1));
od;
M := CMat(rows);
ns := NullspaceMat(M);
Print("      nullity over GF(", prime, ") = ", Length(ns), "\n\n");

Print("      Character visibility check for the full conjugation module\n");
permChar := List(classReps, g -> Number(OmegaSet, u -> u^g = u));
irr := Irr(ct);
mult := List(irr, chi -> ScalarProduct(ct, permChar, chi));
classOfH := List(Elements(H), h -> PositionProperty(classes, c -> h in c));
hFix := List(irr, chi -> Sum(classOfH, t -> chi[t]) / Size(H));
bad := Filtered([1..Length(irr)], t -> mult[t] > 0 and hFix[t] = 0);
Print("         irreducible constituents invisible to the quotient = ", bad, "\n");
Print("         certificate data [quotient nullity, invisible constituents] = ",
      [Length(ns), bad], "\n\n");

Print("   Split (5,3)-quotient for A_(3)+24I\n");
rep := y53;
H := Centralizer(G, rep);
OmegaSet := Concatenation(C531, C532);
rels := C3;
orbs := OrbitsDomain(H, OmegaSet, OnPoints);
Print("      |Omega| = ", Length(OmegaSet), ", |C_(3)| = ", Length(rels),
      ", |C_G(x)| = ", Size(H), ", quotient size = ", Length(orbs), "\n");
idx := NewDictionary(rep, true);
for a in [1..Length(orbs)] do
    for y in orbs[a] do AddDictionary(idx, y, a); od;
od;
rows := [];
for a in [1..Length(orbs)] do
    row := List([1..Length(orbs)], t -> 0);
    x := orbs[a][1];
    for r in rels do
        y := r^-1*x;
        b := LookupDictionary(idx, y);
        if b <> fail then row[b] := RemInt(row[b] + 1, prime); fi;
    od;
    row[a] := RemInt(row[a] + 24, prime);
    Add(rows, CVec(row, prime, 1));
od;
M := CMat(rows);
ns := NullspaceMat(M);
Print("      nullity over GF(", prime, ") = ", Length(ns), "\n\n");

Print("      Character visibility check for the full conjugation module\n");
permChar := List(classReps, g -> Number(OmegaSet, u -> u^g = u));
irr := Irr(ct);
mult := List(irr, chi -> ScalarProduct(ct, permChar, chi));
classOfH := List(Elements(H), h -> PositionProperty(classes, c -> h in c));
hFix := List(irr, chi -> Sum(classOfH, t -> chi[t]) / Size(H));
bad := Filtered([1..Length(irr)], t -> mult[t] > 0 and hFix[t] = 0);
Print("         irreducible constituents invisible to the quotient = ", bad, "\n");
Print("         certificate data [quotient nullity, invisible constituents] = ",
      [Length(ns), bad], "\n\n");

Print("Note: the rank computations are done on centralizer-orbit quotients.\n");
Print("      The character visibility checks show that, for these A9 modules,\n");
Print("      no nonzero eigenspace constituent can be hidden from the quotient.\n");
Print("Done.\n");

LogTo();
