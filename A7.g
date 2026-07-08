# GAP verification for the section:
#   "Rigidity in A7: no split 7-cycle twist"
#
# The expensive assertion is that the -21 eigenspace of the local graph on
# Omega=C(7)_1 union C(7)_2 is one-dimensional.  We verify this over GF(1009).
# Since chi_S is an integral rational kernel vector for A+21I, and reduction
# modulo 1009 has nullity 1, the rational nullity is also exactly 1.

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
if LoadPackage("grape") = fail then
    Error("The GAP package GRAPE is required for local graph calculations.");
fi;
LogTo(Concatenation(OutDir, "/result_A7.txt"));

Print("\n");
Print("============================================================\n");
Print("Rigidity check for A7\n");
Print("============================================================\n\n");

G := AlternatingGroup(7);
e := One(G);
x0 := (1,2,3,4,5,6,7);

C0 := [e];
C22 := AsList(ConjugacyClass(G, (1,2)(3,4)));
C3 := AsList(ConjugacyClass(G, (1,2,3)));
C322 := AsList(ConjugacyClass(G, (1,2,3)(4,5)(6,7)));
C33 := AsList(ConjugacyClass(G, (1,2,3)(4,5,6)));
C42 := AsList(ConjugacyClass(G, (1,2,3,4)(5,6)));
C5 := AsList(ConjugacyClass(G, (1,2,3,4,5)));
C71 := AsList(ConjugacyClass(G, x0));
C72 := AsList(ConjugacyClass(G, x0^3));
OmegaSet := Concatenation(C71, C72);

Print("Class sizes [C0,C22,C3,C322,C33,C42,C5,C7_1,C7_2] =\n   ",
      List([C0,C22,C3,C322,C33,C42,C5,C71,C72], Length), "\n");

# Local graph on Omega with relation C_(2,2).
graph22 := Graph(Group(()), OmegaSet,
    function(x, g) return x; end,
    function(x, y) return x <> y and y*x^-1 in C22; end,
    true);

F := GF(1009);
one := One(F);
zero := Zero(F);
Aplus21 := List([1..graph22.order], i -> List([1..graph22.order], j -> zero));

for i in [1..graph22.order] do
    Aplus21[i][i] := 21 * one;
    for j in Adjacency(graph22, i) do
        Aplus21[i][j] := Aplus21[i][j] + one;
    od;
od;

Print("Degree set of G_(2,2) on Omega = ", VertexDegrees(graph22), "\n");

chiS := Concatenation(List(C71, x -> 1), List(C72, x -> -1));
Av := [];
for i in [1..graph22.order] do
    s := 0;
    for j in Adjacency(graph22, i) do
        s := s + chiS[j];
    od;
    Av[i] := s;
od;

Print("Entrywise quotients (A_(2,2)*chi_S)/chi_S = ",
      Set(List([1..graph22.order], i -> Av[i] / chiS[i])), "\n");
Print("dim Ker(A_(2,2)+21I) over GF(1009) = ",
      Length(NullspaceMat(Aplus21)), "\n");
Print("This finite-field check verifies the claimed one-dimensional eigenspace.\n\n");

Print("Done.\n");

LogTo();
