# GAP verification for the section:
#   "Rigidity for A4"
#
# The script checks the class data and the spectral assertions used to prove
# that no nontrivial split-class twist occurs for A4.

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
if LoadPackage("grape") = fail then
    Error("The GAP package GRAPE is required for local graph calculations.");
fi;
LogTo(Concatenation(OutDir, "/result_A4.txt"));

Print("\n");
Print("============================================================\n");
Print("Rigidity checks for A4\n");
Print("============================================================\n\n");

G := AlternatingGroup(4);
e := One(G);

C0 := [e];
C22 := AsList(ConjugacyClass(G, (1,2)(3,4)));
C31 := AsList(ConjugacyClass(G, (1,2,3)));
C32 := AsList(ConjugacyClass(G, (1,3,2)));
OmegaSet := Concatenation(C31, C32);

Print("Class sizes [C0,C22,C3_1,C3_2] = ",
      List([C0,C22,C31,C32], Length), "\n");

# The graph G_Omega has vertex set Omega and edge x--y when y*x^-1 is in Omega.
graphOmega := Graph(Group(()), OmegaSet,
    function(x, g) return x; end,
    function(x, y) return x <> y and y*x^-1 in OmegaSet; end,
    true);
A := List([1..graphOmega.order], i -> List([1..graphOmega.order], j -> 0));
for i in [1..graphOmega.order] do
    for j in Adjacency(graphOmega, i) do
        A[i][j] := 1;
    od;
od;

Print("G_Omega degree set = ", VertexDegrees(graphOmega), "\n");

K44 := Graph(Group(()), [1..8],
    function(x, g) return x; end,
    function(x, y) return (x <= 4 and y > 4) or (x > 4 and y <= 4); end,
    true);
Print("G_Omega is K_{4,4}? ", IsIsomorphicGraph(graphOmega, K44), "\n");

cp := CharacteristicPolynomial(Rationals, Rationals, A);
Print("Characteristic polynomial of A_Omega = ", cp, "\n");
Print("Factors of characteristic polynomial = ", Factors(cp), "\n");

chiS := Concatenation(List(C31, x -> 1), List(C32, x -> -1));
Av := A * chiS;
Print("A_Omega * chi_S = ", Av, "\n");
Print("dim Ker(A_Omega + 4I) = ",
      Length(NullspaceMat(A + 4 * IdentityMat(Length(OmegaSet)))), "\n\n");

Print("Done.\n");

LogTo();
