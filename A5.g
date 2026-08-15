# GAP verification for the section:
#   "Rigidity for A5"
#
# The script checks the class data and the spectral assertions used to prove
# that no nontrivial split-class twist occurs for A5.

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));
if LoadPackage("grape") = fail then
    Error("The GAP package GRAPE is required for local graph calculations.");
fi;
LogTo(Concatenation(OutDir, "/result_A5.txt"));

Print("\n");
Print("============================================================\n");
Print("Rigidity checks for A5\n");
Print("============================================================\n\n");

G := AlternatingGroup(5);
e := One(G);

C0 := [e];
C22 := AsList(ConjugacyClass(G, (1,2)(3,4)));
C3 := AsList(ConjugacyClass(G, (1,2,3)));
C51 := AsList(ConjugacyClass(G, (1,2,3,4,5)));
C52 := AsList(ConjugacyClass(G, (1,3,5,2,4)));
OmegaSet := Concatenation(C51, C52);

Print("Class sizes [C0,C22,C3,C5_1,C5_2] = ",
      List([C0,C22,C3,C51,C52], Length), "\n");

# Graph on Omega with relation C_(2,2).
graph22 := Graph(Group(()), OmegaSet,
    function(x, g) return x; end,
    function(x, y) return x <> y and y*x^-1 in C22; end,
    true);
A := List([1..graph22.order], i -> List([1..graph22.order], j -> 0));
for i in [1..graph22.order] do
    for j in Adjacency(graph22, i) do
        A[i][j] := 1;
    od;
od;

Print("G_(2,2) degree set = ", VertexDegrees(graph22), "\n");

chiS := Concatenation(List(C51, x -> 1), List(C52, x -> -1));
Av := A * chiS;
Print("A_(2,2) * chi_S = ", Av, "\n");

cp := CharacteristicPolynomial(Rationals, Rationals, A);
Print("Characteristic polynomial of A_(2,2) = ", cp, "\n");
Print("Factored characteristic polynomial = ", Factors(cp), "\n");
Print("dim Ker(A_(2,2) + 5I) = ",
      Length(NullspaceMat(A + 5 * IdentityMat(Length(OmegaSet)))), "\n\n");

Print("Done.\n");

LogTo();
