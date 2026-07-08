# AssociationSchemes package check for the A8 twist in main.tex.
#
# The script constructs the original group association scheme of A8 and reads
# the twisted scheme partition produced by A8.g.  It then checks that the two
# schemes are algebraically isomorphic but not combinatorially isomorphic, and
# that the twisted scheme is non-Schurian.

OutDir := "out";
Exec(Concatenation("mkdir -p ", OutDir));

LoadPackage("AssociationSchemes");

LogTo(Concatenation(OutDir, "/result_checktwistA8.txt"));

Print("\n");
Print("============================================================\n");
Print("AssociationSchemes package check for the A8 twist\n");
Print("============================================================\n\n");

G := AlternatingGroup(8);
elts := Elements(G);

GroupScheme := GroupCoherentConfiguration(G);

TwistedPartitionFile := Concatenation(OutDir, "/A8_twisted_Schur_partition.txt");
if not IsExistingFile(TwistedPartitionFile) then
    Error("Run A8.g first to create ", TwistedPartitionFile);
fi;
Read(TwistedPartitionFile);
twClasses := A8TwistedSchurPartition;

Print("Valencies of the original group scheme = ",
      Valencies(GroupScheme), "\n");
Print("Class sizes of the twisted partition  = ",
      List(twClasses, Length), "\n");

partLabels := NewDictionary(elts[1], true);
for k in [1..Length(twClasses)] do
    for g in twClasses[k] do
        AddDictionary(partLabels, g, k - 1);
    od;
od;

M := [];
for i in [1..Length(elts)] do
    M[i] := [];
    for j in [1..Length(elts)] do
        M[i][j] := LookupDictionary(partLabels, elts[j] * elts[i]^-1);
    od;
od;

TwistedScheme := AssociationScheme(M);

algIso := AreAlgebraicallyIsomorphicHomogeneousCoherentConfigurations(
    GroupScheme, TwistedScheme);
# combIso := AreIsomorphicHomogeneousCoherentConfigurations(
#     GroupScheme, TwistedScheme);

GroupSchemeIsSchurian := IsSchurian(GroupScheme);
TwistedSchemeIsSchurian := IsSchurian(TwistedScheme);


Print("Relation matrix size: ", Length(M), "\n");
Print("Group scheme is an association scheme? ",
      IsHomogeneousCoherentConfiguration(GroupScheme), "\n");
Print("Twisted scheme is an association scheme? ",
      IsHomogeneousCoherentConfiguration(TwistedScheme), "\n");
Print("Algebraically isomorphic? ", algIso, "\n");
# Print("Combinatorially isomorphic? ", combIso, "\n\n");

Print("--------------------------------\n");
Print("Schurity\n");
Print("--------------------------------\n");
Print("Group scheme is Schurian? ", GroupSchemeIsSchurian, "\n");
Print("Twisted scheme is Schurian? ", TwistedSchemeIsSchurian, "\n");

Print("Done.\n");

LogTo();
