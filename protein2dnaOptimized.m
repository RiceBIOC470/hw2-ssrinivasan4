%converting proteins to dna sequence

function dnaseq=protein2dnaOptimized(proteinseq)
codons=readtable('codons.csv');
am_acid_unique=unique(codons.AmAcid);
amino_dict1=containers.Map;
for k=1:length(am_acid_unique)
    all_index=find(contains(codons.AmAcid, am_acid_unique(k)));
    max_index=find(codons.x_1000(all_index)==max(codons.x_1000(all_index)));
    new_dict=containers.Map(codons.AmAcid(all_index(max_index)), codons.Codon(all_index(max_index)));
    amino_dict1=[amino_dict1; new_dict];
end
dnaseq=''
for i=1:3:length(proteinseq)/3
        dna1=amino_dict1(proteinseq(i:i+2));
        dnaseq=strcat(dnaseq,dna1);
end
end
