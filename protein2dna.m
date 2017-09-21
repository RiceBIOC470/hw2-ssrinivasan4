%converting proteins to dna sequence

function dnaseq=protein2dna(proteinseq)
codons=readtable('codons.csv');
am_acid_unique=unique(codons.AmAcid);
amino_dict=containers.Map;
for k=1:length(am_acid_unique)
    all_index=find(contains(codons.AmAcid, am_acid_unique(k)));
    first_index=all_index(1);
    new_dict=containers.Map(codons.AmAcid(first_index), codons.Codon(first_index));
    amino_dict=[amino_dict; new_dict];
end
dnaseq=''
for i=1:3:length(proteinseq)/3
        dna1=amino_dict(proteinseq(i:i+2));
        dnaseq=strcat(dnaseq,dna1);
end
end