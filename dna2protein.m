function proteinseq = dna2protein(dnaseq,frame)
codons=readtable('codons.csv');
codon_dict=containers.Map(codons.Codon, codons.AmAcid);
proteinseq=''
if ismember(frame, 1:3)
    for i=frame:3:length(dnaseq)-2
        prot1=codon_dict(dnaseq(i:i+2));
        proteinseq=strcat(proteinseq,prot1);
    end
else
    disp(['Please enter valid frame of 1, 2 or 3']);
end
end

% Input a dna sequence and a reading frame and returns the corresponding
% protein sequence. 


    
    
    