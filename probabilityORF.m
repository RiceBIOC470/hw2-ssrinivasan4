function [prob] = probabilityORF(N, N_ORF)
letters=["A", "T", "C", "G"];
nums = randi(numel(letters),[1 N]);
dnaseq = strjoin(letters(nums), "");
good_length= '';
orf_min=[];
start_pos=strfind(upper(dnaseq), "ATG");
stop_pos=[strfind(upper(dnaseq), "TAA") strfind(upper(dnaseq), "TGA") strfind(upper(dnaseq), "TAG")];
for i = 1:length(start_pos)
    ORFlength=stop_pos(stop_pos > start_pos(i)) - start_pos(i);
    for j = 1:length(ORFlength)
        if mod(ORFlength(j), 3) == 0
            good_length = [good_length, ORFlength(j)];
            good_length= min(good_length)+3; %finding the minimum per start codon - this is the first stop codon it hits
        end
    end
    orf_min=[orf_min, good_length];
end
numer=sum(orf_min > N_ORF);
denom=length(orf_min);
prob=[numer/denom];
disp(['The probability of seeing ORFs longer than ' int2str(N_ORF) ' base pairs for a total sequence length of '  int2str(N) ' is ' num2str(prob)]);
end