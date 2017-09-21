function [ORFlength, start_pos, stop_pos] = findORF(dnaseq)
%Function to find the length of the longest open reading frame of a
%sequences called dnaseq
start_pos=strfind(upper(dnaseq), 'ATG');
stop_pos=[strfind(upper(dnaseq), 'TAA') strfind(upper(dnaseq), 'TGA') strfind(upper(dnaseq), 'TAG')];

if ~isempty(start_pos) && ~isempty(stop_pos)
    firststop = zeros(1, length(start_pos));
    
for ii = 1:length(start_pos)
    ORFlength = stop_pos - start_pos(ii);
    good_length = 1e8;
    good_index = 0;
    for jj = 1:length(ORFlength)
        if ORFlength(jj) > 0
                mod(ORFlength(jj),3) == 0;
                ORFlength(jj) < good_length;
           good_length = ORFlength(jj);
           good_index = jj;
        end
    end
    if good_index > 0
        firststop(ii) = stop_pos(good_index);
    else
        firststop(ii) = start_pos(ii);
        end
    end
    ORFsize = firststop - start_pos + 3;
    [ORFmaxlength, ind_max] = max(ORFsize);
if ORFmaxlength > 0
    disp(['Length of the longest ORF is ' int2str(ORFmaxlength) ' base pairs long. Start codon at ' int2str(start_pos(ind_max)) '. Stop codon at ' int2str(firststop(ind_max))]);
else
    disp(['No open reading frame found']);
end
end
