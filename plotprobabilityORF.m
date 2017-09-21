function mat = probabilityORF(N_ORF)
mat=[];
N_ORF=50;
for v=50:N_ORF:1000
    counter=0;
    for l = 1:1000
        dnaseq = randdnaseq(v);
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
        ORFsizes = firststop - start_pos + 3;
        
        if sum(ORFsizes > 3) > 0  && sum(ORFsizes > N_ORF) > 0
                    counter = counter + 1;
    
            counter=counter+1;
        end
    end
           
    end
    prob=[counter/1000];
    mat1=[v prob];
    mat=[mat; mat1];
end
plot(mat(:,1), mat(:,2));
title('Probability of Observing ORF > bp');
xlabel('Length of Sequence');
ylabel('Probability');
end