function randomSeq = randdnaseq(N)
letters=["A", "T", "C", "G"];
nums = randi(numel(letters),[1 N]);
randomSeq = strjoin(letters(nums), "");
randomSeq=char(randomSeq);
% returns a random dna sequence of length N
