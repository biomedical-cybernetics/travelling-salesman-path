loaded = readmatrix('./sample-data/ParallelLines.txt');
embedding = loaded(:, 1:2);
communities = loaded(:, 3);
variant = 'ldps';
permutations = 1000;

[cpsPermutations, cpsPermutationsMetadata] = CommunitySeparability(embedding, communities, variant, 'permutations', permutations);

measures = fieldnames(cpsPermutations);
for ix=1:numel(measures)
    disp(measures{ix});
    disp(cpsPermutations.(measures{ix}));
end
