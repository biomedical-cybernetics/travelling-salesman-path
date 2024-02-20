loaded = readmatrix('./sample-data/reduced-halfkernel.txt');
embedding = loaded(:, 1:2);
communities = loaded(:, 3);
positives = generatePositiveClasses(communities);
variant = 'tsp';
enablePermutations = true;

%TODO: refactor input arguments
[cpsIndices, metadata] = CommunitySeparability(embedding, communities, positives, variant, enablePermutations);

disp(cpsIndices);
disp(metadata);