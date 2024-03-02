loaded = readmatrix('./sample-data/Halfkernel.txt');
embedding = loaded(:, 1:2);
communities = loaded(:, 3);
variant = 'tsps';

[cpsIndices, cpsIndicesMetadata] = CommunitySeparability(embedding, communities, variant);

disp(cpsIndices);