loaded = readmatrix('./sample-data/reduced-halfkernel.txt');
embedding = loaded(:, 1:2);
communities = loaded(:, 3);
variant = 'tsps';

[cpsIndices, cpsIndicesMetadata] = CommunitySeparability(embedding, communities, variant);

disp(cpsIndices);