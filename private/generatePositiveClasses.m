function positives = generatePositiveClasses(communities)
    % valid from R2019a
    [counter, list] = groupcounts(communities);

    % valid from R2014b
    %list = unique(communities);
    %counter = histcounts(communities);

    [~,idx] = max(counter);
    positives = list;
    positives(idx) = [];
end