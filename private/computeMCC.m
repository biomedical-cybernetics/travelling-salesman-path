function mcc = computeMCC(labels, scores, positiveClass)
    % total amount of positive labels
    totalPositive = sum(strcmp(labels, positiveClass));
    % total amount of negative labels
    totalNegative = sum(~strcmp(labels, positiveClass));

    % identifying the negative class
    negativeClass = unique(labels(~strcmp(labels, positiveClass)));

    % sort the scores and obtained the sorted indices
    [~, idxs] = sort(scores);

    % sort the original labels according to the sorted scores
    trueLabels = labels(idxs);

    for in=1:2
        % since the position of the groups is unknown
        % we take as comparison a perfect segregation in both ways
        % positive group in the left side, and negative group in the right side,
        % and viseversa
        switch in
            case 1
                predictedLabels = [repmat({positiveClass},totalPositive,1);repmat(negativeClass,totalNegative,1)];
            case 2
                predictedLabels = [repmat(negativeClass,totalNegative,1);repmat({positiveClass},totalPositive,1)];
        end

        % clasifiers
        TP = sum(ismember(predictedLabels, positiveClass) & ismember(trueLabels, positiveClass));
        TN = sum(ismember(predictedLabels, negativeClass) & ismember(trueLabels, negativeClass));
        FP = sum(ismember(predictedLabels, positiveClass) & ismember(trueLabels, negativeClass));
        FN = sum(ismember(predictedLabels, negativeClass) & ismember(trueLabels, positiveClass));

        if ((TP == 0 && FP == 0) || (TN == 0 && FN == 0))
            mccVal{in} = 0;
        else
            % compute the Matthews Correlation Coefficient (MCC)
            mccVal{in} = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
        end
    end
    % select the best MCC side
    mcc = max([mccVal{:}]);
end
