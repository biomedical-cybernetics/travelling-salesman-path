function mroc = computeMROC(labels, scores, positiveClass)
    % binary labels (positive class as TRUE and negative class as FALSE)
    binaryLabels = strcmp(labels, positiveClass);

    measures = prediction_evaluation(scores, binaryLabels);
    mroc = measures.auc_mroc;

    if (mroc >= 0.5)
        return;
    end

    flippedScores = 2 * mean(scores) - scores;
    measures = prediction_evaluation(flippedScores, binaryLabels);
    mroc = measures.auc_mroc;
end
