function [auc, aupr] = computeAUCAUPR(labels, scores, positiveClass)
    [~,~,~,auc] = perfcurve(labels, scores, positiveClass);
    if auc < 0.5
        auc = 1 - auc;
        flippedScores = 2 * mean(scores) - scores;
        aupr = auprEvaluation(labels, flippedScores, positiveClass);
    else
        aupr = auprEvaluation(labels, scores, positiveClass);
    end
end

function aupr = auprEvaluation(labels, scores, positiveClass)
    [rec,prec,~,~] = perfcurve(labels, scores, positiveClass, 'xCrit', 'reca', 'yCrit', 'prec');
    % rec is the recall, prec is the precision.
    % the first value of rec (at recall 0) is NaN (by definition, PRECISION = TP / (FP + TP))
    % at recall 0 we have PRECISION = 0/(0 + 0) (we don't have neither TP nor FP)
    % if at the first point of recall (prec(2)) the precision is 1, the NaN is changed
    % for 1, in the contrary case (in the first point we have a FP), the NaN is changed for 0
    if prec(2) == 1
        prec(1) = 1;
    else
        prec(1) = 0;
    end
    aupr = trapz(rec,prec);
end
