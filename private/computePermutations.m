function permutations = computePermutations(variant, previousResults, metadata, communities, positives, totalPermutations)
    totalPairwiseCombinations = numel(metadata);

    aucValues = zeros(totalPermutations, totalPairwiseCombinations);
    auprValues = zeros(totalPermutations, totalPairwiseCombinations);
    mccValues = zeros(totalPermutations, totalPairwiseCombinations);

    for ix=1:totalPairwiseCombinations
        meta = metadata(ix);

        communitiesGroupA = communities(ismember(communities, meta.communityNameGroupA));
        communitiesGroupB = communities(ismember(communities, meta.communityNameGroupB));
        communitiesMembership = [communitiesGroupA;communitiesGroupB];
        permutedCommunitiesMembership = randomizeCommunities(communitiesMembership, totalPermutations);
        currentPositiveClass = extractPositiveClass(communitiesMembership, positives);
        scores = meta.scores;

        [aucValues(:, ix), auprValues(:, ix), mccValues(:, ix)] = computeCpsMeasures(permutedCommunitiesMembership, scores, currentPositiveClass, totalPermutations);
    end

    % Permutations
    permutations = struct();

    correctedAUCValues = mean(aucValues, 2) ./ (1 + std(aucValues, [], 2));
    aucOriginalVariantValue = previousResults.auc;
    permutations.auc = setPermutationResults(correctedAUCValues, aucOriginalVariantValue, totalPermutations);

    correctedAUPRValues = mean(auprValues, 2) ./ (1 + std(auprValues, [], 2));
    auprOriginalVariantValue = previousResults.aupr;
    permutations.aupr = setPermutationResults(correctedAUPRValues, auprOriginalVariantValue, totalPermutations);

    correctedMCCValues = mean(mccValues, 2) ./ (1 + std(mccValues, [], 2));
    mccOriginalVariantValue = previousResults.mcc;
    permutations.mcc = setPermutationResults(correctedMCCValues, mccOriginalVariantValue, totalPermutations);
end

function [aucValues, auprValues, mccValues] = computeCpsMeasures(permutedCommunitiesMembership, scores, currentPositiveClass, totalPermutations)
    aucValues = zeros(totalPermutations, 1);
    auprValues = zeros(totalPermutations, 1);
    mccValues = zeros(totalPermutations, 1);

    parfor jx=1:totalPermutations
        permutedCommunities = permutedCommunitiesMembership{jx};

        % AUC & AUPR
        [aucValues(jx), auprValues(jx)] = computeAUCAUPR(permutedCommunities, scores, currentPositiveClass);

        % MCC
        mccValues(jx) = computeMCC(permutedCommunities, scores, currentPositiveClass);
    end
end

function currentPositiveClass = extractPositiveClass(communitiesMembership, positives)
    % TODO: Refactor, extract in a separate method.
    % think about swapping the positive and negative classes.
    % selecting the possitive class
    currentPositiveClass = NaN;
    for o=1:numel(positives)
        if any(ismember(communitiesMembership,positives{o}))
            currentPositiveClass = positives{o};
            break;
        end
    end

    if ~isnan(currentPositiveClass)
        return;
    end

    error('positive class cannot be extracted from communities membership');
end

function randomized = randomizeCommunities(communities, totalPermutations)
    randomized = cell(1, totalPermutations);
    totalCommunities = numel(communities);
    for ix=1:totalPermutations
        rng(ix, 'twister');
        positions = randperm(totalCommunities);
        randomized{ix} = communities(positions);
    end
end

function permutations = setPermutationResults(permutedValues, originalValue, totalPermutations)
    permutations = struct();
    permutations.originalValue = originalValue;
    permutations.permutations = permutedValues;
    permutations.pvalue = (sum(permutedValues >= originalValue) + 1) / (totalPermutations + 1);
    permutations.mean = mean(permutedValues);
    permutations.max = max(permutedValues, [], 'all');
    permutations.min = min(permutedValues, [], 'all');
    permutations.standardDeviation = std(permutedValues);
    permutations.standardError = std(permutedValues) / sqrt(totalPermutations);
end
