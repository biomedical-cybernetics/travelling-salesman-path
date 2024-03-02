function [measures, metadata] = CommunitySeparability(embedding, communities, variant, varargin)
    p = inputParser;
    addRequired(p, 'embedding',...
        @(x) assert(~isempty(x) && isnumeric(x) && size(x, 2) >= 2, 'It must be a matrix of at least two dimensions.'));
    addRequired(p, 'communities',...
        @(x) assert(~isempty(x) && (iscellstr(x) || isnumeric(x)), 'It must be a numeric array or a cell array of strings.'));
    addRequired(p, 'variant',...
        @(x) assert(any(ismember(x, {'tsps', 'cps', 'ldps'})), 'It must be one of: tsps, cps, or ldps.'));
    addParameter(p, 'positives', [],...
        @(x) assert(~isempty(x) && (iscellstr(x) || isnumeric(x)), 'It must be a numeric array or a cell array of strings.'));
    addParameter(p, 'permutations', 0,...
        @(x) assert(isnumeric(x) && x >= 0, 'It must be numeric and equal or higher than zero.'));
    parse(p, embedding, communities, variant, varargin{:});
    inputs = p.Results;

    embedding = inputs.embedding;
    communities = inputs.communities;
    variant = inputs.variant;
    positives = inputs.positives;
    permutations = inputs.permutations;

    % load runtime settings
    settings

    % sanity check
    communities = convertArrayItemsToString(communities);
    if isempty(positives)
        positives = generatePositiveClasses(communities);
    end
    positives = convertArrayItemsToString(positives);

    uniqueCommunities = unique(communities);
    totalCommunities = length(uniqueCommunities);

    % grouping data according to communities
    communitiesClustered = cell(1, totalCommunities);
    dataClustered = cell(1, totalCommunities);
    for k=1:totalCommunities
        idx = ismember(communities, uniqueCommunities{k});
        communitiesClustered{k} = communities(idx);
        dataClustered{k} = embedding(idx,:);
    end

    pairwiseGroupCombinations = nchoosek(1:totalCommunities, 2);
    totalPairwiseCombinations = size(pairwiseGroupCombinations, 1);

    aucValues = cell(1, totalPairwiseCombinations);
    auprValues = cell(1, totalPairwiseCombinations);
    mccValues = cell(1, totalPairwiseCombinations);
    %mRocValues = cell(1, totalPairwiseCombinations);

    metadata = struct();

    for l=1:totalPairwiseCombinations
        idxGroupA = pairwiseGroupCombinations(l, 1);
        dataGroupA = dataClustered{idxGroupA};
        communitiesGroupA = communitiesClustered{idxGroupA};
        communityNameGroupA = uniqueCommunities{idxGroupA};
        metadata(l).communityNameGroupA = communityNameGroupA;
        metadata(l).dataGroupA = dataGroupA;

        idxGroupB = pairwiseGroupCombinations(l, 2);
        dataGroupB = dataClustered{idxGroupB};
        communitiesGroupB = communitiesClustered{idxGroupB};
        communityNameGroupB = uniqueCommunities{idxGroupB};
        metadata(l).communityNameGroupB = communityNameGroupB;
        metadata(l).dataGroupB = dataGroupB;

        scores = [];
        switch variant
        case {'cps'}
            projectedPoints = centroidBasedProjection(dataGroupA, dataGroupB, 'median');
            if ~isempty(projectedPoints)
                scores = convertPointsToOneDimension(projectedPoints);
            end
        case {'ldps'}
            pairwiseData = [dataGroupA;dataGroupB];
            pairwiseCommunities = [communitiesGroupA;communitiesGroupB];
            projectedPoints = ldaBasedProjection(pairwiseData, pairwiseCommunities);
            if ~isempty(projectedPoints)
                scores = convertPointsToOneDimension(projectedPoints);
            end
        case {'tsps'}
            pairwiseData = [dataGroupA;dataGroupB];
            metadata(l).pairwiseData = pairwiseData;

            pairwiseCommunities = [communitiesGroupA;communitiesGroupB];
            metadata(l).pairwiseCommunities = pairwiseCommunities;

            bestTour = tspBasedProjection(pairwiseData, runtimeSettings);
            metadata(l).bestTour = bestTour;

            if ~isempty(bestTour)
                [scores, sourceNodes, targetNodes, edgeWeights, cutStartNode, cutEndNode] = convertTourToOneDimension(bestTour, pairwiseData, pairwiseCommunities);
                metadata(l).sourceNodes = sourceNodes;
                metadata(l).targetNodes = targetNodes;
                metadata(l).edgeWeights = edgeWeights;
                metadata(l).cutStartNode = cutStartNode;
                metadata(l).cutEndNode = cutEndNode;
            end
        otherwise
            error('Community separability variant ''%s'' not supported', variant);
        end

        metadata(l).scores = scores;

        if isempty(scores)
            aucValues{l} = 0;
            auprValues{l} = 0;
            mccValues{l} = 0;
            %mRocValues{l} = 0;
            continue;
        end

        % sample membership
        communitiesMembership = [communitiesGroupA;communitiesGroupB];

        % TODO: Refactor, extract in a separate method.
        % think about swapping the positive and negative classes.
        % selecting the possitive class
        for o=1:length(positives)
            if any(ismember(communitiesMembership,positives{o}))
                currentPositiveClass = positives{o};
                break;
            end
        end

        % AUC & AUPR
        [aucValues{l}, auprValues{l}] = computeAUCAUPR(communitiesMembership,scores,currentPositiveClass);

        % MCC
        mccValues{l} = computeMCC(communitiesMembership,scores,currentPositiveClass);

        % mAUC-ROC
        %mRocValues{l} = computeMROC(communitiesMembership, scores, currentPositiveClass);
    end

    % compile all values from the different groups' combinations
    allAUCROCvalues = [aucValues{:}];
    allAUCPRvalues = [auprValues{:}];
    allMCCvalues = [mccValues{:}];
    %allmRocValues = [mRocValues{:}];

    % Corrected
    correctedAUC = mean(allAUCROCvalues) / (1 + std(allAUCROCvalues));
    correctedAUPR = mean(allAUCPRvalues) / (1 + std(allAUCPRvalues));
    correctedMCC = mean(allMCCvalues) / (1 + std(allMCCvalues));
    %correctedmRoc = mean(allmRocValues) / (1 + std(allmRocValues));

    measures.auc = correctedAUC;
    measures.aupr = correctedAUPR;
    measures.mcc = correctedMCC;
    %measures.mroc = correctedmRoc;

    if (permutations > 1)
        measures = computePermutations(variant, measures, metadata, communities, positives, permutations);
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sub-functions
function bestroute = tspBasedProjection(pairwiseData, runtimeSettings)
    if ~isfield(runtimeSettings, 'concordePath')
        error('Concorde path not specified in runtime settings.');
    end

    % sanity check to avoid overriding TSP input and output files when executing Concorde.
    while 1
        shortUUID = regexp(char(java.util.UUID.randomUUID.toString),'\w+(?:-)\w+$','match');
        fileName = strrep(shortUUID{:}, '-', '');

        filePath = strcat(runtimeSettings.tempPath, fileName);
        fileTSP = strcat(filePath, '.tsp');
        fileSol = strcat(filePath, '.sol');

        if (~isfile(fileTSP) && ~isfile(fileSol))
            break;
        else
            warning('file path ''%s'' already exists. Retrying to create a new file.', filePath);
            pause(0.25);
        end
    end

    totalNodes = size(pairwiseData, 1);

    % sanity check to avoid "rounding-up problem" when using EUC_2D
    % https://stackoverflow.com/questions/27304721/how-to-improve-the-quality-of-the-concorde-tsp-solver-am-i-misusing-it

    absMean = abs(median(pairwiseData, 'all'));
    digits = 3;
    invertedMagnitude = 10^(2-digits+floor(log10(absMean)));
    offset = 1;
    while 1
        scalingFactor = 10^(abs(log10(invertedMagnitude))+offset);
        scaledEmbedding = pairwiseData*scalingFactor;
        maxValue = max(scaledEmbedding, [], 'all');
        maxDistance = max(pdist2(scaledEmbedding, scaledEmbedding, 'euclidean'), [], 'all');
        if (maxValue < intmax && maxDistance < 32768)
            break;
        end
        offset = offset - 1;
    end

    % prepare TSP file
    fileID = fopen(fileTSP, 'w');
    fprintf(fileID, 'NAME : %s\n', 'CPS Concorde');
    fprintf(fileID, 'COMMENT : Scaling factor %d\n', scalingFactor);
    fprintf(fileID, 'TYPE : %s\n', 'TSP');
    fprintf(fileID, 'DIMENSION : %d\n', totalNodes);
    fprintf(fileID, 'EDGE_WEIGHT_TYPE : %s\n', 'EUC_2D');
    fprintf(fileID, 'NODE_COORD_SECTION\n');
    for ix=1:totalNodes
        fprintf(fileID, '%d %f %f\n', ix, scaledEmbedding(ix, 1), scaledEmbedding(ix, 2));
    end
    fprintf(fileID, 'EOF');
    fclose(fileID);

    command = [runtimeSettings.concordePath ' -s 40 -x -o ' fileSol ' ' fileTSP];
    [status, cmdOut] = system(command);
    if (status ~= 0 && status ~= 255)
        error('Error executing Concorde: %s', cmdOut);
    end

    try
        fileID = fopen(fileSol, 'r');
        bestroute = fscanf(fileID, '%d');
        fclose(fileID);
    catch
        warning('Concorde solution ''%s'' not found. Concorde output: %s. Returning empty tour!', fileSol, cmdOut);
        bestroute = [];
        return;
    end

    % the first values outputted by concorde is the number of nodes
    % thus, we can remove it
    bestroute(1) = [];
    % concorde will start counting the nodes from 0; thus, we need
    % to increase the nodes numering by 1
    bestroute = bestroute + 1;
    % reshaping
    bestroute = bestroute';

    % clean up files
    if isfile(fileTSP)
        delete(fileTSP);
    end
    if isfile(fileSol)
        delete(fileSol);
    end
    % sometimes, the solution file is created in the root directory
    % and not cleaned up after the execution
    tempFileSol = strcat(runtimeSettings.rootPath, fileName, '.sol');
    if isfile(tempFileSol)
        delete(tempFileSol);
    end
end

function [scores, inputNodes, targetNodes, weights, startNode, endNode] = convertTourToOneDimension(bestTour, pairwiseData, pairwiseCommunities)
    inputNodes = bestTour';
    targetNodes = [inputNodes(2:end);inputNodes(1)];

    weights = zeros(numel(inputNodes),1);
    for ix=1:numel(inputNodes)
        weights(ix) = norm(pairwiseData(inputNodes(ix), :) - pairwiseData(targetNodes(ix), :));
    end

    % weighted Graph by Euclidean distances
    H = graph(inputNodes, targetNodes, weights);

    startNode = NaN;
    endNode = NaN;
    [S, I] = sort(H.Edges.Weight, 'desc');
    for ix=1:numel(S)
        E = H.Edges(I(ix), :);
        L1 = pairwiseCommunities(E.EndNodes(1));
        L2 = pairwiseCommunities(E.EndNodes(2));
        if ~strcmp(L1, L2)
            startNode = E.EndNodes(1);
            endNode = E.EndNodes(2);
            break;
        end
    end
    if (isnan(startNode) || isnan(endNode))
        error('Cannot find best tour cut!');
    end
    H = rmedge(H, startNode, endNode);

    sPath = shortestpath(H, startNode, endNode);
    scores = zeros(numel(sPath), 1);
    for ix=1:numel(sPath)
        node = sPath(ix);
        [~,D] = shortestpath(H,node,startNode);
        scores(node) = D;
        if (node == startNode)
            continue;
        end
    end
end

function projection = centroidBasedProjection(dataGroupA, dataGroupB, centerFormula)
    if (~strcmp(centerFormula,'mean') && ~strcmp(centerFormula,'median') && ~strcmp(centerFormula,'mode'))
        warning('Invalid center formula ''%s''; median will be applied instead!', centerFormula);
        centerFormula = 'median';
    end

    switch centerFormula
        case 'mean'
            centroidA = mean(dataGroupA,1);
            centroidB = mean(dataGroupB,1);
        case 'median'
            centroidA = median(dataGroupA,1);
            centroidB = median(dataGroupB,1);
        case 'mode'
            centroidA = modeDistribution(dataGroupA);
            centroidB = modeDistribution(dataGroupB);
        otherwise
            error('You must select either mean, median, or mode to calculate the centroids');
    end
    if centroidA == centroidB
        warning('Pairwise groups have the same centroid: no line can be traced between them. Returning empty projection!');
        projection = [];
        return;
    end

    centroidsLine = createLineBetweenCentroids(centroidA,centroidB);
    pairwiseData = [dataGroupA; dataGroupB];

    projection = zeros(size(pairwiseData, 1), size(pairwiseData, 2));
    for ox=1:size(pairwiseData, 1)
        projection(ox, :) = projectPointOnLine(pairwiseData(ox, :),centroidsLine);
    end
end

function projection = ldaBasedProjection(pairwiseData, pairwiseSamples)
    try
        Mdl = fitcdiscr(pairwiseData, pairwiseSamples, 'DiscrimType', 'linear');
    catch ME
        warning(ME.identifier, 'Cannot create LDA projection: %s. Returning empty projection!', ME.message)
        projection = [];
        return;
    end

    [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma);
    lambda = diag(LAMBDA);
    [~, SortOrder] = sort(lambda, 'descend');
    W = W(:, SortOrder);
    mu = mean(pairwiseData);
    % projecting data points onto the first discriminant axis
    centered = bsxfun(@minus, pairwiseData, mu);
    projection = centered * W(:,1) * transpose(W(:,1));
    projection = bsxfun(@plus, projection, mu);
end

function modeDist = modeDistribution(distance)
    modeDist = [];
    for ix=1:size(distance,2)
        [f,xi] = ksdensity(distance(:,ix));
        [~,ind] = max(f);
        modeDist = horzcat(modeDist, xi(ind));
    end
end

function centroidsLine = createLineBetweenCentroids(point1, point2)
    centroidsLine = [point1;point2];
end

function projectedPoint = projectPointOnLine(point, line)
    A = line(1,:);
    B = line(2,:);

    AP = point - A;
    AB = B - A;

    projectedPoint = A + dot(AP,AB)/dot(AB, AB) * AB;
end

function V = convertPointsToOneDimension(points)
    % select the point 0 (min value of an axis [where the values are different from one another])
    for ix=1:size(points,2)
        if length(unique(points(:,ix))) ~= 1
            startPoint = points(points(:,ix) == min(points(:,ix)),:);
            break;
        end
    end

    V = 0;
    for ix=1:size(points,2)
        V = V + (points(:,ix)-min(startPoint(:,ix))).^2;
    end

    V = sqrt(V);
end
