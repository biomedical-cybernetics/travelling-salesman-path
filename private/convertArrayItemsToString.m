function convertedArray = convertArrayItemsToString(inputArray)
    if iscellstr(inputArray)
        convertedArray = inputArray;
        return;
    end

    if isnumeric(inputArray)
        convertedArray = arrayfun(@num2str, inputArray, 'UniformOutput', false);
        return;
    end

    error('The inputted array must be a numeric array or a cell array of strings.');
end
