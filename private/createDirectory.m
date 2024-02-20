function createDirectory(dirPath)
    if (exist(dirPath,'dir')==7)
        return;
    end

    try
        mkdir(dirPath);
    catch
        error('The directory ''%s'' cannot be created.', dirPath);
    end
end