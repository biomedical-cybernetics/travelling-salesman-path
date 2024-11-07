runtimeSettings = struct();

runtimeSettings.rootPath = strcat(pwd, '/');
runtimeSettings.tempPath = strcat(runtimeSettings.rootPath, '.temp/');
createDirectory(runtimeSettings.tempPath);

concordeAppName = 'concorde';
if ispc
    cygwin32Path = 'C:\cygwin32\bin\';
    runtimeSettings.concordePath = strcat(cygwin32Path, 'concorde.exe');
elseif isunix
    runtimeSettings.concordePath = strcat(runtimeSettings.rootPath, 'bin/', 'concorde');
end

if ~isfile(runtimeSettings.concordePath)
    error('Concorde software was not found in the specified directory ''%s''. TSP computations require Concorde software. Please, read the ''REQUIREMENTS.md'' file for more details.', runtimeSettings.concordePath);
end