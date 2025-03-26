function install_logger_addon()
    %% Check if installed
    toolboxes = matlab.addons.toolbox.installedToolboxes;
    if ~isempty(toolboxes)
        tbnames = {toolboxes.Name}';
    else
        tbnames={};
    end

    if ismember("Advanced Logger for MATLAB",tbnames)
        disp('Add-on already installed!')
        return
    end
    
    %% Download add-on

    websave('mlogger.mltbx','https://nl.mathworks.com/matlabcentral/mlc-downloads/downloads/fd9733c5-082a-4325-a5e5-e7490cdb8fb1/43beb4a8-75fe-4275-be01-232a19e1fb9e/packages/mltbx')

    %% Install and test add-on
    %IGNORE:: Warning: Package directories not allowed in MATLAB path
    try
        matlab.addons.install('mlogger.mltbx',true)
        disp('Ignore the warning: Package directories not allowed in MATLAB path')
        logObj = mlog.Logger('testLog');
        warning(logObj,'test warning')
        disp('Add-on installed correctly!')
        delete(logObj)
        return
    catch err
        error('Add-on not installed correctly, check download link or manually install from the add-on installer GUI. \nMatlab error message:%s', err.message)
    end
end