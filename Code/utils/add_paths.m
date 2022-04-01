function add_paths()
    if ispc
        addpath(genpath('..\..\..\MissingData'));
        disp('Paths added')
    else
        disp('ERROR: Not running on PC. Fill in function')
    end
end
