function add_paths()
    if ispc
        addpath(genpath('..\..\..\Main'));
        disp('Paths added')
    else
        disp('ERROR: Not running on PC. Fill in function')
    end
end
