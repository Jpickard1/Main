function dp = data_path()
    current_path = pwd;
    len = length(current_path)
    current_path(len-9:end) = [];
    dp = string(current_path) + "Data/";
    disp(dp);
end
