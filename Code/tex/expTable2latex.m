function expTable2latex(T, filename)
% ----------------------------------------------------------------------- %
% Function table2latex(T, filename) converts a given MATLAB(R) table into %
% a plain .tex file with LaTeX formatting.                                %
%                                                                         %
%   Input parameters:                                                     %
%       - T:        MATLAB(R) table. The table should contain only the    %
%                   following data types: numeric, boolean, char or string.
%                   Avoid including structs or cells.                     %
%       - filename: (Optional) Output path, including the name of the file.
%                   If not specified, the table will be stored in a       %
%                   './table.tex' file.                                   %
% ----------------------------------------------------------------------- %
%   Version: 1.1                                                          %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    09/10/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
    % Error detection and default parameters
    if nargin < 2
        filename = 'table.txt';
        fprintf('Output path is not defined. The table will be written in %s.\n', filename); 
    elseif ~ischar(filename)
        error('The output file name must be a string.');
    end
    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(T), error('Input must be a table.'); end
    
    % Parameters
    n_col = size(T,2);
    col_spec = ['|'];
    for c = 1:n_col, col_spec = [col_spec 'c|']; end
    % col_names = strjoin(T.Properties.VariableNames, ' & ');
    row_names = T.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec];
        col_names = ['& ' col_names];
    end
    
    % Writing header
    fileID = fopen(filename, 'w');
    fprintf(fileID, '\\begin{tabular}{%s}\n', col_spec);
    fprintf(fileID, '\\hline \n');
    col_names = T.Properties.VariableNames;
    str = col_names{1};
    for i=2:length(col_names)
        str = str + " & " + col_names(i);
    end
    fprintf(fileID, str);
    fprintf(fileID, '\\\\ \n');
    fprintf(fileID, '\\hline \n');
    fprintf(fileID, '\\hline \n');    
    % Writing the data
    try
        for row = 1:size(T,1)
            temp = cell(1,n_col);
            for col = 1:n_col
                S = T{row,col};
                temp{1,col} = char(string(S));
            end
            if ~isempty(row_names)
                temp = [row_names{row}, temp];
            end
            fprintf(fileID, '%s \\\\ \n', strjoin(temp, ' & '));
            clear temp;
            fprintf(fileID, '\\hline \n');
        end
    catch
        error('Unknown error. Make sure that table only contains chars, strings or numeric values.');
    end
    
    % Closing the file
    fprintf(fileID, '\\end{tabular}');
    fclose(fileID);
end