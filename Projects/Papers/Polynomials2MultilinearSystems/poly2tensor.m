function A = poly2tensor(p)
    k = p.maxD;
    Am = p.Am;
    % Am = flip(Am,2);
    n = size(Am,1);
    A = sptensor([],[],(n+1)*ones(1,k+1));
    s = 1;
    for d=k:-1:0
        % Get matrix of dth degree polynomial terms
        AmD = Am(:,s:s+n^d-1);
        s = s + n^d;
        % Get tensor of dth degree polynomial terms
        Ad = polyD2tensor(AmD);
        
        if length(size(Ad)) > 1
            % Create command to copy the subtensors
            cmd = "A(";
            ct = k+1;
            for i=1:d+1; cmd = cmd + ":,"; ct = ct - 1; end
            while ct > 0; cmd = cmd + "n+1,"; ct = ct - 1; end
            cmd = char(cmd);
            cmd = cmd(1:length(cmd)-1);
            cmd = string(cmd);
            cmd = cmd + ") = Ad;";
            % disp(cmd);
            % Copy values of Ad into the subtensor structure of A
            eval(cmd);
        else
            ct = k;
            cmd = "A(i";
            while ct > 0; cmd = cmd + ",n+1"; ct = ct - 1; end
            cmd = cmd + ")=Ad(i);";
            for i=1:n; eval(cmd); end
            % disp(cmd);
        end
        % disp(tensor(A));
    end
end

function Ad = polyD2tensor(AmD)
    % get the number of variables and degree of these terms
    n = size(AmD,1);
    d = round(log(size(AmD,2)) / log(n));
    idxs = n * ones(1,d+1);
    if numel(idxs) > 1
        AmD = sptensor(AmD);
        AdN = reshape(AmD, idxs);
    else
        AdN = AmD;
    end
    Ad = sptensor([],[],(n+1)*ones(1,d+1));
    idxs = "1:n";
    for i=2:d+1; idxs = idxs + ",1:n"; end
    cmd = "Ad(" + idxs + ") = sptensor(AdN);";
    eval(cmd);

%     % create sparse tensor object
%     Ad = sptensor([],[],(n+1)*ones(1,d+1));
% 
%     % find polynomial terms from mvpoly
%     [r,c,v] = find(AmD);
% 
%     for i=1:length(r)
%         % map r and c of Am to the idxs of A
%         idxs        = zeros(1,d+1);
%         idxs(1)     = r(i);
%         cidxs       = kroneckerIndices(n,d,n^d-c(i));
%         idxs(2:end) = cidxs;
% 
%         % place value of polynomial term into the sparse tensor
%         Ad(idxs)     = v(i);
%     end
end