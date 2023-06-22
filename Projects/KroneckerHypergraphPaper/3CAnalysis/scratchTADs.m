%% Scratch TADs
%
%   Idea: we enforce a Kronecker structure to the polymer by supposing that
%   there are 5 compartments, and given the fracel globule premis we model
%   internal and between compartment interactions


path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\KroneckerHypergraphPaper\3CAnalysis\HipHop2HG\';
filename = "adjList_HIGH_1_0.01.txt";
filepath = path2data + filename;

E = readmatrix(filepath);

B = zeros(5,5,5);
I = zeros(1000, 1000, 1000);
for e=1:size(E,1)
    disp(e);
    edge = E(e,:);
    r = range(edge);
    % Internal edge
    if r <= 1000
        edge = mod(edge - 1, 1000) + 1;
        pedge = perms(edge);
        for p=1:6
            I(pedge(p,1), pedge(p,2), pedge(p,3)) = I(pedge(p,1), pedge(p,2), pedge(p,3)) + 1;
        end
    % Between compartment edge
    else
        edge = floor(edge / 1000) + 1;
        pedge = perms(edge);
        for p=1:6
            B(pedge(p,1), pedge(p,2), pedge(p,3)) = B(pedge(p,1), pedge(p,2), pedge(p,3)) + 1;
        end
    end
end

%% Plot colored polymer
filePath = "C:\Users\picka\Documents\my_projects\DBTM\Main\Code\reproductions\Hip-Hop\HiP-HoP_Pax6_FullConformations\OFF\conf.1.DNA";
T = readLAMMPSoutput(filePath);
cords = getCords(T);
figure;
for i = 1:size(B,1)
    vxc = (i-1) * 1000 + [1:1000];
    plot3(cords(vxc,1), cords(vxc,2), cords(vxc,3), '.-'); hold on;
end
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Polymer');


%% Analysis of above for 15x3

SIMPARMS = ["OFF","ON","HIGH"];
nsamps = 200;
BC = cell(nsamps,3);
path2data = 'C:\Users\picka\Documents\my_projects\DBTM\Main\Projects\KroneckerHypergraphPaper\3CAnalysis\HipHop2HG\';
for s=1:3
    for i=1:nsamps
        
        filename = "adjList_" + SIMPARMS(s) + "_" + string(i) + "_0.01.txt";
        filepath = path2data + filename; disp(filepath);
        
        E = readmatrix(filepath);
        
        B = zeros(5,5,5);
        for e=1:size(E,1)
            % disp(e);
            edge = E(e,:);
            r = range(ceil(edge/1000)*1000);
            % Internal edge
            if r >= 1000
                % disp(edge);
                edge = floor((edge - 1) / 1000) + 1;
                % disp(edge);
                pedge = perms(edge);
                pedge = unique(pedge,'rows');
                for p=1:size(pedge,1)
                    B(pedge(p,1), pedge(p,2), pedge(p,3)) = B(pedge(p,1), pedge(p,2), pedge(p,3)) + 1;
                end
            end
        end
        BC{i,s} = B;
    end
end

%% Compare tensor distances
D = zeros(nsamps*3,nsamps*3);
for i=1:nsamps*3
    ii = 1;
    if i > nsamps*2
        ii = 3;
    elseif i > nsamps
        ii = 2;
    end
    T1 = BC{mod(i-1,nsamps) + 1, ii};
    for j=1:nsamps*3
        jj = 1;
        if j > nsamps*2
            jj = 3;
        elseif j > nsamps
            jj = 2;
        end
        disp(string(ii) + "," + string(mod(i-1,nsamps) + 1) + " -- " + string(jj) + "," +  string(mod(j-1,nsamps) + 1));
        T2 = BC{mod(j-1,nsamps) + 1, jj};

        Td = abs(T1 - T2);
        D(i,j) = sum(Td,'all');

    end
end
%%
figure; imagesc(D);

Ds = zeros(3,3);
Ds(1,1) = sum(D(1:nsamps,            1:nsamps         ), 'all');
Ds(1,2) = sum(D(1:nsamps,            nsamps+1:2*nsamps), 'all');
Ds(1,3) = sum(D(1:nsamps,            2*nsamps+1:3*nsamps), 'all');
Ds(2,1) = sum(D(nsamps+1:2*nsamps,   1:nsamps         ), 'all');
Ds(2,2) = sum(D(nsamps+1:2*nsamps,   nsamps+1:2*nsamps), 'all');
Ds(2,3) = sum(D(nsamps+1:2*nsamps,   2*nsamps+1:3*nsamps), 'all');
Ds(3,1) = sum(D(2*nsamps+1:3*nsamps, 1:nsamps         ), 'all');
Ds(3,2) = sum(D(2*nsamps+1:3*nsamps, nsamps+1:2*nsamps), 'all');
Ds(3,3) = sum(D(2*nsamps+1:3*nsamps, 2*nsamps+1:3*nsamps), 'all');
disp(Ds)

figure; imagesc(Ds)

%% Map each tensor to a matrix (i.e. cliqe expand) and perform similar analysis
BMC = cell(size(BC));
for s=1:3
    for i=1:nsamps
        BMC{i,s} = tensorprod(BC{i,s}, ones(5,1), 1, 1);
    end
end
D = zeros(nsamps*3,nsamps*3);
for i=1:nsamps*3
    ii = 1;
    if i > nsamps*2
        ii = 3;
    elseif i > nsamps
        ii = 2;
    end
    T1 = BMC{mod(i-1,nsamps) + 1, ii};
    for j=1:nsamps*3
        jj = 1;
        if j > nsamps*2
            jj = 3;
        elseif j > nsamps
            jj = 2;
        end
        disp(string(ii) + "," + string(mod(i-1,nsamps) + 1) + " -- " + string(jj) + "," +  string(mod(j-1,nsamps) + 1));
        T2 = BMC{mod(j-1,nsamps) + 1, jj};
        Td = abs(T1 - T2);
        D(i,j) = sum(Td,'all');
    end
end
figure; imagesc(D)
