addpath(genpath('/home/spmoran/indikar/temp-sean/CellTracking_Current'));%path to 3rd party packages
% track
VideoSize = [1000 1000];

densityMatrixRows = 50;
densityMatrixCols = 50;

mainDir ='/home/spmoran/indikar/temp-sean/CellTracking_Current/hold_czi/';%path to CZI

% get list of subdirectories
subDirs = dir([mainDir '/*.czi']);
%subDirs(1:2) = []; % remove '.' and '..'
subDirs = {subDirs.name};

scene = 1;
numOfChannels = 5; %phase, green, red, blue

hblob = vision.BlobAnalysis( ...
    'AreaOutputPort', false, ...
    'BoundingBoxOutputPort', false, ...
    'OutputDataType', 'single', ...
    'MinimumBlobArea', 20, ...
    'MaximumBlobArea', 1000, ...
    'MaximumCount', 10000);

% Acknowledgement
ackText = ['University of Michigan'];

%hVideo = vision.VideoPlayer;
%hVideo.Name  = 'Results';
%hVideo.Position(1) = round(hVideo.Position(1));
%hVideo.Position(2) = round(hVideo.Position(2));
%hVideo.Position([4 3]) = 30+VideoSize;

%no display compatible video mode


frameCount = int16(1);

% store as cell array (for each file)
% row, col, frame
DensityMatrices = cell(1,length(subDirs));
InfoMatrices = cell(1,length(subDirs));

%for fileNum=1:length(subDirs)
fileNum = 1;
% videowriter for video of marked cells
outputVideo= VideoWriter(['markedCells',subDirs{fileNum}]);
outputVideo.FrameRate = 5;
open(outputVideo);
filename = [mainDir,'/',subDirs{fileNum}];
data = bfopen(filename);

[numRows, numCols] = size(data{1,1}{1,1});
numPlanes = length(data{1,1}(:,1));


disp(numRows);
disp(numCols);

numFrames = round(numPlanes/numOfChannels);
DensityMatrices{fileNum} = zeros(densityMatrixRows, densityMatrixCols, numFrames);
InfoMatrices{fileNum} = cell(1, numFrames);
named_string = 0;
frame_tracking = 1;
frame_step = 4;
recall_table = double.empty(0,6);
memory_table = double.empty(0,6);
time_series_dictionary = containers.Map;
MatchingMatrices{fileNum} = cell(1, round(numFrames/frame_step));


for frame=4:frame_step:numFrames
    % Read input video frame
    %image = step(hvfr);
    %figure, imshow(image);
    A = getFrame(frame,data,numOfChannels);
    % channels
    blueChan = mat2gray(A(:,:,4), [0 2^16-1]); % blue channel
    greenChan = mat2gray(A(:,:,2), [0 2^16-1]); % green channel
    redChan = mat2gray(A(:,:,3), [0 2^16-1]); % red channel
    
    % remove background from channels
    %greenChan = greenChan - imopen(greenChan, strel('disk',20));
    %redChan = redChan - imopen(greenChan, strel('disk',20));
       
    %figure, imshow(greenChan, []);
    %figure, imshow(redChan, []);
    %greenChan = imadjust(greenChan);
    %redChan = imadjust(redChan);
    
    %figure, imshow(image);
    
    bgBlue = imopen(blueChan, strel('disk',15));
    
    %figure, imshow(bg);
    
    I2 = blueChan - bgBlue;
    %figure, imshow(I2);
    
    I3 = imadjust(I2);
    %figure, imshow(I3);
    
    % Determine threshold using Otsu's method
    BW = imbinarize(I3);
    %figure, imshow(BW);
    
    % get info on average intensities
    cc = bwconncomp(BW);
    
    % remove area <50 and area >350
    %for i=1:cc.NumObjects
    %    if (length(cc.PixelIdxList{i}) < 30 || length(cc.PixelIdxList{i}) >250)
    %        BW(cc.PixelIdxList{i}) = 0;
    %    end
    %end
        
    stats = regionprops(cc); % what do these outputs looklike
    
    % store information in an array of
    % x center, y center, area, avg green, avg red, avg blue
    info = zeros(cc.NumObjects,6);
    info(:,1:2) = cat(1,stats.Centroid);
    info(:,3) = cat(1, stats.Area);
    for i=1:cc.NumObjects
        info(i,4) = median(greenChan(cc.PixelIdxList{i}));
        info(i,5) = median(redChan(cc.PixelIdxList{i}));
        info(i,6) = median(blueChan(cc.PixelIdxList{i}));
    end
    
    % compute Delaunay trianglation of centroids
    DT = delaunay(info(:,1), info(:,2));
    TR = triangulation(DT,info(:,1:2));
    %f = figure('visible', 'off');
    f = triplot(TR);
    figurename = sprintf('Figure%d',named_string);
    saveas(f,figurename);%print -djpeg  
    %close(f)   
    %figure, triplot(TR);
    
    [imrows, imcols] = size(BW);
    DensityMatrices{fileNum}(:,:,frame) = getDensityMatrix(info(:,1:2),imrows,imcols,densityMatrixRows,densityMatrixCols);
    InfoMatrices{fileNum}{frame} = info;
   
    scatter_figure_name = sprintf('RedGreen_Figure%d',named_string); 
    g = scatter(info(:,4)/max(info(:,4)), info(:,5)/max(info(:,5)));
    saveas(g,scatter_figure_name);
    %xlabel('green'), ylabel('red');

    Centroid = step(hblob, BW);   % Calculate the centroid
    numBlobs = size(Centroid,1);  % and number of cells.

    cell_matches = double.empty(0,6);
    for i=1:size(info,1)
            if sum(ismember(round(info(i,1:2)),round(Centroid))) == 2 && info(i,3) > 135%filter out small debris
                cell_matches = cat(1,cell_matches, info(i,:));
            else
                continue
            end
    end
%     dlmwrite(sprintf("cell_matches%d.txt",named_string),cell_matches, 'delimiter', '\t');
    
    MatchingMatrices{fileNum}{named_string+1}=cell_matches;
    
    % Display the number of frames and cells.
    frameBlobTxt = sprintf('Frame %d, Count %d', frameCount, numBlobs);
    %blueChan = insertText(blueChan, [1 1], frameBlobTxt, ...
    %    'FontSize', 16, 'BoxOpacity', 0, 'TextColor', 'white');
    %blueChan = insertText(blueChan, [1 size(blueChan,1)], ackText, ...
    %    'FontSize', 10, 'AnchorPoint', 'LeftBottom', ...
    %    'BoxOpacity', 0, 'TextColor', 'white');
    
    % Display video (replaced image with y3)
    image_out = insertMarker(mat2gray(blueChan), Centroid, '*', 'Color', 'red');

    writeVideo(outputVideo, image_out);
    %step(hVideo, image_out);
    
    named_string = named_string +1;
    frameCount = frameCount + 1;
       
end

%release(hVideo);
%end

%assignment and tracking loop
for f=2:frameCount
    trackingMatrix{fileNum}{named_string+1}=zeros(1,2); %construct in main loop, assign after
    D = computeHungarianMatrix(MatchingMatrices{1,1}{1,f-1},MatchingMatrices{1,1}{1,f});%assignment
    Z = matchpairs(D, 1000);
    Z = sortrows(Z,1);
    trackingMatrix{1,1}{1,f-1}=Z;
    disp(size(Z,1))
end


celldisp(MatchingMatrices)
celldisp(trackingMatrix)
%%%Write Cell 13 Tracking structure, just change seekcell to change the
%%%indexed cell to follow
%seekcell = 13;
track13 = double.empty(1,0);
for i=1:frameCount-2
    disp(i)    
    disp(seekcell)
    rowval = find(trackingMatrix{1,1}{1,i}(:,1) == seekcell);%each of these compares two frames
    %so it will always have n-1 frames
    disp("rowval")
    disp(rowval)
    disp("a1:a2")
    disp(trackingMatrix{1,1}{1,i}(rowval,:))
    a(1:2) = trackingMatrix{1,1}{1,i}(rowval,:);
    seekcell = a(2);
    b(1:2) = MatchingMatrices{1,1}{1,i}(a(1),1:2);
    track13 = cat(2,track13,b);
    if i==frameCount-2 %collect last frame
        b(1:2) = MatchingMatrices{1,1}{1,i+1}(a(2),1:2);
        track13 = cat(2,track13,b);
    end
end


%%%Extract Cell 13 tracked images
pseudoFrame = int16(1);
alpha = 20;
blueContainerMtx = zeros(2*alpha+1,10);
greenContainerMtx = zeros(2*alpha+1,10);
redContainerMtx = zeros(2*alpha+1,10);

for frame=4:frame_step:numFrames
	A = getFrame(frame,data,numOfChannels);
    blueChan = mat2gray(A(:,:,4), [0 2^16-1]); % blue channel
    greenChan = mat2gray(A(:,:,2), [0 2^16-1]); % green channel
    redChan = mat2gray(A(:,:,3), [0 2^16-1]); % red channel
	bgBlue = imopen(blueChan, strel('disk',15));
    bgGreen = imopen(greenChan, strel('disk',20));
    bgRed = imopen(redChan, strel('disk',20));
    
    I2b = blueChan - bgBlue;
    I2g = greenChan - bgGreen;
    I2r = redChan - bgRed;
    %figure, imshow(I2);
    
    I3b = imadjust(I2b);     
    I3g = imadjust(I2g);
    I3r = imadjust(I2r);

    centroidXY = round(track13(pseudoFrame*2-1:pseudoFrame*2));%form rowvector and slice by pairs
    
    blueContainerMtx = cat(2,blueContainerMtx,I3b(centroidXY(2)-alpha:centroidXY(2)+alpha,centroidXY(1)-alpha:centroidXY(1)+alpha));
    blueContainerMtx = cat(2,blueContainerMtx, zeros(size(blueContainerMtx,1),7));
    
    greenContainerMtx = cat(2,greenContainerMtx,I3g(centroidXY(2)-alpha:centroidXY(2)+alpha,centroidXY(1)-alpha:centroidXY(1)+alpha));
    greenContainerMtx = cat(2,greenContainerMtx, zeros(size(greenContainerMtx,1),7));
    
    redContainerMtx = cat(2,redContainerMtx,I3r(centroidXY(2)-alpha:centroidXY(2)+alpha,centroidXY(1)-alpha:centroidXY(1)+alpha));
    redContainerMtx = cat(2,redContainerMtx, zeros(size(redContainerMtx,1),7));

%try, make 3 matrices starting at zeros(height, 20) to save red, blue and green. 
%For every iteration, cat(2, starterMatrix, Z), where
%Z = I3b(centroidXY(2)-a:centroidXY(2)+a,centroidXY(1)-a:centroidXY(1)+a)
%and again cat(2, starterMatrix, zeros(height, 20))
%
%
%To get color, then do an elementwise multiplication against a mask of 1's
%blue, 0's red, and 0's green, to get a blue only matrix.
% 
%     f = imshow(I3b(centroidXY(2)-a:centroidXY(2)+a,centroidXY(1)-a:centroidXY(1)+a));
%     figurename = sprintf('BlueWindows_%d',pseudoFrame);
%     saveas(f,figurename);%print -djpeg 
% 
%     f = imshow(I3r(centroidXY(2)-a:centroidXY(2)+a,centroidXY(1)-a:centroidXY(1)+a));
%     figurename = sprintf('RedWindows_%d',pseudoFrame);
%     saveas(f,figurename);%print -djpeg  
% 
%     f = imshow(I3g(centroidXY(2)-a:centroidXY(2)+a,centroidXY(1)-a:centroidXY(1)+a));
%     figurename = sprintf('GreenWindows_%d',pseudoFrame);
%     saveas(f,figurename);%print -djpeg  
    pseudoFrame = pseudoFrame + 1;
end

%just vertical cat(1, ..., ...) instead?

b_f = imshow(blueContainerMtx);
figurename = 'BlueWindows';
saveas(b_f,figurename);%print -djpeg

r_f = imshow(greenContainerMtx);
figurename = 'GreenWindows';
saveas(r_f,figurename);%print -djpeg

g_f = imshow(redContainerMtx);
figurename = 'RedWindows';
saveas(g_f,figurename);%print -djpeg

full_windows = cat(1,blueContainerMtx, zeros(7,size(blueContainerMtx,2)));
full_windows = cat(1,full_windows, greenContainerMtx);
full_windows = cat(1,full_windows, zeros(7,size(blueContainerMtx,2)));
full_windows = cat(1,full_windows, redContainerMtx);
full_name = "full_windows";
full_image = imshow(full_windows);
saveas(full_image, full_name);%print -djpeg


%%%% functions %%%%
function C = computeHungarianMatrix(frameT,frameTp1)
    C = zeros(size(frameT,1),size(frameTp1,1));
    for i=1:size(frameT,1)
        for j=1:size(frameTp1,1)
            C(i,j)=solve_penalty(frameT(i,:),frameTp1(j,:));
        end
    end
end

% return cell density matrix of specified size from array of centroids
function A = getDensityMatrix(centroids, imrows, imcols, rows, cols)
A = zeros(rows, cols);

[cells, ~] = size(centroids);
rowSize = imrows/rows; colSize = imcols/cols;

for i=1:cells
    x = ceil(centroids(i,1)/colSize);
    y = ceil(centroids(i,2)/rowSize);
    A(y,x) = A(y,x)+1;
end
end

% return scene as array of (x,y,t,channel)
function A = getFrame(frame,data,numOfChannels)

[numRows, numCols] = size(data{1,1}{1,1});
%numPlanes = length(data{1,1}(:,1));

A = zeros(numRows, numCols, numOfChannels);
for j=1:numOfChannels
    plane = (frame-1)*numOfChannels + j;
    A(:,:,j) = data{1,1}{plane,1};
end
end

function [combined]=solve_penalty(Frame1, Frame2)
    if Frame1(3) > 600
        a = 0.5;
        b = 2.0;
    elseif Frame1(3) > 800
        a = 0.5;
        b = 10;
    elseif 400<Frame1(3)<500
        a = 2.0;
        b = 0.33;
    else
        a = 1;
        b = 1;
    end
    temp_area = abs(Frame2(3)-Frame1(3));
    temp_distance = sqrt((Frame2(1)-Frame1(1))^2+(Frame2(2)-Frame1(2))^2);
    combined = a*temp_area + b*temp_distance;
end
