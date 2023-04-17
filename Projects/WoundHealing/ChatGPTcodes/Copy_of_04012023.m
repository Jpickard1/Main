% This code was initially written by ChatGPT and then modified by Joshua
%
% PROMPT: Please write MATLAB code to do cell segmentation based on
%         thresholds
%
% Auth: Joshua Pickard and ChatGPT
%       jpic@umich.edu
% Date: April 2, 2023

% v2 = VideoWriter('myFile','Archival');

figure("Position",[10, 10, 810, 810]);
F = cell(1);
for fi=1:length(E)
    f = read(v,fi);
    img = f(100:size(f,1)-100,50:size(f,1)-50,:);
    imagesc(img);
    hold on;
    [y, x] = find(E{fi});
    s = scatter(x, y, '.');
    s.SizeData = 100;
    hold off;
    % pause(0.05)
    
    F{fi} = getframe(gcf);
    % writeVideo(v2, gcf);
end

writerObj = VideoWriter('myVideo.mp4');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F{i};    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

