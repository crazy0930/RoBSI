function [selectedCorners] = ANMS(image,Points)
%% detect and display FAST features
corners = detectFASTFeatures(image, 'MinQuality',0.0001, 'MinContrast', 0.01);
%% apply SSC algorithm
numRetPoints = Points+100; % number of points you need   750
tolerance = 0.1; % tolerance of the number of return points   0.1
rows = size(image,1);
cols = size(image,2);
[~, ind] = sort(corners.Metric,'descend');
cornersSorted = corners(ind);
% get homegenously distributed keypoints
selectedIdx = SSC(double(cornersSorted.Location), numRetPoints, tolerance, cols, rows);
if size(selectedIdx,2) > Points
    selectedIdx=selectedIdx(:,1:Points);
    selectedCorners = cornersSorted((selectedIdx+1)'); % +1 since matlab is one-indexed
else
    selectedCorners = cornersSorted((selectedIdx+1)'); % +1 since matlab is one-indexed
end
end
%% visualizing results
% figure('Name', 'Selected Keypoints');
% imshow(I); hold on;
% plot(selectedCorners);
