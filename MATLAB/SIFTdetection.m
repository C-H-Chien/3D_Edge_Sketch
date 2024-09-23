% run('E:/ENGN1610/Lab7/vlfeat-0.9.21/toolbox/vl_setup');
load('MATData/imageArray_10.mat')
% set the threshold
peak_thresh  = 1;
edge_thresh  = 15;

idx = 1;
siftfeature = [];

for (idx = 1 : 50)
    img_current        = imageArray{idx};
    img_current_single = single(rgb2gray(img_current));
    
    % get the descriptors and features
    [featurexy, des] = vl_sift(img_current_single, 'PeakThresh', peak_thresh,'edgethresh',edge_thresh);
    featurexy = featurexy';

    siftfeature{idx} = featurexy;

%     % match the features using matchsift
%     figure(1);
%     imshow(img_current)
%     hold on;
%     plot(featurexy(:,1), featurexy(:,2), 'mx', 'MarkerSize',5, 'LineWidth', 2);
end