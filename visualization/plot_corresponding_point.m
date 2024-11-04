% Define the reference point in image 48
reference_point = [520.5543 426.4240];
distances = sqrt((points1(:,1) - reference_point(1)).^2 + (points1(:,2) - reference_point(2)).^2);
distance_threshold = 0.005;
filtered_indices = find(distances <= distance_threshold);

filtered_points1 = points1(filtered_indices, :);
filtered_points2 = points2(filtered_indices, :);

img1 = imread(image_file1);
img2 = imread(image_file2);

figure('Position', [100, 100, 1400, 600]);


subplot(1, 2, 1);
imshow(img1);
hold on;
scatter(filtered_points1(:,1), filtered_points1(:,2), 'r.', 'SizeData', 50);
title('Image 48 - Filtered Points Near Reference');

subplot(1, 2, 2);
imshow(img2);
hold on;
scatter(filtered_points2(:,1), filtered_points2(:,2), 'r.', 'SizeData', 50);
title('Image 43 - Corresponding Filtered Points');

savefig('edge_visualization_filtered_near_point.fig');
print('edge_visualization_filtered_near_point', '-dpng', '-r300');
