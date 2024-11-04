thresh_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/threshold_ranges.txt';

fid = fopen(thresh_file, 'r');
data = textscan(fid, 'Image ID: %d\nThresh_ore21_1: %f, Thresh_ore21_2: %f\n\n');
fclose(fid);

img_ids = data{1};
thresh_ore21_1_vals = data{2};
thresh_ore21_2_vals = data{3};

for i = 1:length(img_ids)
    img_id = img_ids(i);
    img_path = sprintf('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/train_img/%d_colors.png', img_id);
    img = imread(img_path);

    figure;
    imshow(img);
    hold on;
    title(sprintf('Threshold Range for Image %d', img_id));

    text(10, 20, sprintf('Thresh Range: [%.2f, %.2f]', thresh_ore21_1_vals(i), thresh_ore21_2_vals(i)), 'Color', 'red', 'FontSize', 12);
    hold off;
end
