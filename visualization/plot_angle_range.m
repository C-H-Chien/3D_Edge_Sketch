
thresh_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/threshold_ranges.txt';
edgels31_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_31.txt';
edgels32_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_32.txt';

fid = fopen(thresh_file, 'r');
data = textscan(fid, 'frame #: %d\nThresh_ore21_1: %f, Thresh_ore21_2: %f\n\n');
fclose(fid);
frame_ids = data{1};
thresh_ore21_1_vals = data{2};
thresh_ore21_2_vals = data{3};

edgels_31_data = dlmread(edgels31_file);
edgels_32_data = dlmread(edgels32_file);

num_frames = length(frame_ids);

rows_per_frame_31 = floor(size(edgels_31_data, 1) / num_frames);
rows_per_frame_32 = floor(size(edgels_32_data, 1) / num_frames);

for i = 1:num_frames
    frame_id = frame_ids(i);
    
    img_path = sprintf('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/train_img/%d_colors.png', frame_id);
    img = imread(img_path);
    
    start_idx_31 = (i-1) * rows_per_frame_31 + 1;
    end_idx_31 = min(i * rows_per_frame_31, size(edgels_31_data, 1));
    
    start_idx_32 = (i-1) * rows_per_frame_32 + 1;
    end_idx_32 = min(i * rows_per_frame_32, size(edgels_32_data, 1));
    
    if start_idx_31 <= end_idx_31
        edgels_31 = edgels_31_data(start_idx_31:end_idx_31, :);
    else
        edgels_31 = [];
    end
    
    if start_idx_32 <= end_idx_32
        edgels_32 = edgels_32_data(start_idx_32:end_idx_32, :);
    else
        edgels_32 = [];
    end
    
    figure;
    imshow(img);
    hold on;
    title(sprintf('Threshold Range and Edgels for Frame %d', frame_id));

    text(10, 20, sprintf('Thresh Range: [%.2f, %.2f]', thresh_ore21_1_vals(i), thresh_ore21_2_vals(i)), 'Color', 'red', 'FontSize', 12);
    
    if ~isempty(edgels_31)
        plot(edgels_31(:, 1), edgels_31(:, 2), 'go', 'MarkerSize', 4, 'DisplayName', 'Edgels 31');
    end
    
    if ~isempty(edgels_32)
        plot(edgels_32(:, 1), edgels_32(:, 2), 'bo', 'MarkerSize', 4, 'DisplayName', 'Edgels 32');
    end
    
    legend;
    hold off;

    pause(0.5);
end
