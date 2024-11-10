close all;
clear;

thresh_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/threshold_ranges.txt';
edgels31_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_31.txt';
edgels32_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_32.txt';

hypothesis_frames = [42, 47];  

fid = fopen(thresh_file, 'r');
data = textscan(fid, 'frame #: %d\nThresh_ore31_1: %f, Thresh_ore31_2: %f\nThresh_ore32_1: %f, Thresh_ore32_2: %f\nEpipole Center Hypo1: (%f, %f)\nEpipole Center Hypo2: (%f, %f)\n\n');
fclose(fid);

frame_ids = data{1};
thresh_ore31_1_vals = data{2};
thresh_ore31_2_vals = data{3};
thresh_ore32_1_vals = data{4};
thresh_ore32_2_vals = data{5};
epipole_center_hypo1_x = data{6};
epipole_center_hypo1_y = data{7};
epipole_center_hypo2_x = data{8};
epipole_center_hypo2_y = data{9};

edgels_31_data = read_edgel_file(edgels31_file);
edgels_32_data = read_edgel_file(edgels32_file);

num_frames = min(length(frame_ids), length(edgels_31_data)); 

for i = 1:num_frames
    frame_id = frame_ids(i);

    if ismember(frame_id, hypothesis_frames)
        continue; 
    end

    % Load image
    img_path = sprintf('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/train_img/%d_colors.png', frame_id);
    img = imread(img_path);
    [img_height, img_width, ~] = size(img);

    % Get edgels for current frame
    edgels_31 = edgels_31_data{i};
    edgels_32 = edgels_32_data{i};

    % Plot image and edgels
    figure;
    imshow(img);
    hold on;
    title(sprintf('Threshold Ranges and Edgels for Frame %d', frame_id));

    % Plot edgels for Hypothesis 1 in red
    if ~isempty(edgels_31)
        plot(edgels_31(:, 1), edgels_31(:, 2), 'ro', 'MarkerSize', 4, 'DisplayName', sprintf('Hypothesis 1 Edges on Frame %d', frame_id));
    end
    
    % Plot edgels for Hypothesis 2 in blue
    if ~isempty(edgels_32)
        plot(edgels_32(:, 1), edgels_32(:, 2), 'bo', 'MarkerSize', 4, 'DisplayName', sprintf('Hypothesis 2 Edges on Frame %d', frame_id));
    end

    % Use frame-specific epipole centers
    center1 = [epipole_center_hypo1_x(i), epipole_center_hypo1_y(i)];
    center2 = [epipole_center_hypo2_x(i), epipole_center_hypo2_y(i)];
    
    % Plot wedge for thresh_ore31 using center1
    color1 = 'r';  % Color for the first wedge
    angle1_31 = deg2rad(thresh_ore31_1_vals(i));
    angle2_31 = deg2rad(thresh_ore31_2_vals(i));
    if isfinite(center1(1)) && isfinite(center1(2))
        drawWedgeAcrossBoundaries(center1, angle1_31, angle2_31, img_width, img_height, color1, 'Hypothesis 1');
    end

    % Plot wedge for thresh_ore32 using center2
    color2 = 'b';  % Color for the second wedge
    angle1_32 = deg2rad(thresh_ore32_1_vals(i));
    angle2_32 = deg2rad(thresh_ore32_2_vals(i));
    if isfinite(center2(1)) && isfinite(center2(2))
        drawWedgeAcrossBoundaries(center2, angle1_32, angle2_32, img_width, img_height, color2, 'Hypothesis 2');
    end

    % Set up legend without unwanted entries
    legend();
    hold off;
    pause(1);
end



function frames_data = read_edgel_file(filename)
    fid = fopen(filename, 'r');
    frames_data = {};  
    current_frame = []; 

    while ~feof(fid)
        line = fgetl(fid);
        
        if isempty(line)
            if ~isempty(current_frame)
                frames_data{end+1} = current_frame; 
                current_frame = [];  
            end
        else
            values = sscanf(line, '%f %f');
            current_frame = [current_frame; values'];
        end
    end
    
    if ~isempty(current_frame)
        frames_data{end+1} = current_frame;
    end

    fclose(fid);
end


function drawWedgeAcrossBoundaries(center, angle1, angle2, img_width, img_height, color, legend_text)
    [x1, y1] = calculateBoundaryPoint(center, angle1, img_width, img_height);
    [x2, y2] = calculateBoundaryPoint(center, angle2, img_width, img_height);
    plot([center(1), x1], [center(2), y1], 'Color', color, 'LineWidth', 1.5, 'DisplayName', legend_text);
    plot([center(1), x2], [center(2), y2], 'Color', color, 'LineWidth', 1.5, 'DisplayName', legend_text);
end


function [x, y] = calculateBoundaryPoint(center, angle, img_width, img_height)
    x_limits = [0, img_width];
    y_limits = [0, img_height];
    
    y_intercept1 = center(2) + (x_limits(1) - center(1)) * tan(angle);
    y_intercept2 = center(2) + (x_limits(2) - center(1)) * tan(angle);
    x_intercept1 = center(1) + (y_limits(1) - center(2)) / tan(angle);
    x_intercept2 = center(1) + (y_limits(2) - center(2)) / tan(angle);

    boundary_points = [x_limits(1), y_intercept1; x_limits(2), y_intercept2; x_intercept1, y_limits(1); x_intercept2, y_limits(2)];
    valid_points = boundary_points(boundary_points(:, 1) >= 0 & boundary_points(:, 1) <= img_width & boundary_points(:, 2) >= 0 & boundary_points(:, 2) <= img_height, :);
    
    distances = pdist2(valid_points, valid_points);
    [~, idx] = max(distances(:));
    [row, col] = ind2sub(size(distances), idx);
    x = [valid_points(row, 1), valid_points(col, 1)];
    y = [valid_points(row, 2), valid_points(col, 2)];

end
