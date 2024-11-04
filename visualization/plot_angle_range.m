thresh_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/threshold_ranges.txt';
edgels31_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_31.txt';
edgels32_file = '/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_32.txt';

hypothesis_frames = [42, 47];  

fid = fopen(thresh_file, 'r');
data = textscan(fid, 'frame #: %d\nThresh_ore21_1: %f, Thresh_ore21_2: %f\n\n');
fclose(fid);
frame_ids = data{1};
thresh_ore21_1_vals = data{2};
thresh_ore21_2_vals = data{3};

edgels_31_data = read_edgel_file(edgels31_file);
edgels_32_data = read_edgel_file(edgels32_file);

for i = 1:length(frame_ids)
    frame_id = frame_ids(i);

    if ismember(frame_id, hypothesis_frames)
        continue; 
    end

    img_path = sprintf('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/train_img/%d_colors.png', frame_id);
    img = imread(img_path);

    edgels_31 = edgels_31_data{i};
    edgels_32 = edgels_32_data{i};

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


function frames_data = read_edgel_file(filename)
    fid = fopen(filename, 'r');
    frames_data = {};  % Initialize cell array to hold each frame's data
    current_frame = []; % Temporary storage for the current frame data

    while ~feof(fid)
        line = fgetl(fid);
        
        if isempty(line)
            if ~isempty(current_frame)
                frames_data{end+1} = current_frame; % Save current frame data
                current_frame = [];  % Reset for the next frame
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
