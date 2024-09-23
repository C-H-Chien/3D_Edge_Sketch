function [imgData, sceneData, Params] = readCurveSyntheticDataset(view, dataPath, selectCurveNum, selectCurveId)
    
    %> import 2D data
    points2D = load([dataPath 'frame_' num2str(view,'%04.f') '-pts-2D.txt']);
    tangent = load([dataPath 'frame_' num2str(view,'%04.f') '-tgts-2D.txt']);

    %> import camera parameters
    pose = load([dataPath 'frame_' num2str(view,'%04.f') '.extrinsic']);
    K = load([dataPath '\calib.intrinsic']);
    
    %> import 3D data
    points3D = load([dataPath '\crv-3D-pts.txt']);
    tangent3D = load([dataPath '\crv-3D-tgts.txt']);
    curveId = load([dataPath '\crv-ids.txt']);
    
    R = pose(1:3,1:3);
    C = transpose(pose(4,:));
    
    %temp = randperm(size(points2D,1));
    
    %> picked 2D points indices
    %pickedPointsIndx = temp(1:numPoints);
    %points2DPick = points2D(pickedPointsIndx,:);
    %tangents2DPick = tangent(pickedPointsIndx,:);
    points2DPick = points2D;
    tangents2DPick = tangent;

    %> 3D points obtained from 2D picked points indices
    %points3DPick = points3D(pickedPointsIndx,:);
    %tangents3DPick = tangent3D(pickedPointsIndx,:);
    points3DPick = points3D;
    tangents3DPick = tangent3D;
    
    %> return to the outputs
    Params.R = R;
    Params.T = -R*C;
    Params.K = K;
    Params.C = C;
    
    if selectCurveNum == 0
        imgData.points = points2DPick;
        imgData.tangents = tangents2DPick;
        sceneData.points = points3DPick;
        sceneData.tangents = tangents3DPick;
        sceneData.curveId = curveId;
    else
        %> choose 3D curve(s) randomly from all 39 curves and import the image and scene data
        %rndSelectCrvId = round(unifrnd(0, 38, [selectCurveNum,1]));
        rndSelectCrvId = selectCurveId;
        
        selectCurvePointsIndx = find(curveId == rndSelectCrvId);
        selectCurvePoints2D = points2DPick(selectCurvePointsIndx, :);
        selectCurveTangents2D = tangents2DPick(selectCurvePointsIndx, :);
        selectCurvePoints3D = points3DPick(selectCurvePointsIndx, :);
        selectCurveTangents3D = tangents3DPick(selectCurvePointsIndx, :);
        selecCurveId = rndSelectCrvId;
        
        imgData.points     = selectCurvePoints2D;
        imgData.tangents   = selectCurveTangents2D;
        sceneData.points   = selectCurvePoints3D;
        sceneData.tangents = selectCurveTangents3D;
        sceneData.curveId  = selecCurveId;
    end
    
%     R2 = pose2(1:3,1:3);
%     C2 = transpose(pose2(4,:));
%     R3 = pose3(1:3,1:3);
%     C3 = transpose(pose3(4,:));
%     R12 = R2 * inv(R1);
%     T12 = R2 * (C1 - C2);
%     R23 = R3 * inv(R2);
%     T23 = R3 * (C2 - C3);

%     P1 = [eye(3,3) [0;0;0]];
%     P2 = [R12 T12];
%     P3 = [R23 * R12 R23 * T12 + T23];
    
%     tangents.tangents2 = tangent2Pick;
%     tangents.tangents3 = tangent3Pick;

%     points2Pick = points2(pickedPointsIndx,:);
%     points3Pick = points3(pickedPointsIndx,:);

%     tangent2Pick = tangent2(pickedPointsIndx,:);
%     tangent3Pick = tangent3(pickedPointsIndx,:);

%     parameters.R12 = R12;
%     parameters.R23 = R23;
%     parameters.R13 = R23 * R12;
%     parameters.T12 = T12;
%     parameters.T23 = T23;
%     parameters.T13 = R23 * T12 + T23;
%     parameters.K = K;
%     parameters.quatR12 = rotm2quat(R12);
%     parameters.quatR23 = rotm2quat(R23);
%     parameters.quatR13 = rotm2quat(R23 * R12);
%     parameters.P1 = P1;
%     parameters.P2 = P2;
%     parameters.P3 = P3;

%     points.points2 = points2Pick;
%     points.points3 = points3Pick;
end