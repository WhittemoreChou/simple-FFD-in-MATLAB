clear;clc;close all;
%% 1. Data Preparation
% Generate evenly spaced points on a 3D plane
x = 200:20:380;
y = 200:20:380;
[X,Y] = meshgrid(x,y);
Z = ones(size(X)); % Flat plane at Z=1;

% Convert to column vectors for FFD processing
sourcePoints = [X(:), Y(:), Z(:)]';

figure
% Plot the source model
scatter3(sourcePoints(1,:), sourcePoints(2,:), sourcePoints(3,:), 'filled');
xlabel('X'), ylabel('Y'), zlabel('Z');
title('Original Points'); axis equal;

%% 2. FFD Model Setup
% Define transformation parameters
translation = [190, 190, -2];  % Translation vector
rotation = [0, 0, 0];          % Rotation angles (degrees)
scale = [200, 200, 5];         % Scaling factors
tform = rigidtform3d(rotation,translation).A*diag([scale 1]);

% Initialize FFD model
controlPointNumbers = [3,3,1];
ffd = FFD(tform,controlPointNumbers);
hold on
% Visualize initial FFD lattice
ffd.draw(false);
axis equal

%% 3. Apply Deformation
% Select specific control point to displace
targetIndex = ffd.GridToStackIndices(2, 1, 1); % Indexes on three axes [i,j,k]
displacements = zeros(size(ffd.OriginalControlPoints));
displacements(3,targetIndex) = 50; % Pull up one point

% Perform deformation
deformedPoints = ffd.deform(sourcePoints, displacements);

figure(2)
% Plot the result model
scatter3(deformedPoints(1,:),deformedPoints(2,:),deformedPoints(3,:),'filled')
hold on
% Plot the FFD model
ffd.draw(true);
axis equal