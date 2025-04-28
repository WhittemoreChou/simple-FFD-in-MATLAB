classdef FFD < handle
    % FFD Class implementing Free-Form Deformation (FFD) for 3D models.
    
    properties (SetAccess = private)
        % Core transformation parameters
        AffineTransform   % 4x4 affine transformation matrix (global to local frame)
        ControlPointNumbers     % 1x3 vector [nx, ny, nz] specifying control points per axis
        LocalFrameOrigin        % 1x3 vector: origin of local frame
        LocalFrameAxes          % 3x3 matrix: local coordinate axes

        % Control point data
        GridToStackIndices  % 3xN matrix: grid indices to stack index
        StackToGridIndices  % Stack index to grid index
        OriginalControlPoints  % 3xN matrix: control points in global frame
        DeformedControlPoints  % 3xN matrix: deformed control points in global frame
    end
    
    methods
        %% Constructor
        function obj = FFD(affineTransform, controlPointNumbers)
            % FFD constructor
            %
            % Inputs:
            %   affineTransform: 4x4 affine transformation matrix (global to local)
            %   controlNumbers: 1x3 vector specifying control points per axis

            % Validate input parameters
            validateattributes(affineTransform, {'double'},{'size',[4 4]});
            validateattributes(controlPointNumbers, {'double'},{'size',[1 3]});

            % Assign properties
            obj.AffineTransform = affineTransform;
            obj.ControlPointNumbers = controlPointNumbers;
            obj.LocalFrameOrigin = affineTransform(1:3,4);
            obj.LocalFrameAxes = affineTransform(1:3,1:3)';

            % Initialize control point grid
            [I, J, K] = ndgrid(1:controlPointNumbers(1), ...
                               1:controlPointNumbers(2), ...
                               1:controlPointNumbers(3));
            obj.StackToGridIndices = [I(:), J(:), K(:)]';  % 3×n矩阵
            obj.GridToStackIndices = reshape(1:prod(controlPointNumbers), controlPointNumbers);

            % Generate original control points
            grid = [(I(:)-1)/(obj.positiveRegulation(controlPointNumbers(1)-1)), ...
                    (J(:)-1)/(obj.positiveRegulation(controlPointNumbers(2)-1)), ...
                    (K(:)-1)/(obj.positiveRegulation(controlPointNumbers(3)-1))];
            obj.OriginalControlPoints = repmat(obj.LocalFrameOrigin, 1, prod(obj.ControlPointNumbers)) ...
                                        + (grid * obj.LocalFrameAxes)';
            obj.DeformedControlPoints = obj.OriginalControlPoints; % Initialize to original
        end
        %% Deformation Method
        function [deformedPoints] = deform(obj, sourcePoints, displacementVectors)
            % Apply FFD deformation to source points
            %
            % Inputs:
            %   sourcePoints: 3xM matrix of points to deform
            %   displacementVectors: 3xN matrix of control point displacements
            %
            % Output:
            %   deformedPoints: 3xM matrix of deformed points
            
            % Validate inputs
            obj.DeformedControlPoints = obj.OriginalControlPoints + displacementVectors;
            % Transform source points to local frame
            localPoints = obj.transformToLocal(sourcePoints);

            % Compute Bernstein basis matrices for each axis
            bernMatrices = obj.computeBernsteinMatrices(localPoints);

            % Compute mapping matrix via tensor product
            mappingMatrix = obj.computeMappingMatrix(bernMatrices);

            % Apply deformation
            deformedPoints = obj.DeformedControlPoints * mappingMatrix';
        end

        function draw(obj, isDeformed)
            % Draw control points and links
            %
            % Inputs:
            %   isDeformed: boolean to show deformed/original points
            %   varargin: optional parameters (e.g., lineColor, pointColor)

            lineColor = [0 0 0]; % Black
            pointColor = [0.8,0,0]; % Red
            pointSize = 6; % Points size
            
            if isDeformed && ~isempty(obj.DeformedControlPoints)
                controlPoints = obj.DeformedControlPoints;                
            else
                controlPoints = obj.OriginalControlPoints;
            end
            
            for j = 1:obj.ControlPointNumbers(2)
                for k = 1:obj.ControlPointNumbers(3)
                    X = controlPoints(1,obj.GridToStackIndices(:,j,k));
                    Y = controlPoints(2,obj.GridToStackIndices(:,j,k));
                    Z = controlPoints(3,obj.GridToStackIndices(:,j,k));
                    plot3(X,Y,Z,'color',lineColor);
                    hold on
                end
            end
        
            for i = 1:obj.ControlPointNumbers(1)
                for k = 1:obj.ControlPointNumbers(3)
                    X = controlPoints(1,obj.GridToStackIndices(i,:,k));
                    Y = controlPoints(2,obj.GridToStackIndices(i,:,k));
                    Z = controlPoints(3,obj.GridToStackIndices(i,:,k));
                    plot3(X,Y,Z,'color',lineColor);
                    hold on
                end
            end
        
            for i = 1:obj.ControlPointNumbers(1)
                for j = 1:obj.ControlPointNumbers(2)
                    X = controlPoints(1,obj.GridToStackIndices(i,j,:));
                    Y = controlPoints(2,obj.GridToStackIndices(i,j,:));
                    Z = controlPoints(3,obj.GridToStackIndices(i,j,:));
                    plot3(X,Y,Z,'color',lineColor);
                    hold on
                end
            end         
            scatter3(controlPoints(1,:), controlPoints(2,:), controlPoints(3,:),pointSize,'filled',MarkerEdgeColor=pointColor,MarkerFaceColor=pointColor);
        end
    end
    methods (Access = private)
        %% Private helper functions
        function [correctedNumber] = positiveRegulation(obj,number)
            if number > 0; correctedNumber = number; else correctedNumber = 1; end
        end
        function localPoints = transformToLocal(obj, globalPoints)
            % Transform points from global to local frame
            % Add homogeneous coordinate
            homogPoints = [globalPoints; ones(1, size(globalPoints,2))];
            transformed = obj.AffineTransform \ homogPoints;
            localPoints = transformed(1:3, :);
        end
        function bernMatrices = computeBernsteinMatrices(obj, localPoints)
            % Compute Bernstein basis matrices for each axis
            bernMatrices = cell(1,3);
            for dim = 1:3
                bernMatrices{dim} = bernsteinMatrix(obj.ControlPointNumbers(dim)-1,localPoints(dim,:));
            end
        end
        function mappingMatrix = computeMappingMatrix(obj, bernMatrices)
            % Compute tensor product of Bernstein matrices
            A = bernMatrices{1}; B = bernMatrices{2}; C = bernMatrices{3};
            triBernTensorProd = A .* permute(B, [1 3 2]) .* permute(C, [1 3 4 2]);
            mappingMatrix = reshape(triBernTensorProd, [], prod(obj.ControlPointNumbers));
        end
    end
end
