function [diff22,closest,distance, currentIdx,startIdx,endIdx] = getClosestInPathWindow(trajectory, robotPos, prevIdx, winsize)
    % trajectory: Nx3 matrix of path points
    % robotPos: 1x3 current position
    % prevIdx: The index found in the last iteration
    % winSize: How many points to look ahead/behind (e.g., 20)

    [dim,numpoints] = size(trajectory);
     
    % 1. Define window boundaries (clamped to trajectory limits)
    startIdx = max(1, prevIdx - fix(winsize/2));
    endIdx   = startIdx+winsize;
    if endIdx > numpoints
        endIdx= numpoints;
        startIdx= numpoints-winsize;
    end
        
    %endIdx   = min(numPoints, prevIdx + idivide(winsize, 2);
    
    % 2. Extract the local search area
    searchWindow = trajectory(:, startIdx:endIdx);
     
    % 3. Find closest point in this local slice (Brute Force is fastest here)
    if dim==2
        diff2 = [searchWindow;zeros(1,winsize+1)] - [robotPos;0];
        diff2=diff2.^2;
        for iji=1:dim
            diff2(dim+1,:)=diff2(dim+1,:)+diff2(iji,:);
        end
        [distance, localIdx] = min(diff2(dim+1,:)); 
        diff22= diff2(:,localIdx);
    end
    
    if dim==3
        diff3 = [searchWindow;zeros(1,winsize+1)] - [robotPos;0];
        diff3=diff3.^2;
        for iji=1:dim
            diff3(dim+1,:)=diff3(dim+1,:)+diff3(iji,:);
        end
        [distance, localIdx] = min(diff3(dim+1,:));
        diff22= diff3(:,localIdx);
    end 
    distance = sqrt(distance);    
    
    %[distance, localIdx] = min(distSq);
    
    % 4. Convert local index back to global trajectory index
    currentIdx = startIdx + localIdx - 1;
    closest = trajectory(:, currentIdx);
end