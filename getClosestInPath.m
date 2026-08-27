function [distance] = getClosestInPath(trajectory, robotPos)
    disty= trajectory-robotPos;
    disty= disty.^2;
    [dim,numpoints] = size(trajectory);
    distz=zeros(1,numpoints);
    for i=1:dim
        distz= distz+disty(:,i);
    end
    [distance, localIdx] = min(distz(:));
end
    
    