function [min_dist, min_dist_ind]=minimum_distance(x,y)
% x: Nx3 matrix (N points, x-y-z)
% y: Mx3 matrix (M points, x-y-z)

% min_dist: Nx1 vector (minimum distance of each element of x from y)
% min_dist_ind: Nx1 vector (index of the y-point that is the closest to x)
min_dist=nan(size(x,1),1);
min_dist_ind=nan(size(x,1),1);

for i=1:size(x,1)
    i_dist=[];
    for j=1:size(y,1)
        i_dist(j,1)=pdist([x(i,:); y(j,:)]);
    end
    [min_d, min_ind]=min(i_dist);
    min_dist(i,1)=min_d;
    min_dist_ind(i,1)=min_ind;
end
