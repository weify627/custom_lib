function vr=rotate_2d(v,theta)
% 3d rotation
%random center point, random rotation angle
v = v';
    
center = repmat([0;0;0], 1, size(v,2));

% define a theta degree counter-clockwise rotation matrix
% theta = 2*pi*rand(1);       % pi/3 radians = 60 degrees
R = [cos(theta) -sin(theta) 0; ...
 sin(theta)  cos(theta) 0; ...
         0           0  1];
% pickxy=randperm(2,1);
% s_xy={'xrotate','yrotate'};
% R = makehgtform('xrotate',theta);%s_xy{pickxy},theta);
% R = R(1:3,1:3);

% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);
z_rotated = vo(3,:);
x_rotated=x_rotated';
y_rotated=y_rotated';
z_rotated=z_rotated';
vr=[x_rotated,y_rotated,z_rotated];
end