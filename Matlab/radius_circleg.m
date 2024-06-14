function [x,y]=radius_circleg(bs,s)
%CIRCLEG        Gives geometry data that describes a unit disk centered at
%               the origin.
%
%  CIRCLEG returns the geometry data for a disk of radius one that is
%  centered at the origin. 
%
% The domain has four boundary segments, one in each of the four quadrants
% of the Cartesian plane. Each segment is parameterized from 0 to 1 in a
% counter-clockwise direction. For example, the first boundary segment
% goes from (x,y)=(1,0) (when the local parameter s=0) to (x,y)=(0,1) (when
% the local parameter s=1). 
% 
%  CIRCLEG obeys the following three conventions:
%
%  NE=CIRCLEG assigns the value 4 to NE, the scalar giving the number of
%  boundary segments. 
%
%  D=CIRCLEG(BS), where BS is a vector of integers between 1 and 4,
%  returns a matrix with one column for each boundary segment in BS. Each
%  column is [0;1;1;0]. This is, 
%         Row 1 contains the start parameter value for each segment.
%         Row 2 contains the end parameter value for each segment.
%         Row 3 contains the number of the left-hand region for each
%         segment. The value 1 corresponds to the interior of the circle.
%         Row 4 contains the number of the right hand region for each
%         segment. The value 0 corresponds to the exterior of the circle.
%
%   [X,Y]=CIRCLEG(BS,S) returns coordinates of boundary points. BS
%   specifies the boundary segments and S the corresponding parameter
%   values. S may be a scalar, vector, or matrix. BS must be either a
%   scalar or have the same dimensions as S. X and Y have the same
%   dimensions as S.
%
% See also PDEGEOM.
% Copyright 1994-2012 The MathWorks, Inc.
global R;

nbs=4; % Number of boundary segments

d=[0 0 0 0   % Start parameter value
   1 1 1 1   % End parameter value
   1 1 1 1   % Left hand region
   0 0 0 0]; % Right hand region

start_angles = [0;pi/2;pi;3*pi/2];
switch nargin
    case 0
        x=nbs; 
    case 1    
        x=d(:,bs);
    case 2
        if isscalar(bs)
            bs = repmat(bs,size(s));
        end
        x = zeros(size(s));
        y = zeros(size(s));
        for i = 1:numel(s)
           segment = bs(i);
           offset = pi/2*s(i);
           angle = start_angles(segment)+offset;
           x(i) = R*cos(angle);
           y(i) = R*sin(angle);            
        end
end