function a = getInteriorPoint(gamma, gammap)
%GETINTERIORPOINT   Compute a point in the interior of a curve.
%   A = GETINTERIORPOINT(gamma, gammap) computes a point A in the interior of
%   the curve gamma, which is negatively oriented.  gamma and gammap are
%   discretizations of gamma(t) and gamma'(t) with equispaced points in
%   [0, 2*pi[.
% 
% Author: Olivier Sète, v 0.1, 14 Septembre 2016.

% Try the center of mass of gamma:
a = mean(gamma);

% Check that a is in the interior of gamma:
if ( isInteriorPoint(a, gamma) )
    return
end

% Check orientation of gamma:
signeadArea = 2*pi*mean(real(gamma) .* imag(gammap));
if ( signeadArea >= 0 )
    error('Orientation of the curve is positive, but must be negative.')
end


% Get an interior point by moving inwards from a boundary point.
% We try different boundary points and (up to) two scalings:
len = length(gamma);
if ( len < 4 )
    startIndices = 1;
else
    startIndices = [floor(len/4), floor(len/2), floor(len*3/4), len];
end

% Find eligible points:
a = [];
for jj = 1:length(startIndices)
    startIndex = startIndices(jj);
    scale = max(abs(gamma(startIndex)), abs(gammap(startIndex)))*(0.1:0.1:1).';
    % Get points interior to gamma:
    anew = inwardPoints(gamma, gammap, startIndex, scale);
    
    % If no points were found, try once more with smaller step-size:
    if ( isempty(anew) )
        scale = scale/10;
        anew = inwardPoints(gamma, gammap, startIndex, scale);
    end
    a = [a; anew];
end

if ( isempty(a) )
    error('No interior points found by this method.')
else
    % Compute the distance to the boundary:
    for jj = 1:length(a)
        distance(jj,1) = min(abs(gamma - a(jj,1)));
    end
    % Pick the point that is furthest from the boundary:
    [~, I] = max(distance);
    a = a(I);
end

end


function out = isInteriorPoint(a, gamma)
%ISINTERIORPOINT   Check if a is in the interior of gamma.

gamma = [gamma; gamma(1)];                      % Last point = first point.
argchange = unwrap(angle(gamma-a));             % Continuous arg(f) along gamma.
winding = (argchange(end)-argchange(1))/(2*pi); % Change of Argument.
out = ( abs(winding) > 0.5 );
end


function a = inwardPoints(gamma, gammap, startIndex, scale)
%INWARDPOINTS   COnstruct points by moving inwards from a boundary point.

% Since gamma is negatively oriented, the inner normal at gamma(t) is 
% gamma(t) - 1i * gammap(t).

a = gamma(startIndex) - 1i*gammap(startIndex) * scale;

% Check which points are in the interior of gamma and pick these:
isIntPt = zeros(size(a));
for j = 1:length(a)
    isIntPt(j,1) = isInteriorPoint(a(j), gamma);
end
a = a(isIntPt == 1);

end
