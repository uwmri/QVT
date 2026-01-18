function B = wraphdf5reader(varargin)
% wraphdf5reader is a wrapper for h5read
% arguments are passed to h5read
% if result from h5read is 3-d, it is flipped along all three dimensions
% and passed back

A = h5read(varargin{:});
d = size(A);
if length(d) == 3
    B = flip(flip(flip(A,1),2),3);
else
    B = A;
    fprintf('WARNING: h5read returned a matrix of dimension %d\n',length(d));
    fprintf('    matrix was passed back unchanged');
end

