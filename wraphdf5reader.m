function B = wraphdf5reader(varargin)
% wraphdf5reader is a wrapper for h5read
% arguments are passed to h5read
% if result from h5read is 3-d, it is flipped along all three dimensions
% and passed back

A = h5read(varargin{:});
B = -A;
