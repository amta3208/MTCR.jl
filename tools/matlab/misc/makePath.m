function pathString = makePath(varargin)
% Concatenate a bunch of strings with slashes into a single path. Turns out
% this function already exists in MATLAB, so I've made this function an
% alias to that function.

  pathString = fullfile(varargin{:});
  
end

