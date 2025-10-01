function s = ms2s(ms)
    % Convert units from milliseconds (ms) to seconds (s)
    if nargin == 1
        s = ms * 1e-3;
    else
        s = 1e-3;
    end
end

