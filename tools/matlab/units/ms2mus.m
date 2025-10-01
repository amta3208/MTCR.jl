function mus = ms2mus(ms)
    % Convert units from milliseconds (ms) to microseconds (mus)
    if nargin == 1
        mus = ms * 1e3;
    else
        mus = 1e3;
    end
end

