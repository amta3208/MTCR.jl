function mus = s2mus(s)
    % Convert units from seconds (s) to microseconds (ms)
    if nargin == 1
        mus = s * 1e6;
    else
        mus = 1e6;
    end
end

