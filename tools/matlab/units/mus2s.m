function s = mus2s(mus)
    % Convert units from microseconds (s) to seconds (ms)
    if nargin == 1
        s = mus * 1e-6;
    else
        s = 1e-6;
    end
end

