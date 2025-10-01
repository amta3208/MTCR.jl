function ms = s2ms(s)
    % Convert units from seconds (s) to milliseconds (ms)
    if nargin == 1
        ms = s * 1e3;
    else
        ms = 1e3;
    end
end

