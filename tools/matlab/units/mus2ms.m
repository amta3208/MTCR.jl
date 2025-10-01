function ms = mus2ms(mus)
    % Convert units from microseconds (mus) to milliseconds (ms)
    if nargin == 1
        ms = mus * 1e-3;
    else
        ms = 1e-3;
    end
end

