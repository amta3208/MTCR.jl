function torr = Pa2torr(Pa)
    % Convert units from Pascals (Pa) to torr (torr)
    if nargin == 1
        torr = Pa * 0.00750062;
    else 
        torr = 0.00750062;
    end
end

