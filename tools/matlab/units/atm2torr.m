function torr = atm2torr(atm)
    % Convert units from atmospheres (atm) to torr (torr)
    if nargin == 1
        torr = atm * 760;
    else 
        torr = 760;
    end
end

