function atm = torr2atm(torr)
    % Convert units from torr (torr) to atmospheres (atm)
    if nargin == 1
        atm = torr / 760;
    else 
        atm = 1 / 760;
    end
end

