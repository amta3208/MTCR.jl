function atm = Pa2atm(Pa)
    % Convert units from Pascals (Pa) to atmospheres (atm)
    if nargin == 1
        atm = Pa / 101325;
    else
        atm = 1 / 101325;
    end
end

