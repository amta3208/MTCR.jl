function Pa = atm2Pa(atm)
    % Convert units from atmospheres (atm) to Pascals (Pa)
    if (nargin == 1)
        Pa = atm * 101325;
    else
        Pa = 101325;
    end
end

