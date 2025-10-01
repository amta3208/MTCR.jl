function Pa = torr2Pa(torr)
    % Convert units from torr (torr) to Pascals (Pa)
    if nargin == 1
        Pa = torr * 133.322;
    else
        Pa = 133.322;
    end
end

