function erg = eV2erg(eV)
    % Convert units from electron volts (eV) to ergs (erg)
    if nargin == 1
        erg = eV / 6.242e+11; % erg
    else
        erg = 1 / 6.242e+11;
    end
end

