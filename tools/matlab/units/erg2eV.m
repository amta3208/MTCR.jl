function eV = erg2eV(erg)
    % Convert units from ergs (erg) to electron volts (eV)
    if nargin == 1
        eV = erg * 6.242e+11; % eV
    else
        eV = 6.242e+11;
    end
end
