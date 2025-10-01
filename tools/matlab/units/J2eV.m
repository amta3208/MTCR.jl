function eV = J2eV(J)
    % Convert units from Joules (J) to electron volts (eV)
    if nargin == 1
        eV = J / (1.602177 * 10^(-19));
    else 
        eV = 1 / (1.602177 * 10^(-19));
    end
end

