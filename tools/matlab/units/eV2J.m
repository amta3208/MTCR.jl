function J = eV2J(eV)
    % Convert units from electron volts (eV) to Joules (J)
    if nargin == 1
        J = eV * 1.602177 * 10^(-19);
    else 
        J = 1.602177 * 10^(-19);
    end
end

