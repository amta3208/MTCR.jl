function J = erg2J(erg)
    % Convert units from ergs (erg) to Joules (J)
    if nargin == 1
        J = erg * 1e-7; % J
    else
        J =  1  * 1e-7; % J
    end
end
