function erg = J2erg(J)
    % Convert units from ergs (erg) to Joules (J)
    if nargin == 1
        erg = J / 1e-7; % erg
    else
        erg = 1 / 1e-7; % erg
    end
end
