function eV = cm2eV(cm)
    % Convert units from inverse cm to electron volts (eV)
    if nargin == 1
        eV = cm / 8065.75;
    else
        eV = 1 / 8065.75;
    end
    
end

