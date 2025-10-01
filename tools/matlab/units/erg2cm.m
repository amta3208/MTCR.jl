function cm = erg2cm(erg)
    % Convert units from erg to inverse cm
    if nargin == 0
      cm =  1  / (cm2eV * eV2erg);
    else
      cm = erg / (cm2eV * eV2erg);
    end
end

