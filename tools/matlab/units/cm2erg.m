function erg = cm2erg(cm)
    % Convert units from inverse cm to erg
    if nargin == 0
      erg = cm2eV * eV2erg;
    else
      erg = cm * cm2eV * eV2erg;
    end
end

