function J = cm2J(cm)
    % Convert units from inverse cm to Joules (J)
    if nargin == 1
        J = eV2J(cm2eV(cm));
    else
        J = eV2J(cm2eV(1));
    end
end

