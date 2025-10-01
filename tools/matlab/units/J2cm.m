function cm = J2cm(J)
    % Convert units from Joules (J) to inverse cm
    if nargin == 1
        cm = J / cm2J;
    else
        cm = 1 / cm2J;
    end
end

