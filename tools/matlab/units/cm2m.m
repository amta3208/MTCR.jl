function m = cm2m(cm)
    % Convert units from cm to m
    if nargin == 1
        m = cm/100;
    else
        m = 1/100;
    end
end

