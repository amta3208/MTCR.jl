function J = cal2J(cal)
    % Convert units from calories (cal) to Joules (J)
    if nargin == 1
        J = cal * 4.184;
    else 
        J = 4.184;
    end
end

