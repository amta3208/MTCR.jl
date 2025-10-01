function cal = J2cal(J)
    % Convert units from Joules (J) to calories (kcal)
    if nargin == 1
        cal = J ./ 4.184;
    else 
        cal = 1 / 4.184;
    end
end

