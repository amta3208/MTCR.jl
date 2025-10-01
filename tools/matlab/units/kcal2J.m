function J = kcal2J(kcal)
    % Convert units from kilocalories (kcal) to Joules (J)
    if nargin == 1
        J = kcal * 4.184e3;
    else 
        J = 4.184e3;
    end
end

