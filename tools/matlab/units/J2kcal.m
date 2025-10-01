function kcal = J2kcal(J)
    % Convert units from Joules (J) to kilocalories (kcal)
    if nargin == 1
        kcal = J ./ 4.184e3;
    else 
        kcal = 1 / 4.184e3;
    end
end

