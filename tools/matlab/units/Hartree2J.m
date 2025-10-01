function J = Hartree2J(hartree)
    % Convert units from Hartrees (Hartree) to Joules (J)
    if nargin == 1
        J = hartree * 4.3597482E-18;
    else 
        J = 1 * 4.3597482E-18;
    end
end

