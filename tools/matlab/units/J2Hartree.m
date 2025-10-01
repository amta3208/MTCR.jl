function hartree = J2Hartree(J)
    % Convert units from Joules (J) to Hartrees (Hartree)
    if nargin == 1
        hartree = J ./ 4.3597482E-18;
    else 
        hartree = 1 / 4.3597482E-18;
    end
end

