function kB = boltz(unitSystem)
% boltz(unitSystem) returns the value of Boltzmann's constant in the given
% unit system, either ''CGS'' or ''SI''.
  switch unitSystem
    case 'SI'
      kB = boltzSI;
    case 'CGS'
      kB = boltzCGS;
    otherwise
      error('Unrecognized unit system string. Only accepts ''CGS'' or ''SI''.')
  end
  
end

