function hb = hbar(unitSystem)
% hbar(unitSystem) returns the value of the reduced Planck's constant in the given
% unit system, either ''CGS'' or ''SI''.
  switch unitSystem
    case 'SI'
      hb = hbarSI;
    case 'CGS'
      hb = hbarCGS;
    otherwise
      error('Unrecognized unit system string. Only accepts ''CGS'' or ''SI''.')
  end
  
end

