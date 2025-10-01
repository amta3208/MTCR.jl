function h = planck(unitSystem)
% planck(unitSystem) returns the value of Planck's constant in the given
% unit system, either ''CGS'' or ''SI''.
  switch unitSystem
    case 'SI'
      h = planckSI;
    case 'CGS'
      h = planckCGS;
    otherwise
      error('Unrecognized unit system string. Only accepts ''CGS'' or ''SI''.')
  end
  
end

