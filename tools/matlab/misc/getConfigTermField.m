function [config, term, field] = getConfigTermField(name,energy)
% Get a state name from a species name and energy (eV). If no name is
% found, the name is simply given as "E = <energy> eV"

  % Default entries if none of the following tests are a match
  config = '';
  term = defaultName; 
  field = defaultField;

  % Go through the species and their "important" states 
  switch name
    case 'O'
      if isEnergy(0.0097)
        assignCTF('2p^4','^3P','3P');
      elseif isEnergy(1.9674)
        assignCTF('2p^4','^1D','1D');
      elseif isEnergy(4.1898)
        assignCTF('2p^4','^1S','1S');
      elseif isEnergy(9.1462)
        assignCTF('2p^33s','^5S^o','5So');
      elseif isEnergy(9.5215)
        assignCTF('2p^33s','^3S^o','3So');
      end
    case 'O2'
      if isEnergy(0.0000)
        assignCTF('','X^3 \Sigma_g^-','X');
      elseif isEnergy(0.9817)
        assignCTF('','a^1 \Delta_g','a');
      elseif isEnergy(1.6360)
        assignCTF('','b^1 \Sigma_g^+','b');
      end
    case 'Ar'
      if isEnergy(0.0000)
        assignCTF('3p^6','^1S','1S');
      elseif isEnergy(11.5485)
        assignCTF('(^2P^o_{3/2})4s','^2[3/2]^o','4s_1p5');
      elseif isEnergy(11.6238)
        assignCTF('(^2P^o_{1/2})4s','^2[1/2]^o','4s_0p5');
      elseif isEnergy(11.7233)
        assignCTF('(^2P^o_{3/2})4p','^2[1/2]','4p_0p5');
      elseif isEnergy(11.8282)
        assignCTF('(^2P^o_{3/2})4p','^2[5/2]','4p_2p5');
      end
    case 'Ar+'
      if isEnergy(0.000)
        assignCTF('3p^5','^2P^o_{3/2}','j_1p5');
      elseif isEnergy(0.1775)
        assignCTF('3p^5','^2P^o_{1/2}','j_0p5');
      end
  end
  
  function TF = isEnergy(testEnergy)
    TF = round(energy,4) == round(testEnergy,4);
  end
      
  function exnm = defaultName
    exnm = sprintf('E = %1.4f eV', energy);
  end

  function fnm = defaultField
    fnm = strrep(sprintf('%s_%1.4f_eV',name,energy),'.','_');
  end

  function fnm = namedField(str)
    fnm = sprintf('%s_%s',name,str);
  end

  function assignCTF(cnm,tnm,fnm)
    if nargin == 2
      config = cnm;
      term = tnm;
    elseif nargin == 3
      config = cnm;
      term = tnm;
      field = namedField(fnm);
    else
      error('assignCTF expects 2 or 3 arguments.');
    end
  end
end

