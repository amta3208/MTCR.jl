function name = spnm2fieldname(spnm)
% Convert a species name (ex: Ar+) to a valid fieldname (ex: Arp)
  name = strrep(spnm,'-','n');
  name = strrep(name,'+','p');
  name = strrep(name,'*','ALL');
end

