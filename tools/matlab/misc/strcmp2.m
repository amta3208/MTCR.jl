function TF = strcmp2(str1,str2,option,ids)
% Compare subset of str1 with str2, with multiple checks to avoid out of
% bounds errors for indexing. 
% Options inclide: 'start' -- check if str1 starts with str2
%                  'end'   -- check if str1 ends with str2
%                  'ids'   -- check if str1(ids) == str2

  % --------------------------------------------------
  % Revert to regular strcmp if no option is specified
  % --------------------------------------------------
  if nargin == 2
    TF = strcmp(str1,str2);
    
  % ----------------------
  % Switch special options
  % ----------------------
  elseif any(nargin == [3,4])
    % Required length check for all options 
    ls1 = length(str1);
    ls2 = length(str2);
    if ls1 < ls2
      TF = 0; return;
    end
    switch option
      % Check if str1 starts with str2
      case 'start'
        TF = strcmp(str1(1:ls2),str2);
      % Check if str1 ends with str2
      case 'end'
        TF = strcmp(str1(end-(ls2-1):end),str2);
      % Check if str1(ids) == str2
      case 'ids'
        if (nargin == 4)
          if (ls1 < max(ids))
            TF = 0;
          else
            TF = strcmp(str1(ids),str2);
          end
        else
          error('Must specify array of indices as fourth argument when option == ''ids''');
        end
      % Error on unrecognized options
      otherwise
        error('Accepted option types to strcmp2 include: ''start'', ''end'', and ''ids''.');
    end
  else
    error('Must pass 2, 3, or 4 arguments to strcmp2. Read function description for more information.');
  end
end

