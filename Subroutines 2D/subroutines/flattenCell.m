function flatArray = flattenCell(nestedArray)
  error(nargchk(1,1,nargin));
  if ~iscell(nestedArray),
    error('Must be a cell array.');
  end
  flatArray{1} = [];
  for i=1:numel(nestedArray),
    if iscell(nestedArray{i}),
      y = flattenCell(nestedArray{i});
      [flatArray{end+1:end+length(y)}] = deal(y{:});
    else
      flatArray{end+1} = nestedArray{i};
  end
end
flatArray(1) = [];