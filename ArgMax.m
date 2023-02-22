function [maxIndexes] = ArgMax(s, dim)
    [Values, Indexes] = max(s, [], dim);
    maxIndexes = (s ./ Values) == 1;
end
