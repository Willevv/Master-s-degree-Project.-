function long = Long_spanlengths(x)
    sl = x(43:end);
    inverted = 1./sl;
    long = sum(inverted);
end
