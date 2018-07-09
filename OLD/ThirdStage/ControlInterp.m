function alphainterp = AlphaInterp(t,AlphaList,tsearch)

% alphainterp = interp1(t,AlphaList,tsearch,'previous');
alphainterp = spline(t,AlphaList,tsearch);
end