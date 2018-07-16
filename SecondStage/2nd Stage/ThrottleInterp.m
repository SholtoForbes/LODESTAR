function alphainterp = ThrottleInterp(t,AlphaList,tsearch)

alphainterp = interp1(t,AlphaList,tsearch,'next');
% alphainterp = spline(t,AlphaList,tsearch);
end