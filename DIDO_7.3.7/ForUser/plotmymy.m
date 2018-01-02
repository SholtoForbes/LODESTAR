function plotmymy(x, y1, y2, z1, s1, z2, s2)
%PLOTMYMY Plots 4 graphs on the same figure with two ordinate scales
%	PLOTMYY(X, Y1, Y2, Z1, S1, Z2, S2)
%	X is the abscissa for all plots.
%	Y1 and Y2 are the ordinates Y1(X), Y2(X)
%	Z1 and Z2 are ordinates in the range of Y1 and Y2 respectively, and
%	S1 and S2 are the stings (e.g. '*') to be used for plotting Z1 and Z2 
%	respectively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by I. Michael Ross.  Part of the DIDO package
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = plotyy(x, y1, x, y2);
axes(ax(1))
hold on; plot(x, z1, s1)
axes(ax(2))
hold on; plot(x, z2, s2)
hold off;
