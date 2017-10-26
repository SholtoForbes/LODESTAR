% visualising third stage cost results
% clear all
ThirdStageData = sortrows(dlmread('thirdstageforvis.dat'));
Atmosphere = dlmread('atmosphere.txt');
% map = [0, 0, 0
%     0, 0, 0
%     .0, 0, 0
%     .2, 0, 0
%     .6, 0, 0.0
%     1.0, 0, 0];
% colormap(map)




colormap jet
scatter3(ThirdStageData(:,3),ThirdStageData(:,4),ThirdStageData(:,5),30,ThirdStageData(:,6), 'filled')

ylabel('Trajectory Angle (deg)')
xlabel('Velocity (m/s)')
zlabel('Altitude (km)')
c=colorbar('northoutside');
ylabel(c,'Payload Mass, kg')



[VGrid,thetaGrid,vGrid] = ndgrid(33000:500:34500,[deg2rad(2) :deg2rad(1): deg2rad(5)],[2700:25:2875]);

PayloadDataInterp = scatteredInterpolant(ThirdStageData(:,3),ThirdStageData(:,4),ThirdStageData(:,5),ThirdStageData(:,6));
PayloadData = PayloadDataInterp(VGrid,thetaGrid,vGrid);
global PayloadGrid
PayloadGrid = griddedInterpolant(VGrid,thetaGrid,vGrid,PayloadData,'spline','linear');


% [meshAlt,meshAngle,meshVel] = meshgrid(unique(ThirdStageData(:,3)),unique(ThirdStageData(:,4)),unique(ThirdStageData(:,5)));
% meshPayload = permute(reshape(ThirdStageData(:,6),[length(unique(ThirdStageData(:,5))),length(unique(ThirdStageData(:,4))),length(unique(ThirdStageData(:,3)))]),[2 3 1]);
 figure(200)

colormap(gray)
C = contourf(rad2deg(thetaGrid(:,:,2)),VGrid(:,:,2)/1000,PayloadData(:,:,1),20,'LineWidth',0.)
% C = contourf(rad2deg(meshAngle(:,:,2)),meshAlt(:,:,2)/1000,meshPayload(:,:,3),12)
% C = contourf(rad2deg(meshAngle),meshAlt/1000,meshPayload,13)

xlabel('Release Angle (deg)')
ylabel('Altitude (km)')
% title('2750 m/s')
% ylim([30 40])
c = colorbar;
c.Label.String = 'Payload (kg)';

% slice(threeDMAT,[1,3,5],[],[])
% colormap hsv


% contourslice(threeDMAT,[1,3,6],[],[],10)


% 
% % getting maximum values
% [max1,ind1] = max(threeDMAT,[],3);
% 
% [max2,ind2] = max(max1,[],2);
% 
% ind  = [ind1(ind2),ind2] % this really doesnt tell me much, as it wants very high altitude all the time
% 
% 
% 
% 
% 
% % plotting density
% 
% 
% figure(4)
% plot(Atmosphere(:,1)/1000,  Atmosphere(:,4)),'Color','k','LineWidth',0.8;
% xlim([20,40])
% 
% xlabel('Altitude (km)')
% ylabel('Density (kg/m^3)')
% 
% 
% figure(5)
% 
% [meshAngle2,meshVel2,meshAlt2] = meshgrid(Angle,Vel,Alt);
% 
% colormap(gray)
% 
% colormap(flipud(colormap))
% hx = slice(meshAngle2,meshVel2,meshAlt2./1000,threeDMAT,[],[2600:100:3000],[])
% 
% % [Xq,Yq,Zq] = meshgrid([0:.5:5],[2600:100:3000],[20000:500:40000]);
% % threeDMATint = interp3(meshAngle2,meshVel2,meshAlt2,threeDMAT,Xq,Yq,Zq);
% % hx = slice(Xq,Yq,Zq/1000,threeDMATint,[],[2600:100:3000],[]);
% 
% caxis([220 340]) 
% c=colorbar('northoutside')
% c.LimitsMode = 'manual'
% c.Limits = [180 320]
% 
% set(hx, 'EdgeColor', 'none');
% 
% % set(hx, 'FaceColor', 'interp'); % use this to smooth
% set(hx, 'FaceColor', 'flat');
% 
% xlabel('Trajectory Angle (deg)','FontSize',17)
% ylabel('Velocity (m/s)','FontSize',17)
% zlabel('Altitude (km)','FontSize',17)
% ylabel(c,'Payload Mass, kg','FontSize',17)
