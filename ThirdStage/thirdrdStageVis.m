% visualising third stage cost results
clear all
ThirdStageData = sortrows(dlmread('thirdstage.dat'));
Atmosphere = dlmread('atmosphere.txt');
% map = [0, 0, 0
%     0, 0, 0
%     .0, 0, 0
%     .2, 0, 0
%     .6, 0, 0.0
%     1.0, 0, 0];
% colormap(map)




colormap jet
scatter3(ThirdStageData(:,1),ThirdStageData(:,2),ThirdStageData(:,3),30,ThirdStageData(:,4), 'filled')

ylabel('Trajectory Angle (deg)')
xlabel('Velocity (m/s)')
zlabel('Altitude (km)')
c=colorbar('northoutside');
ylabel(c,'Payload Mass, kg')



[Alt,Angle,Vel,Payload] = thirdstagemanipulation('thirdstagefine.dat');

[meshAlt,meshAngle,meshVel] = meshgrid(Alt,Angle,Vel);

for i = 1:length(Alt)
    for j = 1:length(Angle)
        for k = 1: length(Vel)
            meshPayload(j,i,k) = interp3(meshAlt,meshAngle,meshVel,Payload,Alt(i),Angle(j),Vel(k));

        end
    end
end

figure(2)
ax1 = subplot(2,2,1);
colormap(gray)
C = contourf(meshAngle(:,:,4),meshAlt(:,:,4)/1000,meshPayload(:,:,10),13)

xlabel('Trajectory Angle (deg)')
ylabel('Altitude (km)')

title('2600 m/s')


ax2 = subplot(2,2,2);
colormap(gray)
C = contourf(meshAngle(:,:,5),meshAlt(:,:,5)/1000,meshPayload2,13)

xlabel('Trajectory Angle (deg)')
ylabel('Altitude (km)')

title('2800 m/s')

ax3 = subplot(2,2,3);
colormap(gray)
C = contourf(meshAngle(:,:,6),meshAlt(:,:,6)/1000,meshPayload3,13)

xlabel('Trajectory Angle (deg)')
ylabel('Altitude (km)')

title('3000 m/s')




ax4 = subplot(2,2,4);
colormap(gray)
C = contourf(meshAngle(:,:,6),meshAlt(:,:,6)/1000,meshPayload3,13)
c=colorbar('southoutside')
c.LimitsMode = 'manual'
c.Limits = [200 320]
xlabel('Trajectory Angle (deg)')
ylabel('Altitude (km)')
ylabel(c,'Payload Mass, kg')

set([ax1 ax2 ax3 ax4],'clim',[200 320]);


% slice(threeDMAT,[1,3,5],[],[])
% colormap hsv


% contourslice(threeDMAT,[1,3,6],[],[],10)



% getting maximum values
[max1,ind1] = max(threeDMAT,[],3);

[max2,ind2] = max(max1,[],2);

ind  = [ind1(ind2),ind2] % this really doesnt tell me much, as it wants very high altitude all the time





% plotting density


figure(4)
plot(Atmosphere(:,1)/1000,  Atmosphere(:,4)),'Color','k','LineWidth',0.8;
xlim([20,40])

xlabel('Altitude (km)')
ylabel('Density (kg/m^3)')


figure(5)

[meshAngle2,meshVel2,meshAlt2] = meshgrid(Angle,Vel,Alt);

colormap(gray)

colormap(flipud(colormap))
hx = slice(meshAngle2,meshVel2,meshAlt2./1000,threeDMAT,[],[2600:100:3000],[])

% [Xq,Yq,Zq] = meshgrid([0:.5:5],[2600:100:3000],[20000:500:40000]);
% threeDMATint = interp3(meshAngle2,meshVel2,meshAlt2,threeDMAT,Xq,Yq,Zq);
% hx = slice(Xq,Yq,Zq/1000,threeDMATint,[],[2600:100:3000],[]);

caxis([220 340]) 
c=colorbar('northoutside')
c.LimitsMode = 'manual'
c.Limits = [180 320]

set(hx, 'EdgeColor', 'none');

% set(hx, 'FaceColor', 'interp'); % use this to smooth
set(hx, 'FaceColor', 'flat');

xlabel('Trajectory Angle (deg)','FontSize',17)
ylabel('Velocity (m/s)','FontSize',17)
zlabel('Altitude (km)','FontSize',17)
ylabel(c,'Payload Mass, kg','FontSize',17)
