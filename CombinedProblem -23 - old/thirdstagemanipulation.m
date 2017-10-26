function [phivals,zetavals,x,y,z,Payload] = thirdstagemanipulation(file)
% This function sorts the payload to orbit data files into a form that is
% easily interpolated in MATLAB

ThirdStageData = dlmread(file);
ThirdStageData = sortrows(ThirdStageData);

phivals = unique(ThirdStageData(:,1));
zetavals = unique(ThirdStageData(:,2));
Vvals = unique(ThirdStageData(:,3));
thetavals = unique(ThirdStageData(:,4));
vvals = unique(ThirdStageData(:,5));

x = Vvals;
y = thetavals;
z = vvals;

n = 1;
for a = 1:length(phivals)
    for b = 1:length(zetavals)
        for i = 1:length(Vvals)
            for j = 1:length(thetavals)
                for k = 1:length(vvals)

                    if ThirdStageData(n,1) == phivals(a) && ThirdStageData(n,2) == zetavals(b) && ThirdStageData(n,3) == x(i) && ThirdStageData(n,4) == y(j) && ThirdStageData(n,5) == z(k) % allows for null values in certain cells
                    Payload(a,b,j,i,k) = ThirdStageData(n,6);
                    n = n+1; 
                    else
                    Payload(a,b,j,i,k) = NaN;
                    end
                    
                end
            end
        end
    end
end

