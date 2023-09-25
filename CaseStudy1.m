load COVIDbyCounty.mat

%3.2

divisionCol = COVIDbycounty.CNTY_CENSUS("DIVISION");

idx = kmeans(divisionCol, 9, Replicates = 10); %clustering the divisions into nine groups, 
             %CNTY_CENSUS(:,3)

figure()
%plot(dates, idx)
plot(divisionCol(idx==1,1),divisionCol(idx==1,2),'r.')


% silhouette(CNTY_CENSUS', idx);
% 
% [idx,C,sumd,D] = kmeans(CNTY_CENSUS', 9, Replicates = 10);
% figure()
% silhouetteData = silhouette(CNTY_CENSUS', idx);
% title("Population");

% 3.3
%Nearest neighbor
%Idx = knnsearch(X,Y)