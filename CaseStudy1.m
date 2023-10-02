load COVIDbyCounty.mat

%% Lab 1

    %% Q1: IDENTIFY MOST POPULATED COUNTY IN EACH DIVISION

S = [];  %creating an empty array to add our max pop val's and corresponding names to.
indexMaxpop = zeros(9,156);
figure('Name', 'Covid Trends for the most populated county in each divisions');
hold on
for div = 1:9      %for numbers one through nine 
    idx = CNTY_CENSUS(:,3)==div; %all the places where division = i, put a 1  
    idxvalues = table2array(idx);  %turns idx into actual 1's and 0's
   
    divpop = CNTY_CENSUS(idxvalues==1, 6);    %wherever idx = 1, grab the 6th column of that row in CNTY CENUS (the pop data)
    divnames =  CNTY_CENSUS(idxvalues==1, 5);  %same thing, except it's col 5 (the name)
    
    maxdivpop = max(divpop);   %get the max pop out of all the counties in the division
    maxval= table2array(maxdivpop);  %turn it into a number so it can be used in an == way.
 
    maxnameindex = CNTY_CENSUS(CNTY_CENSUS.POPESTIMATE2021==maxval,5);  %when the row value = maxval, grab col 5 (the name)
    maxnameindexval = table2array(maxnameindex);

    a  = find(CNTY_CENSUS.POPESTIMATE2021==maxval);

    data = CNTY_COVID(a,:); 
    indexMaxpop(div,:) = data;

    plot(dates, data)  %plotting the dates by the max pop per division trends (cases)

    hold on 
    T = table(maxval,maxnameindexval); %table with the maxval in col 1 and the corresponding name in col 2

    S = vertcat(S,T); % concatenate into a table, it should be 9x2 by the end of the for loop (one row per div)

end

    %% Q2: SEE IF TRAJECTORIES ARE LINEARLY INDEPENDENT BY FINDING ANGLE BETWEEN TRAJECTORIES

anglevals = zeros(9,9);

for i = 1:9
    for j = (i+1):9
        dotProd  = dot(indexMaxpop(i,:),indexMaxpop(j,:));
        angle = acosd(dotProd/(norm(indexMaxpop(i,:))*norm(indexMaxpop(j,:))));
       % if angle~=0
           % fprintf('the angle is not zero--linearly independent');
       % end
        anglevals(i,j) = angle;
    end

end

    %% Q3: NORMALIZING THE DATA

normalizeDataTable = zeros(9,156);
  
for i = 1:9
     normalizeData = indexMaxpop(i,:)/norm(indexMaxpop(i,:));
     normalizeDataTable(i,:) = normalizeData;
end

d1 = normalizeDataTable(1,:);
d2 = normalizeDataTable(2,:);
d3 = normalizeDataTable(3,:);
d4 = normalizeDataTable(4,:);
d5 = normalizeDataTable(5,:);
d6 = normalizeDataTable(6,:);
d7 = normalizeDataTable(7,:);
d8 = normalizeDataTable(8,:);
d9 = normalizeDataTable(9,:);
      
    %% Q4: FIND CASE DATA FOR STL city (row 135), PLUG INTO FORMULA TO OBTAIN VECTOR

STLdataindex  = find(CNTY_CENSUS.CTYNAME=="St. Louis city"); %locating row number of stl in cnty census
c = CNTY_COVID(STLdataindex,:);  %for that row in cnty covid, get all the columns  (this will get the entire row of covid data values at the STL row)

%  ri = c-(c.*di) * di %this is the formula 

r1 = c - (dot(c,d1)*d1);  %finding each r
r2 = c - (dot(c,d2)*d2);
r3 = c - (dot(c,d3)*d3);
r4 = c - (dot(c,d4)*d4);
r5 = c - (dot(c,d5)*d5);
r6 = c - (dot(c,d6)*d6);
r7 = c - (dot(c,d7)*d7);
r8 = c - (dot(c,d8)*d8);
r9 = c - (dot(c,d9)*d9);

magr1 = norm(r1);  %finding the magnitude of each r
magr2 = norm(r2);
magr3 = norm(r3);
magr4 = norm(r4);
magr5 = norm(r5);
magr6 = norm(r6);
magr7 = norm(r7);
magr8 = norm(r8);
magr9 = norm(r9);

magnitudesOfR = table(magr1,magr2,magr3,magr4,magr5,magr6,magr7,magr8, magr9);  %table with all magnitude vals

%% Q5

% Interpret ri and its norm. What do these describe?
% ri shows the deviation between STL data and the biggest counties per
% division's covid case data. In the formula, we project the normalized
% vector, d, onto c, the STL covid case data. Then we multiply it by d to
% make this into a basis expansion, and subtract it from c. In total, this
% gives us the deviation between the biggest counties and STL's case
% data.The NORM of ri gives us a numerical difference between STL and the
% big counties/the magnitude of the ri vector.

% What might they indicate about St. Louis City relative to the 9 census divisions?
%This would mean that STL's data is similiar to some divisions more than
%others. Ex: STL is more different from div2 than it is from div1.

%% 3.2

%% GENERATING THE RANDOM CENTROIDS

c1 = randn(1, 2); %this returns a 1x2 table of random values--this acts as our coordinates.
c2 = randn(1, 2);
c3 = randn(1, 2);
c4 = randn(1, 2);
c5 = randn(1, 2);
c6 = randn(1, 2);
c7 = randn(1, 2);
c8 = randn(1, 2);
c9 = randn(1, 2);

% ALTERNATIVE METHOD
% for (loop 9 times)
%     generate a random number, assign to x,
%     generate a random number, assign to y 
%     Then make that a coordinate vector
%     Also add it to a larger vector
%     Plot on graph

% OR
% x = math.Random(9) %coordinates of the centroids
% y = math.Random(9) 
% 
% plot (x,y)


%% CALCULATING DISTANCE BETWEEN DATA POINT AND CENTROIDS 
%calculate the distance from data point to each [randomly generated] centroid
 %WHICH MEANS   
        % make each randomly generated centroid a vector, within the larger matrix of centroids
        %FOR each data point (each centroid), calculate the distance between

% D = pdist2(c1,datapoint,'euclidean');
% disp(D)


%ALTERNATIVE METHOD -- NESTED FOR LOOP?
% for (each data point)
% 
%     for each centroid


%xmin distance from centroid 

%Assign data point to that centroid

%Recalculate the cluster's centroid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Building kmeans 

choppedData = CNTY_COVID(1:180,:); % only use the first 200 columns for kmeans.
hiddenData = CNTY_COVID(181:225, :);

figure()
plot(choppedData')

figure()
[idx,C,D] = kmeans(choppedData, 9, Replicates = 200);
plot(C')

%Pacific, Mountain, West North Central, West South Central, East North
%Central, East South Central, Middle Atlantic, South Atlantic, New England

%% How do we classify the data
%find the centroid, get the nearest neighbor, wherever that nearest neighbor is, use that as the label for the centroid.
% instead of nearest neighbor, manually find the closest data point to the
% centroid.



%find the distances between the centroids and each point in the cluster.
for idx = 1:9   
       distance = norm(choppedData(idx(i)-C(idx)))  %county1 - centroid
       if mindistance > distance 
          mindistance = distance
       end 
end
