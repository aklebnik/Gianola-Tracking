%Function to determine the nearest neighbors for each particle. Along with
%the index of each relevant particle. Aiming to implement either Tim Still
%or Eric Weeks polydispersity script.

load finalvalues.mat

%%
close all

% out(:,1) = out (:,1) - min (out(:,1)); %move particles to origin, due to mask of cntrd.m
% out(:,2) = out (:,2) - min (out(:,2));

out (:,1:2) = out (:,1:2)/9; %size correction, x,y
out (:,3) = out(:,3)/5; %size correction, z
save ('adjustedvalues', 'out')

sz = length (out);
dist = zeros(sz); %create preloop

cnt = out; %create temporary storage of particle centers
cnt (:,4:end) = []; %get rid of superflous data, to run faster


for i = 1:sz
    %for i = 1:10
    temppeak= cnt(i,:); %get the specific point
    x1 = bsxfun(@minus,cnt,temppeak); %subtract the coordinates from the matrix of point coordinates
    x1 = x1.*x1;   %square the values
    dist(:,i) = sqrt (x1 (:,1) + x1 (:,2) + x1 (:,3)); %find distance
end


[pLoc, pInd] = sort (dist); %sort all the values, save the values and the indices
pLoc = pLoc (2:6,:); %for first seven nearest neighbors
pInd = pInd (2:6,:); %find unique particle number for the seven nearest neighbors

pAvg = mean (pLoc); %find the average particle distance
hist (pAvg,20)
title ('Average distance in um to nearest 7 neighbors')

pMin = pLoc(1,:); %find the smallest particle distance
min (pMin)
figure
hist (pMin,20)
title ('Minimum distance in um to nearest neighbor')

%%

abar = 2.7/2; %microns

x = [0:0.01:3];
y  = normpdf (x, 1.35, 0.135);
y = y/100;

pLoc1 = pLoc(pLoc<=4*abar);
delta = mean (pLoc1) - 2*abar

%%
pRadii = ones (1, size (pLoc,2)) .*abar;
temp = 0;
num =1;

%numel(pLoc(pLoc<2*abar))

for t= 1:15 %iterations
    for i = 1:size (pLoc,2) %number particles
        for j = 1:size (pLoc,1) %number of nearest neighbors
            tempInd = pInd(j,i);       
           overlap(j,i)=  pLoc (j,i) - pRadii(i) - pRadii(tempInd);
        end
    end
    yOverlap = logical(sum(overlap <0));

    pRadii(yOverlap)= (pRadii(yOverlap))*0.95;
    nOverlap = numel (overlap (overlap <0))
%     if nOverlap ==0
%         pRadii = pRadii *1.01;
%     end
end
hist (pRadii,30)

   
end


%%
B =reshape (pLoc,size (pLoc,1)*size(pLoc,2),1);
figure
hist (B,30)
title ('Distances for nearest 7 neighbors')
