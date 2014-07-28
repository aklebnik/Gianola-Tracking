% WRITTEN BY SASHA KLEBNIKOV, 7/1/14
% Based on implementing Ye Xu's cntrd3.m, and trackS.m
% working on 7/2/14

filelist_generator_Sasha
%%
load filenamelist.mat
Params = [8, 31, 2, 1,1];% for Dilute Yodh

field1 = 'mem';  value1 = (0);
field2 = 'goodenough';  value2 = (7);
field3 = 'dim';  value3 = (2);
field4 = 'quiet';  value4 = (0);
field5 = 'mxdisp';  value5 = (10);
field6 = 'method';  value6 = ('normal');

param = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);

th = Params (1);
szx = Params (2);
szz = 20;

im = filenamelist;
z1 = 1;
z2 = length (filenamelist);
imstack = 0;
mxdisp = 10;
maxdisp = mxdisp;


tic
pks = [];
for i = z1:z2
    img=imread(im(i,:));
    img = double(img);
    im_filter=bpass(img,2,szx); %form: img data, minimum size, maximum size
    pka=pkfnd(im_filter,th,szx);
    pkb=cntrd(im_filter,pka,szx+2);
    if ~isempty(pkb)
        pks = [pks; [pkb(:,1:4) i*ones(size(pkb,1),1)]];
    end
end
toc

%track particles vertically
trks = trackS(pks,mxdisp,param);


%% Trying to get particle centers
np = max(trks(:,6)) %number of total tracks
temp = 1;

new_trks = zeros(1,6);

for i = 1:np
    
    clear Loc
    clear trki
    clear vec
    
    trki = trks(trks(:,6)==i,:); % select the track we are fitting stuff too
    Slices = size(trki,1); %how many slices is the data?
    vec = [trki(:,5) trki(:,3)]; %create a subset, consisting of z slice number and brightness
    
    if Slices>15
       
        output = peakfit(vec,0,0,round(Slices/18),1,0,10,0,1); % I found this works best for fitting the peaks, still not the best we 
        % can do
        
        %get indices
        Loc = find (trks (:,6)==temp); %where in the massive trks matrix are we?
        Ind = []; %clear previous operations
        
    else
        
        output = peakfit(vec,0,0,1,6);
        
    end
    
    new_trks_temp = zeros(size(output,1),6); % x, y, z, brightness, rg, particle index at peak center in z
    
    particles = size(new_trks,1);
    
    for k = 1:size(output,1)
        
        z_temp = output(k,2); % z-center
        index = round(z_temp); % round to neared slice #
        temp_val = abs(trki(:,5)-index); 
        temp_vec = trki(:,5);   
        [num idx] = min(temp_val); % find closes z-slice to fit center         
        slice_no = idx;        
        new_trks_temp(k,:) = [trki(slice_no,1), trki(slice_no,2), z_temp, trki(slice_no,3), trki(slice_no,4), particles+k];
        % create new vector with x, y, z, brightness, rg, and particle #
        
    end
    
    new_trks = [new_trks;  new_trks_temp]; % make new trks
        
end

new_trks(1,:) = [];

new_trks(:,6) = new_trks(:,6)-ones(size(new_trks,1),1);

[val, binary] = find(new_trks(:,3) > z2); % get rid of fit centers larger than # of slices
new_trks(val,:) = []; % get rid of fit centers larger than # of slices

[val] = find(new_trks(:,3) < 0); % get rid of fit centers less than 0
new_trks(val,:) = []; % get rid of fit centers less than 0

new_trks(:,1) = new_trks(:,1)/9.371; % scale pixels
new_trks(:,2) = new_trks(:,2)/9.371; % scale pixels
new_trks(:,3) = new_trks(:,3)/5; % scale pixels

save finalvalues.mat new_trks

%% Old code, fitting centroids
np = max(trks(:,6)) %number of total particles

out = zeros(np,4);


%find the 3D centroids based on their brightness in each slice
for i = 1:np
    trki = trks(trks(:,6)==i,:);
    out(i,1) = sum(trki(:,1).*trki(:,3))/sum(trki(:,3));
    out(i,2) = sum(trki(:,2).*trki(:,3))/sum(trki(:,3));
    out(i,3) = sum(trki(:,4).*trki(:,3))/sum(trki(:,3));
    out(i,4) = sum(trki(:,3));
    
end

    
    save ('finalvalues','out')
    
    
    %%
    %plot stuff in a 3d area
    %load adjustedvalues.mat
    figure
    hold on
    grid on
    [x,y,z] = sphere;
%     x = x./9.31; %based on dim
%     y = y./9.31; %change these to change size of spheres
%     z = z./5 %20 slides thick, gives smear, not sphere
%     
    
    x1 = new_trks(:,1); %create temporary variables
    y1 = new_trks(:,2);
    z1 = new_trks(:,3);
    
    for i = 1:length (x1)
        surf(x-x1(i), y-y1(i), z-z1(i))
    end
    axis equal
    
    clear x; clear y, clear z, clear x1, clear y1, clear z1
    
    
    
    
