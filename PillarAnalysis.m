cd ('Z:\USERS\Klebnikov')

load finalvalues.mat
cd ('Z:\USERS\Klebnikov\Matlab Codes\My Code')
out = out(1:5,:);
out1=[];
out2 = [];
%out (:,3) = out (:,3)*20;
for i = 1:30

    
    out1  =[out i*ones(size(out,1),1)]; %put a unique modifier to identify which time step this is
    out2 = [out2; out1(:,1:3)+ rand(length (out),3)*20 out1(:, 4:5)]; %does brownian motion. 
end

%hist(mod(out(:,1),1),50); %Pixel Biasing

%out2 (:,3) = out2(:,3)*6; % scaling factor
xyzs = out2;
clear out1
% CALIBRATE AXIS HERE
% 
% TF1 = out2(:,1) <= 800; %use this to delete extraneous points.
% out2(TF1,:) = [];
%% Plot 3D spheres
hold on
grid on
[x,y,z] = sphere;
x = x.*20; %based on dim
y = y.*20; %change these to change size of spheres
z = z.*14; %20 slides thick, gives smear, not sphere

%surf (x,y,z)
x1 = out2(:,1); %create temporary variables
y1 = out2(:,2);
z1 = out2(:,3);

for i = 1:length (x1)
surf(x-x1(i), y-y1(i), z-z1(i))
end
clear x; clear y, clear z, clear x1, clear y1, clear z1
%% Check all dimensions for sub pixel accuracy/biasing
subplot (1,3,1)
hist(mod (out(:,1),1),20)
title('Pixel Bias (x)')
subplot (1,3,2)
hist(mod (out(:,2),1),20)
title('Pixel Bias (y)')
subplot (1,3,3)
hist(mod (out(:,3),1),20)
title('Pixel Bias (z)')

%% Pair Correlation Function, g(r)
%out (:,3) = out(:,3)*9;
sz = length (out);
dist = zeros(sz); %create preloop
binsize = 5;
maxbin = 200;
dim = [1500,1500,1500]; %arbitrarily set

cnt = out;
cnt (:,4:end) = []; %get rid of superflous data, to run faster

for i = 1:sz
%for i = 1:10
    temppeak= cnt(i,:); %get the specific point
    x1 = bsxfun(@minus,cnt,temppeak); %subtract the coordinates from the matrix of point coordinates
    x1 = x1.*x1;   %square the values
    dist(:,i) = sqrt (x1 (:,1) + x1 (:,2) + x1 (:,3)); %find distance
end

binranges = 0:binsize:maxbin; %the set of bins (can be non-linear)
[bincounts] = histc(dist,binranges); %this is the number of points in each bin
bincounts = sum (bincounts,2); %sum all distributions of points

% Now Normalize for size of area/radius
Nbincounts =zeros (1,numel(bincounts)); %create pre-loop
for i = 1:numel(bincounts)
    rin = binsize*(i-0.5); %inner radius of shell
    rout = binsize*(i+0.5); %outer radius of shell
   %Vshell = 4* pi * i*i *binsize*binsize * (rout-rin);
    %Vshell = 2*pi *binsize*i*(rout-rin);
    Vshell = 4*pi *binsize^2*i^2 *(rout-rin);
    %Vshell = (4/3)*pi*(rout^3 - rin^3); %volume of shell
    Nbincounts(i)= bincounts(i) * (dim(1)*dim(2)*dim(3)/sz/sz) * (1/Vshell); %copied off Dan Magagnosc
end

bar(binranges,Nbincounts,'histc')
title ('Radial Distribution Function')
xlabel ('pixels from center')
ylabel ('Normalized frequency')
%%

mxdisp = 40; %yes this is redundant. no don't delete it. 
field1 = 'mem';  value1 = (0); %how many frames a particle can be missing for
field2 = 'goodenough';  value2 = (8); %sometimes needs to be called 'good'.
field3 = 'dim';  value3 = (3);
field4 = 'quiet';  value4 = (0);
field5 = 'mxdisp';  value5 = (40);
field6 = 'method';  value6 = ('normal'); %normal is better

param = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);

cd ('Z:\USERS\Klebnikov\Matlab Codes\My Code') %go to code location.
trks2 = trackS(out2,mxdisp,param); %find tracks, in 3D

plot3 (trks2 (:,1), trks2 (:,2), trks2 (:,3),'x')
grid on
title ('tracking 20 particles over 9 time steps')

figure
rms3D %RMS and individual tracks script. 
