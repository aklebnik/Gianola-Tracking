%testing out Eric Weeks Particle tracking script
%Sasha Klebnikov, 4/6/14
%last edited 17/6/14 for tests taken from our microscope

% 1) set stuff up
% 2) Find Basic Centers
% 3) Find subpixel center
% 4) Radial Distribution Function
% 5) Packing Fraction
% 6) Final plot
% 6) Other temporary code
%profile on

%% Set stuff up
clear; close all;
Brightness = 8; %set arbitrary minimum brightness, keep it pretty low
a = double(imread('stack3 1_1_page_0071.tif')); %512x512 pixels, some dark spots.
colormap ('gray'), imagesc(a)
Voronoi =0; %Turn Voronoi on/off - slow!
Fast = 0; %if 1, goes much faster.

if Fast ==0
    %ask to display intermediate plots
    plotsQ=questdlg('Do you want to display intermediate graphs?');
    if strcmp(plotsQ, 'Yes')==1
        plots = 1;
        disp ('showing intermediate plots')
    elseif strcmp (plotsQ, 'Cancel')==1 %terminate code
        disp ('you hate Sashas code...')
        return
    else
        plots =0;
    end
else
    plots = 0;
end


%% Start Code

%Crop the image
% crop=questdlg('Do you want to crop the image?');
% if strcmp(crop, 'Yes')==1
%     [a,RECT]=imcrop(a); %crop image
% end

%Set Pwidth
colormap ('gray'), imagesc(a(1:50, 1:50))
prompt={'Pixels (even is better)'};
answer=inputdlg(prompt, 'Pwidth of the largest Particles',1,{'28'});
Pwidth = str2double (answer{1});
%Pwidth = 50; %set arbitrarily
dim = size (a)-2*Pwidth;


% a1 = a(:,:,1); %if a appears 3D, do this.
% a2 = a(:,:,2);
% a3 = a(:,:,3);
% a = (a1+a2+a3)/3;

%a = 255-a;  %Invert image if necessary.
b = bandpass(a,1,dim,plots,Pwidth); %create a new variable, b, which has been filtered
%first term is data set (a),
%second (min wavelength) should default to 1, but for Pwidth >40, try 2,3 or higher.
%third is the average particle diameter

%% Find basic centers
%b=b*255;
pk = pkfnd(b,Brightness,Pwidth); %second value is minimum particle intensity
%third value is slightly larger than third value of bpass


if plots ==1
    close all
    figure(1)
    pk2(:,2) = length(b)-pk(:,2); %invert y values, as image center is top left, not bottom left
    plot (pk(:,1),pk2(:,2),'o') %plot just the particle size
    
    figure(2)
    colormap('gray'), imagesc(b);
    hold on
    plot (pk(:,1),pk(:,2),'x')
end

%example of filter
%b = b-25;
%b(b<0) = 0; %this is the most simplistic offset system you can use.
%note, pkfnd gives same answer afterwords. which makes sense, as the
%maximum value is still there

%% Find better centers
cnt = cntrd(b,pk,Pwidth+1); %if you loose a lot of data (ie size (pk) > size (cnt), reduce third term
%if getting annoying integer error, make sure Pwidth+num == odd. alter num
%to suit.

if plots ==1
    figure (1)
    plot (cnt(:,1),cnt(:,2),'o')
    
    figure(3)
    colormap('gray'), imagesc(b);
    hold on
    plot (cnt(:,1),cnt(:,2),'x')
end

%% Calculate if points are too close together
if Fast == 0 %if doing intermediate steps
    sz = length (cnt);
    dist = zeros(sz); %create preloop
    cnt2 = cnt;
    cnt2 (:,3:end) = []; %get rid of superflous data, to run faster
    
    
    for i = 1:sz
        %for i = 1:10
        temppeak= cnt2(i,:); %get the specific point
        x1 = bsxfun(@minus,cnt2,temppeak); %subtract the coordinates from the matrix of point coordinates
        x1 = x1.*x1;   %square the values
        dist(:,i) = sqrt (x1 (:,1) + x1 (:,2)); %find distance
    end
    
    dist (dist > (Pwidth/2)) = 0;
    numel (nonzeros(dist)); %tell me how many points are within x pixels.
    
    clear cnt2; clear sz; clear x1; clear temppeak; clear i
end
%% Voronoi Construction - takes FOREVER
if Voronoi ==1
    figure (2)
    colormap('cool'), imagesc(b);
    hold on
    plot (cnt(:,1),cnt(:,2),'x')
    
    center_positions = cnt (:,1:2);
    if plots ==1
        [V,E,C,bdry_cells]=vertex_edges_cells1(center_positions,'ON');
    else
        [V,E,C,bdry_cells]=vertex_edges_cells1(center_positions,'OFF');
    end
end

%% Radial Distribution Function
if Fast ==0
    [binranges, Nbincounts]=Radial_Dist_Fun(cnt, dim, plots,1,60);
    % cnt - contains the x,y coordinates of points in the first two rows
    % dim - dimensions of window in pixels
    % plots - if 1, will plot lots of intermediate graphs
    % binsize - how large each bin will be (pixels)
    % maxbin - what distance do we explore up to?
    
end

%% Packing Fraction (very rough)
if Fast ==0
    Parea =length (pk) *pi*(Pwidth/2)^2;
    PackingFraction = Parea/(dim(1)*dim(2)); %note, this uses pk, not cnt, to get more accurate readings.
end

%what about using max (Nbincounts) to find the maximum point of g(r), and
%then subtracting 10% to use as packing fraction.

%this is because the optimum radius for filtering may not be the actual
%width.

%% make some pretty graphs, to check accuracy
if plots ==1
    close(3)
end


figure(1)
subplot(2,3,1)
colormap('gray'), imagesc(b);
hold on
plot (cnt(:,1),cnt(:,2),'bx','Markersize',2)
title('filtered image')
axis equal

subplot(2,3,2)
cnt2(:,2) = length(b)-cnt(:,2); %invert y values, as image center is top left, not bottom left
plot (cnt(:,1),cnt2(:,2),'o') %plot just the particle size
title('centers of particles only')
axis equal

subplot(2,3,3)
plot (sqrt(cnt(:,4)),cnt(:,3), 'x')
xlabel ('radius of gyration')
ylabel('brightness')

subplot(2,3,4)
hist(mod(cnt(:,1),1),30);
title('modulus of location (should be flat)')

if Fast ==0
    subplot(2,3,5)
    bar(binranges,Nbincounts,'histc')
    title ('Radial Distribution Function')
    xlabel ('pixels from center')
    ylabel ('Normalized frequency')
end

subplot(2,3,6)
b1 = b (Pwidth:40+Pwidth, Pwidth:40+Pwidth);
colormap ('gray'), imagesc(b1);
title ('close up of quality')
axis equal

if Fast == 0
    figure (2)
    axis equal
    colormap('gray'), imagesc(a);
    hold on
    plot (cnt(:,1),cnt(:,2),'rx','Markersize',10)
    title('filtered image with points and Voronoi')
    if Voronoi ==1
        for ic=1:size(E,1)
            line([V(E(ic,1),1) V(E(ic,2),1)],[V(E(ic,1),2) V(E(ic,2),2)],'color',[1 0 0]);
        end;
        drawnow;
    end
end

%% finish up
clear prompt, clear b1, clear pk2, clear cnt2, clear Nbincounts, clear binranges,
if Fast ==0
    Pwidth
    Brightness
   % PackingFraction
    LostParticles = length (pk)-length (cnt)
    FoundParticles = length (cnt)
    ParticlesTooClose = numel (nonzeros(dist)) %how many centers are inside one particle
end

%profile viewer
%% Other code ideas
% close all
% colormap('gray'), imagesc(c);

% for i = 1:30
%     n = i*100 + 2500
% c = b>n;
%   c(c<1) = [];
%   d(i) =length(c)
% end
% plot (d)

%% to find minimum particle intensity
% subplot (2,3,4)
%  plot (max(b)) % and see the highest value in each column
%
% cs = b(:); %create a temporary testing matrix, concatonated into a column vector
% cs (cs<(Brightness/2)) = []; %get rid of lower values, to better see curve
% subplot (2,3,5)
% plot(sort(cs)) %plot results
% subplot (2,3,6)
% semilogx(sort(cs)) %another option

%Not using this, as its much less effective than simply observing by eye



