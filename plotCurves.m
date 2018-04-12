clear all; close all;clc;format short e;

% I have saved all the truncated unweighted branches into a single matlab
% structure.
% haz2016.brnach{#branch}.C is the branch values with first two column being 
% lat and long, following by the annual exceedance rates in next 20 column 
% and the last column is the weight of the curve. 
% haz2016.brnach{#branch}.N gives the name of the branch
% haz2016.USGSmean is the mean hazard published by USGS
% haz2016.refPGA are the reference PGAs associated to column 3 to 22;

load('2016Model');

% constructing the mean hazard
estmean = haz2016.brnach{1}.C;
estmean(:,3:23) = 0;
for i = 1:length(haz2016.brnach);
    temp = haz2016.brnach{i}.C;
    temp(:,3:22) = temp(:,3:22).*temp(1,23); 
    estmean(:,3:22) = estmean(:,3:22) + temp(:,3:22); 
end

%% plotting hazard curves for a single point
% lat = 39 ; lon = -84; ttl= 'Kentucky-Ohaio Boreder';
% lat = 35.5 ; lon = -83; ttl= 'Eastern TN';
% lat = 36.60 ; lon = -89.55; ttl= 'New-Madrid';
lat = 32.8 ; lon = -96.8; ttl= 'Dallas';
% lat = 35.45; lon = -97.5; ttl= 'Oklahoma City';

station = find(haz2016.USGSmean(:,1) == lat & haz2016.USGSmean(:,2) == lon) 

figure
%% plotting individual branches
for lp = 1:length(haz2016.brnach)-1
    
an = haz2016.brnach{lp}.C;
anR1 = an(station,3:22);
anR1 = anR1';
%%%%%%%%%%%%%%
anR1 = anR1.*(haz2016.brnach{lp}.C(1,23)); %% weighting
%%%%%%%%%%%%%%%
anR11 = anR1(:) > 10e-11;
s = nonzeros(anR11);
if length(s) > 3
ed = length(s) + 2;
hold on 
C = strsplit(haz2016.brnach{lp}.N,'_'); 
tf1 = strcmp(char(C(4)),'Mx6p0'); tf2 = strcmp(char(C(5)),'CEUS'); tf3 = strcmp(char(C(2)),'PI');
ft = tf1*tf2;

if tf3 == 1
h1 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'--k','LineWidth',1)
elseif tf3 == 0  
h2 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'-+b','LineWidth',1)
end
hold on 
if ft == 1   
h3 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),':k*','LineWidth',1,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63])
elseif ft == 0
h4 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),':ro','LineWidth',1,'MarkerSize',6.5,'MarkerEdgeColor','k','MarkerFaceColor','k')
end
end
hold on 
end

%plotting USGS mean
anR1 = haz2016.USGSmean(station,3:22);
anR1 = anR1';
anR11 = anR1(:) > 10e-11;
s = nonzeros(anR11);
ed = length(s) + 2;
h6 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'-g^','LineWidth',2.5,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','c')

%plotting estimated mean
anR1 = estmean(station,3:22);
anR1 = anR1';
anR11 = anR1(:) > 10e-5;
s = nonzeros(anR11);
ed = length(s) + 2;
h7 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'--ro','LineWidth',2.5,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on 

title(ttl)
grid on;
set(gca,'LineWidth',2)
ylabel('Annual Frequency of Exceedances');
xlabel('Peak Ground Acceleration (g)');

lgd = legend([h1 h2 h3 h4 h6 h7],{'2016.Potentialy Indiced','2016.Potentialy Natureal','2016.Mx6.CEUS','2016.Rest of Branches',...
    '2016-IUSGS mean','2016.Estimated Mean)'}, ...
    'Location','northeast','Orientation','vertical');

xlim([0.005 2.5])
set(gca,'fontsize',21)
