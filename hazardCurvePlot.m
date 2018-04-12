clear all; close all;clc;format short e;

%% PGA ref 2016
PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];

%% loading hazard curves from the 2016 model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('hazardCurve16.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
z5 = zeros(length(w1),1);
% dis, lat, lon, PGAREF(1) ... PGAREF(20)
forcast16 = [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,z5];
% removing western hazard curves.
forcast16(find(forcast16(:,2) < -115.000),:) = [];

%% loading the adaptive model. 2016
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('./branch-adp/Llenos_max_NSH_CEUS.pga.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
%     w3(find(w3 < 0.0001)) = 0.0;    w8(find(w8 < 0.0001)) = 0.0;    w13(find(w13 < 0.0001)) = 0.0;    w17(find(w17 < 0.0001)) = 0.0;
%     w4(find(w4 < 0.0001)) = 0.0;    w9(find(w9 < 0.0001)) = 0.0;    w13(find(w13 < 0.0001)) = 0.0;    w18(find(w18 < 0.0001)) = 0.0;
%     w5(find(w5 < 0.0001)) = 0.0;    w10(find(w10 < 0.0001)) = 0.0;    w14(find(w14 < 0.0001)) = 0.0;    w19(find(w19 < 0.0001)) = 0.0;
%     w6(find(w6 < 0.0001)) = 0.0;    w11(find(w11 < 0.0001)) = 0.0;    w15(find(w15 < 0.0001)) = 0.0;    w20(find(w20 < 0.0001)) = 0.0;
%     w7(find(w7 < 0.0001)) = 0.0;    w12(find(w12 < 0.0001)) = 0.0;    w16(find(w16 < 0.0001)) = 0.0;    w21(find(w21 < 0.0001)) = 0.0;
%     w22(find(w22 < 0.0001)) = 0.0;
%% no weight
Z = zeros(size(w1)); Z = Z + 0.5;
adaptiveNW = [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,...
           w16,w17,w18,w19,w20,w21,w22,Z];
    
clear w3;clear w4;clear w5;clear w6;clear w7;clear w8;
clear w9;clear w10;clear w11;clear w12; clear w13;clear w14;clear w15;
clear w16;clear w17;clear w18;clear w19;clear w20;clear w21; clear w22;

%% loading informed branches from logic tree. 
% reading branch names and their weights into a structure
[bn] = textread('weights_20162.txt','%s');
bns = char(bn);
aa = 0; sumW =0;
for vv = 1:length(bns)
    aa = aa+1;
CC = strsplit(char(bns(vv,:)),'_');
 if length(CC) == 9
sn = strcat(char(CC(1)),'_',char(CC(2)),'_',char(CC(3)),'_',char(CC(4)),'_',...
    char(CC(5)),'_',char(CC(6)),'_',char(CC(7)),'_',char(CC(8)));
sw = strcat(char(CC(9)));
 elseif length(CC) == 10
sn = strcat(char(CC(1)),'_',char(CC(2)),'_',char(CC(3)),'_',char(CC(4)),'_',...
    char(CC(5)),'_',char(CC(6)),'_',char(CC(7)),'_',char(CC(8)),'_',char(CC(9)));
sw = strcat(char(CC(10)));
 end
sn = sprintf('%s', sn);
WW{aa}.N = sn;
WW{aa}.W = sw;
sumW = sumW + str2num(sw) 
end
clear aa; clear bn; clear bns; clear sn;clear sw; clear vv; clear Z; clear z5
% cfator = 0.5/sumW % for normalizing weights

% [bw] = textread('weights_2016.txt','%f');
all = dir('branch-infv');branches = all(4:length(all));cc = 0;a =0;CK=0;
for bnum = 1:numel(branches);
    bnum
    C = strsplit(branches(bnum).name,'_'); CC = strsplit(char(C(8)),'.');  
    s = strcat(char(C(2)),'_',char(C(3)),'_',char(C(6)),'_',char(C(7)),'_',char(CC(1)),'_',char(C(4)));
    ss = strcat(char(C(1)),'_',char(C(2)),'_',char(C(3)),'_',char(C(4)),'_',char(C(5)),'_',char(C(6)),...
        '_',char(C(7)),'_',char(CC(1)));
    ss = sprintf('%s', ss);
    
    % finding the weight of the branch
    b = NaN;
    for zz = 1:length(WW)
       tf = strcmp(WW{zz}.N,ss);
       if tf == 1; b = str2num(WW{zz}.W);end
    end
    b; if isnan(b) == 1; 
        ss
    end
    
%     b = b.* cfator;
    CK = CK + b;
    % checking the cumulative weights
    Z = zeros(size(w1)); Z = Z + b;
    % reading the branch
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('./branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    % truncating
%     w3(find(w3 < 0.0001)) = 0.0;    w8(find(w8 < 0.0001)) = 0.0;    w13(find(w13 < 0.0001)) = 0.0;    w17(find(w17 < 0.0001)) = 0.0;
%     w4(find(w4 < 0.0001)) = 0.0;    w9(find(w9 < 0.0001)) = 0.0;    w13(find(w13 < 0.0001)) = 0.0;    w18(find(w18 < 0.0001)) = 0.0;
%     w5(find(w5 < 0.0001)) = 0.0;    w10(find(w10 < 0.0001)) = 0.0;    w14(find(w14 < 0.0001)) = 0.0;    w19(find(w19 < 0.0001)) = 0.0;
%     w6(find(w6 < 0.0001)) = 0.0;    w11(find(w11 < 0.0001)) = 0.0;    w15(find(w15 < 0.0001)) = 0.0;    w20(find(w20 < 0.0001)) = 0.0;
%     w7(find(w7 < 0.0001)) = 0.0;    w12(find(w12 < 0.0001)) = 0.0;    w16(find(w16 < 0.0001)) = 0.0;    w21(find(w21 < 0.0001)) = 0.0;
%     w22(find(w22 < 0.0001)) = 0.0;
    % putting together the branches with their associated weifghts
    brnw = [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,Z];
%     br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
%           b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
%     branchHazW{cc} = br; %% with weight
    branchHazNW{bnum}.C = brnw; %% without weight
    branchHazNW{bnum}.N = s; %% name
%     branchNmW{cc} = s; %% name
end
CK
branchHazNW{bnum+1}.C = adaptiveNW;
branchHazNW{bnum+1}.N = 'Llenos_max_NSH_CEUS';

DC = 0;
for k = 1 : length(branchHazNW)
    k
    tmp = branchHazNW{k}.C;  DC = DC + tmp(1,23)
end

% % putting together all info for 2016 model
% PGA2017.brnach = branchHazNW;
% PGA2017.USGSmean = forcast16;
% PGA2017.refPGA = PGAREF;

haz2016.brnach = branchHazNW;
haz2016.USGSmean = forcast16;
haz2016.refPGA = PGAREF;


save('haz2016','-v7.3');

load('haz2016');

% constructing the mean hazard
estmean = branchHazNW{1}.C;
estmean(:,3:23) = 0;
for i = 1:length(branchHazNW);
    temp = branchHazNW{i}.C;
    temp(:,3:22) = temp(:,3:22).*temp(1,23); 
    estmean(:,3:22) = estmean(:,3:22) + temp(:,3:22); 
end


%% loading informed branches from logic tree. 
% reading branch weights.
[bw] = textread('weights_2017.txt','%f');

all = dir('17branch-infv');
branches = all(4:length(all));
cc = 0;



%% Percentile
tpl = branchHazNW{1};
perct15 = branchHazNW{1};
perct85 = branchHazNW{1};
perct = branchHazNW{1};

for i = 1:length(branchHazW{1});
for ii = 3 : 22
    vartmp = [];
    for k = 1:length(branchHazW);
    temp = branchHazW{k};  
%     vartmp = [vartmp, temp(i,ii)];
% %     temp = branchHazW{k};
    if length(nonzeros(temp(i,:))) > 3; vartmp = [vartmp, temp(i,ii)]; end
    end  
    perct15(i,ii) = prctile(nonzeros(vartmp),15);
    perct85(i,ii) = prctile(nonzeros(vartmp),85);
    perct(i,ii) = prctile(nonzeros(vartmp),85)-prctile(nonzeros(vartmp),15);
end
end
    
%% plotting hazard curves for a single point
% lat = 39 ; lon = -84; % Kentucky-Ohaio Boreder
% lat = 35.5 ; lon = -83; % Eastern TN
% lat = 36.60 ; lon = -89.55; % New-Madrid

lat = 32.8 ; lon = -96.8; % Dallas
% lat = 35.45; lon = -97.5; % Oklahoma City
station = find(PGA2017.USGSmean(:,1) == lat & PGA2017.USGSmean(:,2) == lon) 

figure
% plotting individual branches
for lp = 1:length(PGA2017.brnach)-1
    
an = PGA2017.brnach{lp}.C;
anR1 = an(station,3:22);
anR1 = anR1';
%%%%%%%%%%%%%%%
% anR1 = anR1.*haz2016.brnach{lp}.C(1,23); %% weighting
%%%%%%%%%%%%%%%%
anR11 = anR1(:) > 10e-11;
s = nonzeros(anR11);
if length(s) > 3
ed = length(s) + 2;
hold on 
C = strsplit(PGA2017.brnach{lp}.N,'_'); 
tf1 = strcmp(char(C(4)),'Mx6p0'); tf2 = strcmp(char(C(5)),'CEUS'); tf3 = strcmp(char(C(2)),'PI');
ft = tf1*tf2;

if tf3 == 1
h1 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),'--r','LineWidth',2)
elseif tf3 == 0  
h2 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),'--b','LineWidth',2)
end
hold on 
if ft == 1   
h3 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),':k^','LineWidth',2,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63])
elseif ft == 0
h4 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),':yo','MarkerSize',7.5,'MarkerEdgeColor','k','MarkerFaceColor','r')
end
end
hold on 
end

% plotting adaptive curve 2016
anR1 = PGA2017.brnach{length(PGA2017.brnach)}.C(station,3:22);
anR1 = anR1';
anR11 = anR1(:) > 10e-11;
s = nonzeros(anR11);
ed = length(s) + 2;
h5 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),'-ms','LineWidth',2,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')

% plotting USGS mean
anR1 = PGA2017.USGSmean(station,3:22);
anR1 = anR1';
anR11 = anR1(:) > 10e-11;
s = nonzeros(anR11);
ed = length(s) + 2;
h6 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),'-g*','LineWidth',2,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')

% % plotting informed curve 2016
% anR1 = informedNW(station,:);
% anR1 = anR1';
% anR11 = anR1(3:22) > 10e-11;
% si = nonzeros(anR11);
% edi = length(si) + 2;
% h6 = plot((PGAREF(1:length(si))),(anR1(3:edi)),'-ks','LineWidth',2,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')

% plotting estimated mean
anR1 = estmean(station,3:22);
anR1 = anR1';
anR11 = anR1(:) > 10e-5;
s = nonzeros(anR11);
ed = length(s) + 2;
h7 = plot((PGA2017.refPGA(1:length(s))),(anR1(1:length(s))),'--ks','LineWidth',2,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on 



% perct
anR1 = perct15(i,:);
anR1(isnan(anR1)) = 0;
anR1 = anR1';
anR11 = anR1(3:22);
s = nonzeros(anR11);
ed = length(s) + 2;
h8 = plot((PGAREF(1:length(s))),(anR1(3:ed)),'--r','LineWidth',2.7,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')

% perct
anR1 = perct85(i,:);
anR1(isnan(anR1)) = 0;
anR1 = anR1';
anR11 = anR1(3:22);
s = nonzeros(anR11);
ed = length(s) + 2;
h9 = plot((PGAREF(1:length(s))),(anR1(3:ed)),'--g','LineWidth',2.7,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')





% ,'MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]

% hold on
% % plotting individual branches
% for i = 1: length(branchHazW17)
% an = branchHazW17{i};
% anR1 = an(station,:);
% anR1 = anR1';
% anR11 = anR1(3:22) > 10e-11;
% s = nonzeros(anR11);
% ed = length(s) + 2;
% hold on 
% 
% C = strsplit(branchNmW{i},'_'); 
% tf1 = strcmp(char(C(4)),'Mx6p0'); tf2 = strcmp(char(C(5)),'CEUS'); tf3 = strcmp(char(C(2)),'PI')
% ft = tf1*tf2;
% 
% if ft == 1 & tf3 ==1
% h8 = plot((PGAREF(1:length(s))),(anR1(3:ed)),'-b^','LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor','b')
% elseif ft == 0 & tf3 ==1  
% h9 = plot((PGAREF(1:length(s))),(anR1(3:ed)),':b^','LineWidth',1,'MarkerSize',7.5)
% elseif ft == 1 & tf3 ==0  
% f1 = plot((PGAREF(1:length(s))),(anR1(3:ed)),'-r^','LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor','r')
% elseif ft == 0 & tf3 ==0  
% f2 = plot((PGAREF(1:length(s))),(anR1(3:ed)),':r^','LineWidth',1,'MarkerSize',7.5)
% end
% 
% hold on 
% end
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% hold on 
% 
% % plotting adaptive curve 2017
% anR1 = adaptiveNW17(station,:);
% anR1 = anR1';
% anR11 = anR1(3:22) > 10e-11;
% s = nonzeros(anR11);
% ed = length(s) + 2;
% f3 = plot((PGAREF(1:length(s))),(anR1(3:ed)),'-md','LineWidth',1.7,'MarkerSize',7.5)
% 
% % plotting informed curve 2017
% anR1 = informedNW17(station,:);
% anR1 = anR1';
% anR11 = anR1(3:22) > 10e-11;
% si = nonzeros(anR11);
% edi = length(si) + 2;
% f4 = plot((PGAREF(1:length(si))),(anR1(3:edi)),'-kd','LineWidth',1.7,'MarkerSize',7.5)
% 
% % plotting 2016 model 2017
% anR1 = forcast17(station,:);
% anR1 = anR1';
% anR11 = anR1(4:23) > 10e-11;
% s = nonzeros(anR11);
% ed = length(s) + 3;
% f5 = plot((PGAREF(1:length(s))),(anR1(4:ed)),'-gd','LineWidth',1.7,'MarkerSize',7.5)

title('East TN No-Weight')
grid on;
set(gca,'LineWidth',2)
ylabel('Annual Frequency of Exceedances');
xlabel('Peak Ground Acceleration (g)');
% lgd = legend([h1 h2 h3 h4 h3 h4 h5 h6 h7 h8 h9 f2 f3 f4 f5],{'2016-PI_Mx6_CEUS','2016_PI','2016_Mx6_CEUS', ...
%     '2016_rest of branches','2016_Adaptive','2016_Informed','2016_Final','2017_PI_Mx6_CEUS','2017_PI', ...
%     '2017_rest of branches','2017_Adaptive','2017_Informed','2017_Final'}, ...
%     'Location','northeast','Orientation','vertical');
% 
lgd = legend([h1 h2 h3 h4 h5 h6 h7],{'2016.PI','2016.PN','2016.Mx6.CEUS','2016.Rest of Branches','2016.Adaptive',...
    '2016-Informed','2016.Final(mean hazard)'}, ...
    'Location','northeast','Orientation','vertical');

lgd = legend([h8 h9 h7],{'15-Percentile','85-Percentile','Mean hazard)'}, ...
    'Location','northeast','Orientation','vertical');

% lgd.FontSize = 18;
% lgd.TextColor = 'k';

xlim([0.005 2.5])

set(gca,'fontsize',21)

figure

% plotting 2016 model 2016
anR1 = forcast16(station,:);
anR1 = anR1';
anR11 = anR1(4:23) > 10e-11;
s = nonzeros(anR11);
ed = length(s) + 3;
h7 = plot((PGAREF(1:length(s))),(anR1(4:ed)),'-k*','LineWidth',2.7,'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','c')



set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on 
%% reformating standard dv files of Q Canada
% [w1,w2,w3,w4,w5] = textread('U.txt','%f %f %f %f %f');
% delete('onePercent_1YrModel.pga.1pc1.txt')
% fileID = fopen('onePercent_1YrModel.pga.1pc1.txt','w');
% for k = 1 : length(w1)
%     fprintf(fileID,'%3.3f %2.3f %f\n',w1(k,:),w2(k,:),w5(k,:)/100);
% end
% fclose(fileID);