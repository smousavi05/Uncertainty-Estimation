clear all; close all;clc;format short e;

PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];

%% loading the adaptive model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
textread('../../../data/beroza/mmousavi/17branch-adp/Llenos_max_NSH_CEUS_RS.5hz.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
        
%% with weight 
adaptive = [w1,w2,(0.5).*w3,(0.5).*w4,(0.5).*w5,(0.5).*w6,(0.5).*w7,(0.5).*w8,...
        (0.5).*w9,(0.5).*w10,(0.5).*w11,(0.5).*w12,(0.5).*w13,(0.5).*w14,(0.5).*w15,...
        (0.5).*w16,(0.5).*w17,(0.5).*w18,(0.5).*w19,(0.5).*w20,(0.5).*w21,(0.5).*w22];   

%% loading informed branches from logic tree. 
% reading branch weights.
[bw] = textread('weights_2017.txt','%f');

all = dir('../../../data/beroza/mmousavi/17branch-infv5Hz');
branches = all(4:length(all));
cc = 0;

% with weight
for bnum = 1:numel(branches);
    b = bw(bnum);
       cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../../data/beroza/mmousavi/17branch-infv5Hz/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
        b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    aa = br(:,3:22);
    aa(find(aa <=0.0001)) = 0.0; % truncatinng the curves. 
    br(:,3:22) = aa;
    branchHazW{cc} = br;
end

%% Bootstrapping
% constructing informed model
Nsampling = 1000;
for k = 1: Nsampling; k
r = randi([2 171],1,floor(171*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
    temp = branchHazW{iii};
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22);
    modelR(:,3:22) = informedNW(:,3:22) + adaptive(:,3:22);  
end
   bootmodel{k} = modelR;
end

mean = branchHazW{1}; stboot = branchHazW{1};ciboot = branchHazW{1};
for i = 1:length(informedNW);
    i
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:1000
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end

%% loading hazard curves from the 2016 model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('hazardCurve17_5Hz.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
z5 = zeros(length(w1),1);
% dis, lat, lon, PGAREF(1) ... PGAREF(20)
forcast16 = [z5,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,z5,z5];
% removing western hazard curves.
forcast16(find(forcast16(:,3) < -115.000),:) = [];


%% extracting the hazard form curves
haz = zeros(length(branchHazW{1}),3);
ucr = zeros(length(branchHazW{1}),3);
ucrP = zeros(length(branchHazW{1}),3);
temp2 = branchHazW{1};
for j = 1 : length(mean)
   H = mean(j,3:22);
   U = stboot(j,3:22);
   
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
          yi = PGAREF(index1);
          ui = U(index1);
          
   elseif closestValues1 > 0.0101
      index2 = index1 + 1;
      closestValues2 = H(index2); 
      
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
   y = [pga1;pga2];
   x = [closestValues1;closestValues2];
   yi = interp1(x,y,0.0101,'linear'); % 'nearest' 'linear
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear 
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
       else
          index2 = index1; 
          index1 = index2 - 1;
          closestValues1 = H(index1);
          closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          y = [pga1;pga2];
          x = [closestValues1;closestValues2];
          yi = interp1(x,y,0.0101,'linear'); % 'nearest' 'linear
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear 
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = ui;
   
   ucrP(j,1) = temp2(j,2);
   ucrP(j,2) = temp2(j,1);
   ucrP(j,3) = ui./yi;
   
end
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;

delete('one17Percent_1YrModel.5Hz.1pc1.txt')
fileID = fopen('one17Percent_1YrModel.5Hz.1pc1.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f\n',haz(k,:));
end
fclose(fileID);

delete('one17Percent_1YrModel.5Hz.1pc2.txt')
fileID = fopen('one17Percent_1YrModel.5Hz.1pc2.txt','w');
for k = 1 : length(ucr)
fprintf(fileID,'%3.3f %2.3f %f\n',ucr(k,:));
end
fclose(fileID);

delete('one17Percent_1YrModel.5Hz.1pc3.txt')
fileID = fopen('one17Percent_1YrModel.5Hz.1pc3.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f\n',ucrP(k,:));
end
fclose(fileID);


delete('2017_5Hz_Boot_StD.txt')
fileID = fopen('2017_5Hz_Boot_StD.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',stboot(k,:));
end
fclose(fileID);
