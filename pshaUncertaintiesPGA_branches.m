clear all; close all;clc;

PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];

%% loading the adaptive model. 
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('../../data/beroza/mmousavi/branch-adp/Llenos_max_NSH_CEUS.pga.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

% %% no weight
% adaptiveNW = [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,...
%            w16,w17,w18,w19,w20,w21,w22];
%        
%% with weight 
adaptive = [w1,w2,(0.5).*w3,(0.5).*w4,(0.5).*w5,(0.5).*w6,(0.5).*w7,(0.5).*w8,...
        (0.5).*w9,(0.5).*w10,(0.5).*w11,(0.5).*w12,(0.5).*w13,(0.5).*w14,(0.5).*w15,...
        (0.5).*w16,(0.5).*w17,(0.5).*w18,(0.5).*w19,(0.5).*w20,(0.5).*w21,(0.5).*w22];   

%% loading informed branches from logic tree. 
% reading branch weights.
[bw] = textread('weights_2016.txt','%f');


%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');
    ft = tf1*tf2;
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std.
fileID = fopen('2016_branch.txt','w');
MxNSH_CEUS = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',MxNSH_CEUS);

%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');
    ft = tf1*tf2;
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
Mx6p0_GAIL2 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_GAIL2);

%%%%%%%%%%%%%%%
%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');
    ft = tf1*tf2;
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
Mx6p0_GAIL5 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_GAIL5);

%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');
    ft = tf1*tf2;
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
Mx6p0_CEUS = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_CEUS);
%%%%%%%%%%%%%%%%%%%

%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
 
    ft = strcmp(R,'R10')
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
R_10 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',R_10);
%%%%%%%%%%%%%%%%%%%%%

%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
 
    ft = strcmp(R,'R20')
%      ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
R_20 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',R_20);
%%%%%%%%%%%%%%%%
%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
 
    ft = strcmp(char(C(2)),'01')
    if ft == 1
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
D_1 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',D_1);
%%%%%%%%%%%%%%

%% with weight
cc = 0;
all = dir('../../data/beroza/mmousavi/branch-infv');
branches = all(4:length(all));
for bnum = 1:numel(branches); 
    b = bw(bnum);
    C = strsplit(branches(bnum).name,'_'); 
    CC = strsplit(char(C(8)),'.');
    R = char(C(6)); Max = char(C(7)); GMPE = char(CC(1)); 
 
    ft = strcmp(char(C(2)),'01')
    if ft == 0
    cc = cc+1
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread(sprintf('../../data/beroza/mmousavi/branch-infv/%s',branches(bnum).name),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    br = [w1,w2,b.*w3,b.*w4,b.*w5,b.*w6,b.*w7,b.*w8,b.*w9,b.*w10,b.*w11,b.*w12,b.*w13,...
    b.*w14,b.*w15,b.*w16,b.*w17,b.*w18,b.*w19,b.*w20,b.*w21,b.*w22];
    branchHazW{cc} = br;
       end
end

%% Bootstrapping
% constructing informed model
for k = 1:500;
r = randi([2 length(branchHazW)],1,floor(length(branchHazW)*0.9));
informedNW = branchHazW{1};
for i = 1:length(r);
    iii = r(i);
for ii = 3 : 22
    temp = branchHazW{iii};
    informedNW(:,ii) = informedNW(:,ii) + temp(:,ii);
    modelR(:,ii) = informedNW(:,ii) + adaptive(:,ii);
end   
end
   bootmodel{k} = modelR;
end
% load('informedNW')

mean = branchHazW{1}; stboot = branchHazW{1};
for i = 1:length(informedNW);
    for ii = 3 : 22
        vartmp = [];
        for iii = 1:500
            tempp = bootmodel{iii};
            vartmp = [vartmp, tempp(i,ii)];
        end   
    mean(i,ii) = sum(vartmp)/length(vartmp);
    stboot(i,ii) = std(vartmp);
    end   
end


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
          zi = (ui)./0.0101;
          
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
   zi = (ui)./0.0101;   
   
   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
          zi = (ui)./0.0101;
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
          zi = (ui)./0.0101;  
       end
   end
 
   haz(j,1) = temp2(j,2);
   haz(j,2) = temp2(j,1);
   haz(j,3) = yi; 
   
   ucr(j,1) = temp2(j,2);
   ucr(j,2) = temp2(j,1);
   ucr(j,3) = zi;
   
  
end
ucrP(:,3) = ucrP(:,3)./max(ucrP(:,3));
haz(isnan(haz)) = 0;
ucr(isnan(ucr)) = 0;
ucrP(isnan(ucr)) = 0;


% mean std
D_LT = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',D_LT);
fclose(fileID);