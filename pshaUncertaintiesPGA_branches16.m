clear all; close all;clc;

PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];

load('2016Model');

%% Bootstraping 
Nsampling = 500;

%%%%%%%%%%%%%%%%%%%  Mx6p0 GAIL2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');
    ft = tf1*tf2;
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling; 
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
Mx6p0_GAIL2 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_GAIL2);

%%%%%%%%%%%%%%%%%%%  Mx6p0 GAIL5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');
    ft = tf1*tf2;
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling; 
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
Mx6p0_GAIL5 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_GAIL5);

%%%%%%%%%%%%%%%%%%%  Mx6p0 CEUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');
    ft = tf1*tf2;
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling;
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
   index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
Mx6p0_CEUS = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',Mx6p0_CEUS);

%%%%%%%%%%%%%%%%%%%  Mx6p0 R10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    ft = strcmp(R,'R10');
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling;
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
R_10 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',R_10);

%%%%%%%%%%%%%%%%%%%  Mx6p0 R20 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    ft = strcmp(R,'R20');
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling; 
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
R_20 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',R_20);

%%%%%%%%%%%%%%%%%%%  Mx6p0 D_1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    ft = strcmp(char(C(2)),'01')  
    if ft == 1
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling; 
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.
   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
D_1 = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',D_1);


%%%%%%%%%%%%%%%%%%%  Mx6p0 D_LT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with weight
cc=0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    ft = strcmp(char(C(2)),'01')  
    if ft == 0
    cc = cc+1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22).*(tmp(1,23));
    branchHazW{cc} = tmp;
    end
end

%% Bootstraping 
for k = 1: Nsampling; 
WW = 0;
r = linspace(1,length(branchHazW),length(branchHazW));
ns = floor(length(branchHazW)*0.9);
y = datasample(r,ns); 
informedNW = branchHazW{y(1)};
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = branchHazW{y(i)};
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = 4.982- WW; y2 = datasample(y,1); 
    temp2 = branchHazW{y2};
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    for EE = 3 : 22
        vartmp = [];
        for EEE = 1:Nsampling
            tempp = bootmodel{EEE}; vartmp = [vartmp, tempp(E,EE)];
        end   
    mean(E,EE) = sum(vartmp)/length(vartmp);
    stboot(E,EE) = 1.96*(std(vartmp));
    end   
end

%% extracting the hazard form curves
haz = zeros(length(stboot),3);
ucr = zeros(length(stboot),3);
ucrP = zeros(length(stboot),3);
temp2 = haz2016.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
       ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
    x = [closestValues1;closestValues2];
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
           ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          x = [closestValues1;closestValues2];
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); 
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

% mean std
D_LT = sum(ucr(:,3))./length(ucr)
fprintf(fileID,'%f\n',D_LT);
fclose(fileID);