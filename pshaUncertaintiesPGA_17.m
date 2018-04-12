clear all; close all;clc; format short e;

PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];

load('PGA2017');

TW = 0;% finding the total weight
for k = 1 : length(PGA2017.brnach)
    k
    tmp = PGA2017.brnach{k}.C;  TW = TW + tmp(1,23)
end

%% Bootstraping 
Nsampling = 1000;
for k = 1: Nsampling; 
    pk = 100.*(k/Nsampling);
    X = [num2str(pk),' % of the Bootstraping is done.'];
    disp(X)
WW = 0;
r = linspace(1,length(PGA2017.brnach),length(PGA2017.brnach));
ns = floor(length(PGA2017.brnach)*0.9);
y = datasample(r,ns); 
informedNW = PGA2017.brnach{y(1)}.C;
bm = informedNW; WW = informedNW(1,23);
informedNW(:,3:22).*WW;
for i = 2:length(y);
    temp = PGA2017.brnach{y(i)}.C;
    WT = temp(1,23); WW = WW + WT;
    informedNW(:,3:22) = informedNW(:,3:22) + temp(:,3:22).*WT;
end  
    RW = TW - WW; y2 = datasample(y,1); 
    temp2 = PGA2017.brnach{y2}.C;
    bm(:,3:22) = informedNW(:,3:22)+temp2(:,3:22).*RW;
    bootmodel{k} = bm;
end

informedNW(:,3:23)=0; mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    EP = (E./length(informedNW))*100;
    X = [num2str(EP),' % of the Standard Deviation calculation is done.'];
    disp(X)
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
temp2 = PGA2017.brnach{y(1)}.C;
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
          yi = PGAREF(index1); ui = U(index1);
   elseif closestValues1 > 0.0101
      index2 = index1 + 1; closestValues2 = H(index2); 
   pga1 = PGAREF(index1); % finding the associated PGAs 
   pga2 = PGAREF(index2);
   y = [pga1;pga2]; x = [closestValues1;closestValues2];
   yi = interp1(x,y,0.0101,'linear'); % 'nearest' 'linear
   u = [U(index1);U(index2)];
   ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear

   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1); ui = U(index1);
       else
          index2 = index1; index1 = index2 - 1;
          closestValues1 = H(index1); closestValues2 = H(index2);
          pga1 = PGAREF(index1); % finding the associated PGAs 
          pga2 = PGAREF(index2);
          y = [pga1;pga2]; x = [closestValues1;closestValues2];
          yi = interp1(x,y,0.0101,'linear'); % 'nearest' 'linear
          u = [U(index1);U(index2)];
          ui = interp1(x,u,0.0101,'linear'); % 'nearest' 'linear
       end
   end
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); haz(j,3) = yi; 
   ucr(j,1) = temp2(j,2); ucr(j,2) = temp2(j,1); ucr(j,3) = ui;   
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

delete('one17Percent_1YrModel.pga.1pc1.txt')
fileID = fopen('one17Percent_1YrModel.pga.1pc1.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f\n',haz(k,:));
end
fclose(fileID);

delete('one17Percent_1YrModel.pga.1pc2.txt')
fileID = fopen('one17Percent_1YrModel.pga.1pc2.txt','w');
for k = 1 : length(ucr)
fprintf(fileID,'%3.3f %2.3f %f\n',ucr(k,:));
end
fclose(fileID);

% mean std
sum(ucr(:,3))./length(ucr)

delete('2017_PGA_Boot_StD.txt')
fileID = fopen('2017_PGA_Boot_StD.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',stboot(k,1:22));
end
fclose(fileID);

