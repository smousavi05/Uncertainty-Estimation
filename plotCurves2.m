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

%%%%%%%%%%%%%%%%%%%  CEUS_Mx6p0_R10_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_Mx6p0_R10_01 = haz2016.brnach{1}.C;
CEUS_Mx6p0_R10_01(:,3:23) = 0;
CH = 0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1; 
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_Mx6p0_R10_01(:,3:22) = CEUS_Mx6p0_R10_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{1} = CEUS_Mx6p0_R10_01;

%%%%%%%%%%%%%%%%%%%  CEUS_Mx6p0_R10_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_Mx6p0_R10_02 = haz2016.brnach{1}.C;
CEUS_Mx6p0_R10_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_Mx6p0_R10_02(:,3:22) = CEUS_Mx6p0_R10_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{2} = CEUS_Mx6p0_R10_02;

%%%%%%%%%%%%%%%%%%%  CEUS_Mx6p0_R20_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_Mx6p0_R20_01 = haz2016.brnach{1}.C;
CEUS_Mx6p0_R20_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_Mx6p0_R20_01(:,3:22) = CEUS_Mx6p0_R20_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{3} = CEUS_Mx6p0_R20_01;
%%%%%%%%%%%%%%%%%%%  CEUS_Mx6p0_R20_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_Mx6p0_R20_02 = haz2016.brnach{1}.C;
CEUS_Mx6p0_R20_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_Mx6p0_R20_02(:,3:22) = CEUS_Mx6p0_R20_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{4} = CEUS_Mx6p0_R20_02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% GAIL2_Mx6p0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Mx6p0_R10_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL2_Mx6p0_R10_01 = haz2016.brnach{1}.C;
GAIL2_Mx6p0_R10_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL2_Mx6p0_R10_01(:,3:22) = GAIL2_Mx6p0_R10_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{5} = GAIL2_Mx6p0_R10_01;
%%%%%%%%%%%%%%%%%%%  GAIL2_Mx6p0_R10_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL2_Mx6p0_R10_02 = haz2016.brnach{1}.C;
GAIL2_Mx6p0_R10_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL2_Mx6p0_R10_02(:,3:22) = GAIL2_Mx6p0_R10_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{6} = GAIL2_Mx6p0_R10_02;

%%%%%%%%%%%%%%%%%%%  GAIL2_Mx6p0_R20_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL2_Mx6p0_R20_01 = haz2016.brnach{1}.C;
GAIL2_Mx6p0_R20_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL2_Mx6p0_R20_01(:,3:22) = GAIL2_Mx6p0_R20_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{7} = GAIL2_Mx6p0_R20_01;
%%%%%%%%%%%%%%%%%%%  GAIL2_Mx6p0_R20_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL2_Mx6p0_R20_02 = haz2016.brnach{1}.C;
GAIL2_Mx6p0_R20_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL2');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL2_Mx6p0_R20_02(:,3:22) = GAIL2_Mx6p0_R20_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{8} = GAIL2_Mx6p0_R20_02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% GAIL5_Mx6p0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Mx6p0_R10_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL5_Mx6p0_R10_01 = haz2016.brnach{1}.C;
GAIL5_Mx6p0_R10_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL5_Mx6p0_R10_01(:,3:22) = GAIL5_Mx6p0_R10_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{9} = GAIL5_Mx6p0_R10_01;
%%%%%%%%%%%%%%%%%%%  GAIL5_Mx6p0_R10_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL5_Mx6p0_R10_02 = haz2016.brnach{1}.C;
GAIL5_Mx6p0_R10_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL5_Mx6p0_R10_02(:,3:22) = GAIL5_Mx6p0_R10_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{10} = GAIL5_Mx6p0_R10_02;

%%%%%%%%%%%%%%%%%%%  GAIL5_Mx6p0_R20_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL5_Mx6p0_R20_01 = haz2016.brnach{1}.C;
GAIL5_Mx6p0_R20_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1; 
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL5_Mx6p0_R20_01(:,3:22) = GAIL5_Mx6p0_R20_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{11} = GAIL5_Mx6p0_R20_01;
%%%%%%%%%%%%%%%%%%%  GAIL5_Mx6p0_R20_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GAIL5_Mx6p0_R20_02 = haz2016.brnach{1}.C;
GAIL5_Mx6p0_R20_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'Mx6p0'); tf2 = strcmp(GMPE,'GAIL5');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1; 
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    GAIL5_Mx6p0_R20_02(:,3:22) = GAIL5_Mx6p0_R20_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{12} = GAIL5_Mx6p0_R20_02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CEUS_MxNSH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  MxNSH_R10_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R10_01 = haz2016.brnach{1}.C;
CEUS_MxNSH_R10_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R10_01(:,3:22) = CEUS_MxNSH_R10_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{13} = CEUS_MxNSH_R10_01;
%%%%%%%%%%%%%%%%%%%  CEUS_MxNSH_R10_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R10_02 = haz2016.brnach{1}.C;
CEUS_MxNSH_R10_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R10');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R10_02(:,3:22) = CEUS_MxNSH_R10_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{14} = CEUS_MxNSH_R10_02;
%%%%%%%%%%%%%%%%%%%  CEUS_MxNSH_R10_LT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R10_LT = haz2016.brnach{1}.C;
CEUS_MxNSH_R10_LT(:,3:23) = 0;
for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf4 = strcmp(Cat,'LT');
    ft = tf1*tf2*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R10_LT(:,3:22) = CEUS_MxNSH_R10_LT(:,3:22) + tmp(:,3:22); 
    end
end
bran{15} = CEUS_MxNSH_R10_LT;
%%%%%%%%%%%%%%%%%%%  CEUS_MxNSH_R20_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R20_01 = haz2016.brnach{1}.C;
CEUS_MxNSH_R20_01(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'01');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R20_01(:,3:22) = CEUS_MxNSH_R20_01(:,3:22) + tmp(:,3:22); 
    end
end
bran{16} = CEUS_MxNSH_R20_01;
%%%%%%%%%%%%%%%%%%%  CEUS_MxNSH_R20_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R20_02 = haz2016.brnach{1}.C;
CEUS_MxNSH_R20_02(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'02');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R20_02(:,3:22) = CEUS_MxNSH_R20_02(:,3:22) + tmp(:,3:22); 
    end
end
bran{17} = CEUS_MxNSH_R20_02; 
%%%%%%%%%%%%%%%%%%%  CEUS_MxNSH_R20_36 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEUS_MxNSH_R20_36 = haz2016.brnach{1}.C;
CEUS_MxNSH_R20_36(:,3:23) = 0;

for bnum = 1:numel(haz2016.brnach)-1; 
    C = strsplit(haz2016.brnach{bnum}.N,'_'); 
    Cat =char(C(1)); R = char(C(3)); Max = char(C(4)); GMPE = char(C(5)); 
    tf1 = strcmp(Max,'MxNSH'); tf2 = strcmp(GMPE,'CEUS');tf3 = strcmp(R,'R20');tf4 = strcmp(Cat,'36');
    ft = tf1*tf2*tf3*tf4;
    if ft == 1
    CH = CH + 1;
    tmp = haz2016.brnach{bnum}.C;
    tmp(:,3:22) = tmp(:,3:22).*(tmp(1,23));
    CEUS_MxNSH_R20_36(:,3:22) = CEUS_MxNSH_R20_36(:,3:22) + tmp(:,3:22); 
    end
end
bran{18} = CEUS_MxNSH_R20_36;



% 
% %% constructing the mean hazard
% estmean = bran{1};
% for i = 2:length(bran);  
%     temp = bran{i};
%     estmean(:,3:22) = estmean(:,3:22) + temp(:,3:22); 
% end
% temp = haz2016.brnach{length(haz2016.brnach)}.C;
% temp(:,3:22) = temp(:,3:22).*temp(1,23);
% estmean(:,3:22) = estmean(:,3:22) + temp(:,3:22); 
% 
% 
% %% plotting hazard curves for a single point
% % lat = 39 ; lon = -84; ttl= 'Kentucky-Ohaio Boreder';
% % lat = 35.5 ; lon = -83; ttl= 'Eastern TN';
% % lat = 36.60 ; lon = -89.55; ttl= 'New-Madrid';
% % lat = 32.8 ; lon = -96.8; ttl= 'Dallas';
% lat = 35.45; lon = -97.5; ttl= 'Oklahoma City';
% 
% station = find(haz2016.USGSmean(:,1) == lat & haz2016.USGSmean(:,2) == lon) 
% 
% figure
% %% plotting individual branches
% for vv = 1:length(bran)
% an = bran{vv}; anR1 = an(station,3:22);
% anR11 = anR1(:) > 10e-11; s = nonzeros(anR11);
% h1 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'-k','LineWidth',1.5); hold on 
% end
% 
% %plotting USGS mean
% anR1 = haz2016.USGSmean(station,3:22);
% anR1 = anR1'; anR11 = anR1(:) > 10e-11;
% s = nonzeros(anR11); ed = length(s) + 2;
% h6 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'-m^','LineWidth',2.5,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','c')
% 
% %plotting estimated mean
% anR1 = estmean(station,3:22);
% anR1 = anR1'; anR11 = anR1(:) > 10e-5;
% s = nonzeros(anR11);ed = length(s) + 2;
% h7 = plot((haz2016.refPGA(1:length(s))),(anR1(1:length(s))),'--go','LineWidth',2.5,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','r')
% 
% set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');hold on 
% title(ttl);grid on;set(gca,'LineWidth',2);ylabel('Annual Frequency of Exceedances');
% xlabel('Peak Ground Acceleration (g)');xlim([0.005 2.5]);ylim([10e-6 100]);set(gca,'fontsize',21)
% 



%% Bootstraping 
Nsampling = 40;
for k = 1: Nsampling;
    pk = 100.*(k/Nsampling);
    X = [num2str(pk),' % of the Bootstraping is done.'];disp(X)
    
y = randi([2 length(bran)],1,floor(length(bran)*0.9));
ini = bran{y(1)};
for i = 6:length(y);
    temp = bran{y(i)};
    ini(:,3:22) = ini(:,3:22) + temp(:,3:22);
end  
% adding adaptive curve
temp = haz2016.brnach{length(haz2016.brnach)}.C;
temp(:,3:22) = temp(:,3:22).*temp(1,23);
ini(:,3:22) = ini(:,3:22) + temp(:,3:22);
    bootmodel{k} = ini;
end

informedNW= bran{y(1)}; informedNW(:,3:23)=0;mean = informedNW; stboot = informedNW;
for E = 1:length(informedNW);
    EP = (E./length(informedNW))*100;
    X = [num2str(EP),' % of the Standard Deviation calculation is done.'];disp(X)
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
PGAREF=[0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;0.26900E-01;...
0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;0.20300E+00;...
0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;0.15200E+01;...
0.22000E+01;0.33000E+01];


haz = zeros(length(bran),3);
ucr = zeros(length(bran),3);
ucrP = zeros(length(bran),3);
temp2 = bran{1};
for j = 1 : length(mean)
   H = mean(j,3:22); U = stboot(j,3:22);
   
   [c index1] = min(abs(H-0.0101));
   closestValues1 = H(index1); % finding the first minimum.

   if closestValues1 == 0.0101;
          yi = PGAREF(index1);
          ui = U(index1);
%           zi = (ui)./0.0101;
          
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
%    zi = (ui)./0.0101;   
%     ui = (U(index1)+U(index2))./2;   

   elseif closestValues1 < 0.0101
       if index1 == 1
          yi = PGAREF(index1);
          ui = U(index1);
%           zi = (ui)./0.0101;
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
%           zi = (ui)./0.0101;  
%           ui = (U(index1)+U(index2))./2;
       end
   end
 
   haz(j,1) = temp2(j,2); haz(j,2) = temp2(j,1); haz(j,3) = yi; 
   ucr(j,1) = temp2(j,2); ucr(j,2) = temp2(j,1); ucr(j,3) = ui;
   
end
haz(isnan(haz)) = 0; ucr(isnan(ucr)) = 0; 

delete('onePercent_1YrModel.pga.1pc1.txt')
fileID = fopen('onePercent_1YrModel.pga.1pc1.txt','w');
for k = 1 : length(haz)
fprintf(fileID,'%3.3f %2.3f %f\n',haz(k,:));
end
fclose(fileID);

delete('onePercent_1YrModel.pga.1pc2.txt')
fileID = fopen('onePercent_1YrModel.pga.1pc2.txt','w');
for k = 1 : length(ucr)
fprintf(fileID,'%3.3f %2.3f %f\n',ucr(k,:));
end
fclose(fileID);
