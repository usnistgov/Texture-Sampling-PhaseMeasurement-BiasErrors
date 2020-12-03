%% Read and analyze data from: 
% Micromechanical Response Quantification using High-energy X-rays during Phase Transformations in Additively Manufactured 17-4 Stainless Steel

% https://data.mendeley.com/datasets/3mddz99wsr/1

% plotting convention:
plotx2east

%% Load data from .mat file

%load('./ChessData/Wrought_Specimen_Pole_Figure_Data.mat')
% 10 data points
% fully martensitic, ignore all austenite points
% Starts along Z direction

load('./ChessData/AM_Specimen_Pole_Figure_Data.mat')
% 12 datapoints on loading
% starts along X direction

%choose which load step to use 1 - undeoformed, 12 - end of test
Load_step=1

%The ?nS? is the list of unit vectors for each location on the pole figure in x, y, and z
%The ?i_list? is the list of the associated intensity values for each point on the pole figure.
% There are 6 total sets, represented by the last number in ?val(:,:,1)?, each set is associated with a diffraction ring. With 1 being the inner most diffraction ring. According to our indexing per https://ars.els-cdn.com/content/image/1-s2.0-S0921509319306203-gr3_lrg.jpg, the rings should correspond to:
% 1: austenite {111}
% 2: martensite {110}
% 3: austenite {200}
% 4: martensite {200}
% 5: austenite {220}
% 6: martensite {211}
%  

%% Create folders to output
A_name = strcat("Chess_AM_A_step",num2str(Load_step))
M_name = strcat("Chess_AM_M_step",num2str(Load_step))
mkdir("ChessData", A_name)
mkdir("ChessData", M_name)

%% convert to mtex format

% Use PoleFigure Constructor, not a supported data type
% https://mtex-toolbox.github.io/PoleFigure.PoleFigure.html

% crystal symmetry, from Phan, Thien Q., Felix H. Kim, and Darren C. Pagan. 2019.
%?Micromechanical Response Quantification Using High-Energy X-Rays during Phase 
% Transformations in Additively Manufactured 17-4 Stainless Steel.? 
% Materials Science and Engineering: A 759 (June): 565?73. https://doi.org/10.1016/j.msea.2019.05.017.

CrySym = {... 
  crystalSymmetry('m-3m', [3.520 3.520 3.520], 'mineral', 'austenite', 'color', 'light blue'),
  crystalSymmetry('m-3m', [2.934 2.934 2.934], 'mineral', 'martensite', 'color', 'green')};

% Choose between Orthotropic (orthorhombic) sample symmetry 
% and Triclinic (no) sample symmetry
% Triclinic set currently
SpecSym={specimenSymmetry('-1'),
    specimenSymmetry('mmm')};

h_austenite=Miller({1,1,1},{2,0,0},{2,2,0}, CrySym{1});
h1=Miller({1,1,1}, CrySym{1});
h3=Miller({2,0,0}, CrySym{1});
h5=Miller({2,2,0}, CrySym{1});

h_martensite=Miller({1,1,0},{2,0,0},{2,1,1}, CrySym{2});
h2=Miller({1,1,0}, CrySym{2});
h4=Miller({2,0,0}, CrySym{2});
h6=Miller({2,1,1}, CrySym{2});

%% Create the diffraction vectors
%nS(1,:,1) % returns the first row of x,y,z values for pole figure 1
% need to take transpose for vector3d to work properly
r1vectors=vector3d(nS(:,:,1).');
r2vectors=vector3d(nS(:,:,2).');
r3vectors=vector3d(nS(:,:,3).');
r4vectors=vector3d(nS(:,:,4).');
r5vectors=vector3d(nS(:,:,5).');
r6vectors=vector3d(nS(:,:,6).');

%% Plot the vectors
% % test function to check plotting conventions and data import
% 
% % 36*72= 2592 points.  But length is 2736.. =>38
% % Looks like all points greater than 2592 are x,y,z=0,0,1
% figure(1)
% n=0;% (0 to 35 -> 36 rotation steps around y)
% range=[n*72+1:n*72+72]; % 72 segments per line, (360/5° azimuth) seems to start with a slight rotation
% r1vectors(range);
% plot(r1vectors(range),'marker','s','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r')
% 
% n=17% (0 to 35 -> 36 rotation steps around y)
% range=[n*72+1:n*72+72]; % 72 segments per line, (360/5° azimuth) seems to start with a slight rotation
% r1vectors(range);
% annotate(r1vectors(range),'marker','o','MarkerSize',6,'MarkerEdgeColor','g','MarkerFaceColor','g','antipodal')
% 
% n=35;% (0 to 35 -> 36 rotation steps around y)
% range=[n*72+1:n*72+72]; % 72 segments per line, (360/5° azimuth) seems to start with a slight rotation
% r1vectors(range);
% annotate(r1vectors(range),'marker','x','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b','antipodal')
% 
% annotate(xvector,'label',{'X'});
% annotate(yvector,'label',{'Y'});
% annotate(zvector,'label',{'Z'});
% annotate(-xvector,'label',{'-X'});
% annotate(-yvector,'label',{'-Y'});
% annotate(-zvector,'label',{'-Z'});
% %hold off
% % shows an odd extra point on the upper hemisphere

%% plot vectors as filled contour
% % Looks similar to the oversampling plot in paper
% figure(2)
% plot(r1vectors([1:2592]),'contourf','antipodal','smooth' ); CLim(gcm,[0, 8]); mtexColorbar; 
% run('ColorMap8.m');


%% Read in the intensities

I_pf1=i_list(Load_step,:,1);
I_pf2=i_list(Load_step,:,2);
I_pf3=i_list(Load_step,:,3);
I_pf4=i_list(Load_step,:,4);
I_pf5=i_list(Load_step,:,5);
I_pf6=i_list(Load_step,:,6);

%% Create pole figures

pf1 = PoleFigure(h1,r1vectors,I_pf1);
pf2 = PoleFigure(h2,r2vectors,I_pf2);
pf3 = PoleFigure(h3,r3vectors,I_pf3);
pf4 = PoleFigure(h4,r4vectors,I_pf4);
pf5 = PoleFigure(h5,r5vectors,I_pf5);
pf6 = PoleFigure(h6,r6vectors,I_pf6);

pf_austenite=[pf1,pf3,pf5];
pf_martensite=[pf2,pf4,pf6];

figure(11); plot(pf_austenite, 'MarkerSize',8);
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AustenitePFs_counts",".png"))

figure(12); plot(pf_martensite, 'MarkerSize',8);
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensitePFs_counts",".png"))

figure(13); plot(pf_austenite.normalize, 'MarkerSize',8);CLim(gcm,[0, 4]); mtexColorbar;run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AustenitePFs",".png"))

figure(14); plot(pf_martensite.normalize,'MarkerSize',8);CLim(gcm,[0, 4]); mtexColorbar;run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensitePFs",".png"))


%% Calculate summed intensities for each pole figure
% basing code from mean.m in PoleFigureClass
% Not currently used

% mA = zeros(1,pf_austenite.numPF);
% for i = 1:pf_austenite.numPF
%      
%   mA(i) = sum(pf_austenite.allI{i});
%   
% end
% 
% A_IntName = strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","Chess_AM_Austenite_Int_step",num2str(Load_step),".txt")
% writematrix(mA,A_IntName)
% 
% 
% mM = zeros(1,pf_martensite.numPF);
% for i = 1:pf_martensite.numPF
%      
%   mM(i) = sum(pf_martensite.allI{i});
%   
% end
% 
% M_IntName = strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","Chess_AM_Martensite_Int_step",num2str(Load_step),".txt")
% writematrix(mM,M_IntName)


%% Plot normalized and contoured
% by inspection, these look similar to values in the paper
figure(15); plot(normalize(pf_austenite),'contourf','antipodal','smooth');
CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AustenitePFs_contour",".png"))

figure(16); plot(normalize(pf_martensite),'contourf','antipodal','smooth');
CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensitePFs_contour",".png"))


%% ODFs

% See Issue #506, need to normalize pole figures
odf_austenite = calcODF(pf_austenite.normalize, CrySym{1},SpecSym{1}, 'halfwidth',4*degree, 'resolution',4*degree, 'iterMin', 20)
odf_martensite = calcODF(pf_martensite.normalize, CrySym{2},SpecSym{1}, 'halfwidth',4*degree, 'resolution',4*degree, 'iterMin', 20)

figure(21); plot(odf_austenite,'phi2','sections',18,'projection','plain','minmax', 'off',CrySym{1},SpecSym{1});CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AusteniteODF",".png"))

figure(22); plot(odf_martensite,'phi2','sections',18,'projection','plain','minmax', 'off',CrySym{2},SpecSym{1});CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensiteODF",".png"))

%% Export ODFs in .maa format for MAUD
% Not currently used

% 
% %disp('Output ODFs to .maa and mtex')
% disp('Output ODFs to .maa')
% 
% arg= cell(1,2);
% arg{1}= 'resolution';
% arg{2}= degree ;
% 
% csplot = crystalSymmetry('m-3m');
% ssplot = specimenSymmetry('triclinic');
% S3G = getClass(arg,'SO3Grid',regularSO3Grid(csplot,ssplot,arg));
% 
% %SavePaths
% savepath='XPCFiles';
% mkdir(savepath)
% %bname='CHESS-AM-step1-A';
% %v = eval(odf_austenite,S3G,arg);
% bname='CHESS-AM-step1-F';
% v = eval(odf_martensite,S3G,arg);
% 
% file= fopen( strcat(savepath,'/',bname,'.maa'),'w');
% 
% fprintf(file, 'title \n');
% fprintf(file, '7 5.0 \n');
% fprintf(file, '2.86 2.86 2.86 90.00 90.00 90.00 \n');
% 
% %need to add one to the value for Beartex - nope, scale of pole figures 
% %is still messed up
% %between the two programs.  Do we need to rescale?
% 
% for i=1:size(v,3) 
%     for j=1:size(v,2)
%        for k=1:size(v,1) 
%            fprintf(file,'%2.1f ',v(k,j,i));
%            
%        end
%        fprintf(file,'%2.1f ',v(1,j,i)); %repeat for extra term
%        fprintf(file,'  \n');
%     end
%     
%     fprintf(file,'  \n');
% end
% 
% %repeat for extra data block
%     for j=1:size(v,2)
%        for k=1:size(v,1) 
%            fprintf(file,'%2.1f ',v(k,j,1));
%            
%        end
%        fprintf(file,'%2.1f ',v(1,j,1)); 
%        fprintf(file,'  \n');
%     end
%         fprintf(file,'  \n');


%% Recalculated pole figures
r_save = regularS2Grid('resolution',5*degree);

% Choose list of 
%h_austenite_save=Miller({1,1,1},{2,0,0},{2,2,0},{3,1,1},{2,2,2},{4,0,0},{3,3,1},{4,2,0},{4,2,2},{5,1,1},{3,3,3}, CrySym{1});
h_austenite_save=Miller({1,1,1},{2,0,0},{2,2,0}, CrySym{1});
%h_austenite_save=Miller({1,1,1}, CrySym{1});

%h_martensite_save=Miller({1,1,0},{2,0,0},{2,1,1},{3,1,0},{2,2,2},{3,2,1},{4,0,0}, CrySym{2})
h_martensite_save=Miller({1,1,0},{2,0,0},{2,1,1}, CrySym{2})
%h_martensite_save=Miller({1,1,0},{2,0,0},{2,1,1}, CrySym{2})

pdf_rc_a = calcPoleFigure(odf_austenite,h_austenite_save,r_save)
pdf_rc_m = calcPoleFigure(odf_martensite,h_martensite_save,r_save)

figure(31); plot(pdf_rc_a,'contourf','antipodal','smooth');
CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AustenitePFs_RC_contour",".png"))

figure(32); plot(pdf_rc_m,'contourf','antipodal','smooth');
CLim(gcm,[0, 4]); mtexColorbar; 
run('ColorMap4.m');
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensitePFs_RC_contour",".png"))


%% Export recalculated pole figures

%mtex seems unhappy with strings, wants character arrays...

savefile_A=strcat("ChessData/",A_name,"/",A_name);
export(pdf_rc_a,convertStringsToChars(savefile_A),'degree')

savefile_M=strcat("ChessData/",M_name,"/",M_name);
export(pdf_rc_m,convertStringsToChars(savefile_M),'degree')
%% Difference Pole Figures

figure(33);
plotDiff(pf_austenite.normalize,odf_austenite,'contourf','antipodal','smooth' )
mtexColorbar;
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_A_step",num2str(Load_step),"/","AustenitePFs_diff",".png"))


figure(34);
plotDiff(pf_martensite.normalize,odf_martensite,'contourf','antipodal','smooth' )
mtexColorbar;
print('-dpng', '-r300', strcat("ChessData/","Chess_AM_M_step",num2str(Load_step),"/","MartensitePFs_diff",".png"))
      
