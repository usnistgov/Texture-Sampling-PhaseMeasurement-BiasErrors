% Mtex color maps

% Needed to produce test ODF
%clc
%clear
%close all
%% set colomap with blue -1, red +1, white in middle
mtexColorMap([ .05 .05 1;.1 .1 1; .2 .2 1; .3 .3 1; .4 .4 1; .5 .5 1; .6 .6 1;.7 .7 1;.8 .8 1; .9 .9 1;1 .9 .9;...
1 .8 .8;1 .7 .7; 1 .6 .6; 1 .5 .5; 1 .4 .4; 1 .3 .3;1 .2 .2;1 .1 .1;1 .05 .05 ])
cmap = colormap;
setMTEXpref('defaultColorMap',cmap);
setMTEXpref('FontSize',20);
%% test ODF
%{
savepath='.';

HW=20*degree;

cs = crystalSymmetry('m-3m');
ss = specimenSymmetry('orthorhombic');
Cube = orientation('Miller',[0 0 1],[1 0 0],cs,ss);
odf = unimodalODF(Cube,'halfwidth',HW,cs,ss);
figure; plot(odf,'phi2',[45]*degree,'projection','plain','minmax', 'off',cs,ss);CLim(gcm,[0, 4]);mtexColorbar;
export_fig(strcat(savepath,'/','ColorMapExample-phi2-45ODF.tiff'),'-r300') 
%export_fig(strcat('ColorMapExample-phi2-45ODF.tiff'),'-r300') 
%}

%% reset to default Mtex colomaps
%setMTEXpref('defaultColorMap',WhiteJetColorpMap);
