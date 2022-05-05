
ss=specimenSymmetry('mmm');
cs=crystalSymmetry('m-3m');

sampschemes={}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a list of hkls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_ferrite = {Miller(1,1,0,cs),Miller(2,0,0,cs),Miller(2,1,1,cs)} ;

h_austenite = {Miller(1,1,1,cs),Miller(2,0,0,cs),Miller(2,2,0,cs)} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import DS Spiral Hex as PF scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname = '/Users/creuzige/Documents/NIST_Research/GitHub/Texture-Sampling-PhaseMeasurement-BiasErrors/Matlab/Spiral_Histograms/';
fname = [pname 'Spiral.txt'];

pf_single111_Spiral=PoleFigure.load({fname},{Miller({1,1,1}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single200_Spiral=PoleFigure.load({fname},{Miller({2,0,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single220_Spiral=PoleFigure.load({fname},{Miller({2,2,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);

pf_Spiral=[pf_single111_Spiral,pf_single200_Spiral,pf_single220_Spiral]

pf_Spiral.SS=ss;
pf_Spiral.CS=cs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Holden as PF scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname = '/Users/creuzige/Documents/NIST_Research/GitHub/Texture-Sampling-PhaseMeasurement-BiasErrors/Matlab/Spiral_Histograms/';
fname = [pname 'Holden.txt'];
pf_single111_Holden=PoleFigure.load({fname},{Miller({1,1,1}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single200_Holden=PoleFigure.load({fname},{Miller({2,0,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single220_Holden=PoleFigure.load({fname},{Miller({2,2,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);

pf_Holden=[pf_single111_Holden,pf_single200_Holden,pf_single220_Holden]

pf_Holden.SS=ss;
pf_Holden.CS=cs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Single Rizzie Hex as PF scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname = '/Users/creuzige/Documents/NIST_Research/GitHub/Texture-Sampling-PhaseMeasurement-BiasErrors/Matlab/Spiral_Histograms/';
fname = [pname 'Rizzie.txt'];
pf_single111_Rizzie=PoleFigure.load({fname},{Miller({1,1,1}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single200_Rizzie=PoleFigure.load({fname},{Miller({2,0,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single220_Rizzie=PoleFigure.load({fname},{Miller({2,2,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);

pf_Rizzie=[pf_single111_Rizzie,pf_single200_Rizzie,pf_single220_Rizzie]

pf_Rizzie.SS=ss;
pf_Rizzie.CS=cs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Matthies Hex as PF scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname = '/Users/creuzige/Documents/NIST_Research/GitHub/Texture-Sampling-PhaseMeasurement-BiasErrors/Matlab/Spiral_Histograms/';
fname = [pname 'Matthias.txt'];

pf_single111_Matthias=PoleFigure.load({fname},{Miller({1,1,1}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single200_Matthias=PoleFigure.load({fname},{Miller({2,0,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);
pf_single220_Matthias=PoleFigure.load({fname},{Miller({2,2,0}, cs)},'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);

pf_Matthias=PoleFigure.load({fname,fname,fname}, hList,'interface','generic','HEADER',1,'degree','ColumnNames',...
    {'polar angle','azimuth angle','intensity'},'Columns',[1 2 3]);

pf_Matthias=[pf_single111_Matthias,pf_single200_Matthias,pf_single220_Matthias]
pf_Matthias.SS=ss;
pf_Matthias.CS=cs;


%% Historgram of angles

edges=linspace(0,10,40);

omega_Spiral=angle(pf_Spiral.allR{1},pf_Spiral.allR{1}')/degree;
sortedomega_Spiral=sort(omega_Spiral);
% Second row, since first row is zeros
figure(1);
histogram(sortedomega_Spiral(2,:),edges,'Normalization','probability',...
    'FaceColor','k' )
xlabel('Angle to closest adjacent vector') 
ylabel('Relative probability') 
xtickformat('degrees')
median(sortedomega_Spiral(2,:))
exportgraphics(gca,'DSS_hist.png','Resolution',600)

omega_Holden=angle(pf_Holden.allR{1},pf_Holden.allR{1}')/degree;
sortedomega_Holden=sort(omega_Holden);
% Second row, since first row is zeros
figure(2);
histogram(sortedomega_Holden(2,:),edges,'Normalization','probability',...
    'FaceColor','k' )
xlabel('Angle to closest adjacent vector') 
ylabel('Relative probability') 
xtickformat('degrees')
median(sortedomega_Holden(2,:))
exportgraphics(gca,'Holden_hist.png','Resolution',600)

omega_Rizzie=angle(pf_Rizzie.allR{1},pf_Rizzie.allR{1}')/degree;
sortedomega_Rizzie=sort(omega_Rizzie);
% Second row, since first row is zeros
figure(3);
histogram(sortedomega_Rizzie(2,:),edges,'Normalization','probability',...
    'FaceColor','k' )
xlabel('Angle to closest adjacent vector') 
ylabel('Relative probability') 
xtickformat('degrees')
median(sortedomega_Rizzie(2,:))
exportgraphics(gca,'Rizzie_hist.png','Resolution',600)

omega_Matthias=angle(pf_Matthias.allR{1},pf_Matthias.allR{1}')/degree;
sortedomega_Matthias=sort(omega_Matthias);
% Second row, since first row is zeros
figure(4);
histogram(sortedomega_Matthias(2,:),edges,'Normalization','probability',...
    'FaceColor','k' )
xlabel('Angle to closest adjacent vector') 
ylabel('Relative probability') 
xtickformat('degrees')
median(sortedomega_Matthias(2,:))
exportgraphics(gca,'Matthias_hist.png','Resolution',600)


%% Calculate for a series of ODFs

disp('Define Orientations')

Cube = orientation.byMiller([0 0 1],[1 0 0],cs,ss);
Goss = orientation.byMiller([0 1 1],[1 0 0],cs,ss);
Shear = orientation.byMiller([0 0 1],[1 -1 0],cs,ss);
RGoss = orientation.byEuler(0*degree,90*degree,45*degree,cs,ss);

alpha1=orientation.byMiller([1 1 5],[1 -1 0],cs,ss);
alpha2=orientation.byMiller([1 1 3],[1 -1 0],cs,ss);
alpha3=orientation.byMiller([1 1 2],[1 -1 0],cs,ss);
alpha4=orientation.byMiller([2 2 3],[1 -1 0],cs,ss);

o554=orientation.byMiller([5 5 4],[-2 -2 5],cs,ss);

Copper = orientation.byEuler(90*degree,35.264*degree,45*degree,cs,ss);
CopperS=orientation.byEuler(74.49*degree,35.982*degree,54.2175*degree,cs,ss);
S = orientation.byEuler(58.98*degree,36.699*degree,63.435*degree,cs,ss);
BrassS=orientation.byEuler(47.122*degree,40.85*degree,76.718*degree,cs,ss);
Brass = orientation.byEuler(35.264*degree,45*degree,90*degree,cs,ss);

gamma1 = orientation.byEuler(0*degree,54.736*degree,45*degree,cs,ss);
gamma2 = orientation.byEuler(30*degree,54.736*degree,45*degree,cs,ss);

%% Start for loop of different orientaions

%HW=[20*degree]

for HW=[2.5*degree 5*degree 10*degree 15*degree 20*degree 25*degree 30*degree 35*degree 40*degree 45*degree 50*degree]

    psi = deLaValleePoussinKernel('HALFWIDTH',HW);
    fprintf('Halfwidth of %3.1f .\n',HW);


    RssqList=[];
    ODFEvalList=[];
    MinAbsDiffList=[];
    MeanAbsDiffList=[];
    MaxAbsDiffList=[];


    % Define ODFs
    %for i=7 %just make cube orientation plots
    for i=[1:20] %all of the textures

        %% uniform distributions - creates lots of them, but I couldn't find an elegant way to stop that from happening

        if i==1
            %uniform ODF
            disp('Uniform Austenite')
            bname='UniformA';
            phase ='austenite';
            odf=uniformODF(cs,ss);

        elseif i==2
            %uniform ODF
            disp('Uniform Ferrite')
            bname='UniformF';
            phase ='ferrite';
            odf=uniformODF(cs,ss);

            %% fiber distributions
        elseif i==3
            % ferrite alpha fiber 110 || RD
            bname='AlphaFiberF';
            phase ='ferrite';
            %gamma= fibre.gamma(cs);
            %ah=alphaFiber.h
            %ar=alphaFiber.r
            %odf = fibreODF(Miller(times(ah,o),cs),Miller(times(ar,o),cs),'halfwidth',HW)

            %h = Miller(1,1,0,cs);
            %r = xvector;
            %odf = fibreODF(h,r,ss,psi);
            sixth=(1.0/6.0)
            odf1 = (sixth)*unimodalODF(Shear,'halfwidth',HW,cs,ss,psi);
            odf2 = (sixth)*unimodalODF(alpha1,'halfwidth',HW,cs,ss,psi);
            odf3 = (sixth)*unimodalODF(alpha2,'halfwidth',HW,cs,ss,psi);
            odf4 = (sixth)*unimodalODF(alpha3,'halfwidth',HW,cs,ss,psi);
            odf5 = (sixth)*unimodalODF(alpha4,'halfwidth',HW,cs,ss,psi);
            odf6 = (sixth)*unimodalODF(gamma1,'halfwidth',HW,cs,ss,psi);

            odf=odf1+odf2+odf3+odf4+odf5+odf6;

        elseif i==4
            % ferrite fiber - gamma 111 || ND
            bname='GammaFiber2F';
            phase ='ferrite';

            odf1=.5*unimodalODF(gamma1,'halfwidth',HW,cs,ss,psi);
            odf2=.5*unimodalODF(gamma2,'halfwidth',HW,cs,ss,psi);
            odf=odf1+odf2;

        elseif i==5
            % austenite fiber - beta fiber, Brass -> Copper via S
            %Not working well using the defined beta fiber in mtex
            phase ='austenite';
            bname='BetaFiberA';

            odf1 = .2*unimodalODF(Brass,'halfwidth',HW,cs,ss);
            odf2 = .2*unimodalODF(S,'halfwidth',HW,cs,ss);
            odf3 = .2*unimodalODF(Copper,'halfwidth',HW,cs,ss);
            odf4 = .2*unimodalODF(CopperS,'halfwidth',HW,cs,ss);
            odf5 = .2*unimodalODF(BrassS,'halfwidth',HW,cs,ss);

            odf=odf1+odf2+odf3+odf4+odf5;

            %% single orientations

            % austenite single orientations
        elseif i==6
            bname='BrassA';
            phase ='austenite';
            odf = unimodalODF(Brass,'halfwidth',HW,cs,ss);

        elseif i==7
            bname='CubeA';
            phase ='austenite';
            odf = unimodalODF(Cube,'halfwidth',HW,cs,ss);

        elseif i==8
            bname='CopperA';
            phase ='austenite';
            odf = unimodalODF(Copper,'halfwidth',HW,cs,ss);

        elseif i==9
            bname='SA';
            phase ='austenite';
            odf = unimodalODF(S,'halfwidth',HW,cs,ss);

        elseif i==10
            bname='GossA';
            phase ='austenite';
            odf = unimodalODF(Goss,'halfwidth',HW,cs,ss);

            % ferrite single orientation
        elseif i==11
            bname='Gamma1F';
            phase ='ferrite';
            odf = unimodalODF(gamma1,'halfwidth',HW,cs,ss);

        elseif i==12
            bname='Gamma2F';
            phase ='ferrite';
            odf = unimodalODF(gamma2,'halfwidth',HW,cs,ss);

        elseif i==13
            bname='GossF';
            phase ='ferrite';
            odf = unimodalODF(Goss,'halfwidth',HW,cs,ss);

        elseif i==14
            bname='ShearF';
            phase ='ferrite';
            odf = unimodalODF(Shear,'halfwidth',HW,cs,ss);

        elseif i==15
            bname='O554F';
            phase ='ferrite';
            odf = unimodalODF(o554,'halfwidth',HW,cs,ss);

        elseif i==16
            bname='Alpha1F';
            phase ='ferrite';
            odf = unimodalODF(alpha1,'halfwidth',HW,cs,ss);

        elseif i==17
            bname='Alpha2F';
            phase ='ferrite';
            odf = unimodalODF(alpha2,'halfwidth',HW,cs,ss);

        elseif i==18
            bname='Alpha3F';
            phase ='ferrite';
            odf = unimodalODF(alpha3,'halfwidth',HW,cs,ss);

        elseif i==19
            bname='Alpha4F';
            phase ='ferrite';
            odf = unimodalODF(alpha4,'halfwidth',HW,cs,ss);

        elseif i==20
            bname='RGossF';
            phase ='ferrite';
            odf = unimodalODF(RGoss,'halfwidth',HW,cs,ss);
        end



        %%
        SchemeList={pf_Spiral.r,pf_Holden.r,pf_Rizzie.r,pf_Matthias.r};

        %pf_rc = calcPoleFigure(odf_Al,hList,pf_dup_Rizzie.r)
        %odf_rc = calcODF(pf_rc, cs,ss, 'halfwidth',5*degree, 'resolution',5*degree)

        for j=1:length(SchemeList)

            if strcmp('ferrite',phase)
                hList=h_ferrite;
            elseif strcmp('austenite',phase)
                hList=h_austenite;
            else
                disp('Incorrect Phase')
            end

            pf_rc = calcPoleFigure(odf,hList,SchemeList{j});
            SchemeList{j};

            % pf_rc = makePartial(pf_rc); % Comment or uncomment depending on your need
            odf_rc = calcODF(pf_rc, cs,ss, 'halfwidth',5*degree, 'resolution',5*degree);

            odfdiff= odf - odf_rc;

            ODFEvalList = eval(odfdiff, OrList);
            MinAbsDiffList= [MinAbsDiffList min(abs(ODFEvalList))];
            MeanAbsDiffList= [MeanAbsDiffList mean(abs(ODFEvalList))];
            MaxAbsDiffList= [MaxAbsDiffList max(abs(ODFEvalList))];
            RssqList= [RssqList rssq(ODFEvalList)/length(ODFEvalList)];

            %         figure(100*i+j)
            %         plot(odfdiff,'phi2','sections',18,'projection','plain','minmax','off',cs,ss);
            %         CLim(gcm,[-1, 1]);mtexColorbar;
            %         run('ColorMapCoolWarm.m')
        end

    end

    
    A=reshape(MeanAbsDiffList,[4,20])
    writematrix(A,strcat('MeanAbsDiff-HW', num2str(HW,3),'.csv'),'Delimiter','tab');
    B=sum(reshape(MeanAbsDiffList,[4,20]),2)/20
    writematrix(B,strcat('MeanAbsDiff-HW-Ave', num2str(HW,3),'.csv'),'Delimiter','tab');
end


