function [Wwing]=Start(X,acc)

filename = 'test';

%% geometry and load calculations

DisplayOption = 0; 

%% INITIATOR (not an external function)

kink_perc=0.37;

% Design Weight

I.Weight.MTOW =  X(45);         % Maximum Take Off Weight        
I.Weight.ZFW = X(45)-X(46);     % Zero Fuel Weight     
I.n_max = 2.5;

% Wing Geometry

r=X(1);
k=X(1)-0.37*X(2);
t=X(3)-X(2);
a1=(r+k)*0.5*0.37*X(4);
a2=(k+t)*0.5*0.63*X(4);
area=(a1+a2);           % one wing area

I.Wing(1).Area= 2*area; % two wings area (total area)
I.Wing(1).Span= 2*X(4); % total span (both wings)   
I.Wing(1).SectionNumber = 3;
I.Wing(1).AirfoilNumber = 3;

% GENERATE .dat FILE FROM X(i) %%%
xd=linspace(0,1,20);

for i=1:3
    Au=[X(5+(i-1)*12), X(6+(i-1)*12), X(7+(i-1)*12), X(8+(i-1)*12), X(9+(i-1)*12), X(10+(i-1)*12)];
    Al=[X(11+(i-1)*12), X(12+(i-1)*12), X(13+(i-1)*12), X(14+(i-1)*12), X(15+(i-1)*12), X(16+(i-1)*12)];
    [Xtu,Xtl] = D_airfoil2(Au,Al,xd);
    
    [S1,S2]=size(Xtu);
    transition=zeros(S1,S2);
    
    for k=1:S1
        transition(k,:)=Xtu(S1-k+1,:);
    end
    
    Xtu=transition;
    
    switch i
        case 1
             fid = fopen('root.dat', 'w');   
        case 2
             fid = fopen('kink.dat', 'w');   
        case 3
             fid = fopen('tip.dat', 'w');  
    end
    
    xy_coor=[Xtu;Xtl];
    fprintf(fid,'%g %g\r\n',xy_coor');
    fclose(fid);
    
end

I.Wing(1).AirfoilName = {'root' 'kink' 'tip'};
I.Wing(1).AirfoilPosition = [0 kink_perc 1];

% root

I.Wing(1).WingSection.Chord=X(1);
I.Wing(1).WingSection.Xle=0;
I.Wing(1).WingSection.Yle=0;
I.Wing(1).WingSection.Zle=0;
I.Wing(1).WingSection.FrontSparPosition=  0.18 ;
I.Wing(1).WingSection.RearSparPosition= 0.58  ;


% kink

I.Wing(2).WingSection.Chord=X(1)-kink_perc*X(2);
I.Wing(2).WingSection.Xle=X(2)*kink_perc; % from the triangle relations
I.Wing(2).WingSection.Yle=X(4)*kink_perc; % from the triangle relations
I.Wing(2).WingSection.Zle=0;
I.Wing(2).WingSection.FrontSparPosition= 0.18  ;
I.Wing(2).WingSection.RearSparPosition= 0.58 ;

% tip

I.Wing(3).WingSection.Chord=X(3)-X(2);
I.Wing(3).WingSection.Xle=X(2);
I.Wing(3).WingSection.Yle=X(4);
I.Wing(3).WingSection.Zle=0;
I.Wing(3).WingSection.FrontSparPosition= 0.18  ; %lucky number! instruction: 15-20%
I.Wing(3).WingSection.RearSparPosition = 0.58 ;  %instruction: 55-60%


% fuel tanks geometry

I.WingFuelTank.Ystart = 0.1 ; % due to landing gear system for example
I.WingFuelTank.Yend = 0.85;   % from the instruction


% Power plant and landing gear and wing atructure
   
I.PP(1).WingEngineNumber   =   1;  
I.PP(1).EnginePosition = 0.34;      % position of the engine 
I.PP(1).EngineWeight = 2380*1.6;    % (CFM56-5B3) correlation used: W_propulsion=W_engine,dry*1.6


% Material and Structure

I.Material.Wing.UpperPanel.E = 70.1e9;  %upper panel
I.Material.Wing.UpperPanel.rho =2800;
I.Material.Wing.UpperPanel.Sigma_tensile = 295e6;
I.Material.Wing.UpperPanel.Sigma_compressive = 295e6;

I.Material.Wing.LowerPanel.E = 70.1e9;  %lower panel
I.Material.Wing.LowerPanel.rho = 2800;
I.Material.Wing.LowerPanel.Sigma_tensile =295e6;
I.Material.Wing.LowerPanel.Sigma_compressive = 295e6;

I.Material.Wing.FrontSpar.E = 70.1e9; %front spar
I.Material.Wing.FrontSpar.rho = 2800;
I.Material.Wing.FrontSpar.Sigma_tensile = 295e6;
I.Material.Wing.FrontSpar.Sigma_compressive = 295e6;

I.Material.Wing.RearSpar.E = 70.1e9;    %rear spar
I.Material.Wing.RearSpar.rho = 2800;
I.Material.Wing.RearSpar.Sigma_tensile =295e6;
I.Material.Wing.RearSpar.Sigma_compressive = 295e6;


I.Structure.Wing.UpperPanelEfficiency = 1.02 ;  %Assumed - Integral Zed Stiffners - table: F=1.02
I.Structure.Wing.RibPitch =0.5;


%% Aerodynamic analysis

Res=Q3D_Start(X,acc);  

q=0.5*0.364*(295.2*0.82)^2;

MAC=((r-2/3*((r-k)*(0.5*r+k))/(r+k))*a1+(k-2/3*((k-t)*(0.5*k+t))/(k+t))*a2)/(a1+a2); % mean aerodynamic chord

AS.Y = linspace(0,1,20);
AS.L = interp1(Res.Wing.Yst,Res.Wing.ccl*q,AS.Y*I.Wing(1).Span/2,'spline'); %lift distribution
AS.T = interp1(Res.Wing.Yst,Res.Wing.cm_c4.*Res.Wing.chord*q*MAC,AS.Y*I.Wing(1).Span/2,'spline'); % pitching moment distribution


%% creating init file

fid = fopen([filename '.init'], 'wt');   
fprintf(fid,'%g %g\n',I.Weight.MTOW,I.Weight.ZFW);
fprintf(fid,'%g\n',I.n_max);
fprintf(fid,'%g %g %g %g\n',I.Wing(1).Area,I.Wing(1).Span,I.Wing(1).SectionNumber,I.Wing(1).AirfoilNumber);

for i=1:length(I.Wing(1).AirfoilName)
    fprintf(fid,'%g %s\n',I.Wing(1).AirfoilPosition(i),I.Wing(1).AirfoilName{i});
end

for i=1:I.Wing(1).SectionNumber
   fprintf(fid,'%g %g %g %g %g %g\n',I.Wing(i).WingSection.Chord,I.Wing(i).WingSection.Xle,I.Wing(i).WingSection.Yle,I.Wing(i).WingSection.Zle,I.Wing(i).WingSection.FrontSparPosition,I.Wing(i).WingSection.RearSparPosition); 
end

fprintf(fid,'%g %g\n',I.WingFuelTank.Ystart,I.WingFuelTank.Yend);
fprintf(fid,'%g\n',I.PP(1).WingEngineNumber);

for i = 1:I.PP(1).WingEngineNumber
    fprintf(fid,'%g %g\n',I.PP(i).EnginePosition,I.PP(i).EngineWeight);
end

fprintf(fid,'%g %g %g %g\n',I.Material.Wing.UpperPanel.E,I.Material.Wing.UpperPanel.rho,I.Material.Wing.UpperPanel.Sigma_tensile,I.Material.Wing.UpperPanel.Sigma_compressive);
fprintf(fid,'%g %g %g %g\n',I.Material.Wing.LowerPanel.E,I.Material.Wing.LowerPanel.rho,I.Material.Wing.LowerPanel.Sigma_tensile,I.Material.Wing.LowerPanel.Sigma_compressive);
fprintf(fid,'%g %g %g %g\n',I.Material.Wing.FrontSpar.E,I.Material.Wing.FrontSpar.rho,I.Material.Wing.FrontSpar.Sigma_tensile,I.Material.Wing.FrontSpar.Sigma_compressive);
fprintf(fid,'%g %g %g %g\n',I.Material.Wing.RearSpar.E,I.Material.Wing.RearSpar.rho,I.Material.Wing.RearSpar.Sigma_tensile,I.Material.Wing.RearSpar.Sigma_compressive);

fprintf(fid,'%g %g\n',I.Structure.Wing.UpperPanelEfficiency,I.Structure.Wing.RibPitch);
fprintf(fid,'%g\n',DisplayOption);

fclose(fid);

%% creating load file

fid = fopen([filename '.load'], 'wt');  

for i=1:length(AS.Y)
    fprintf(fid,'%g %g %g\n',AS.Y(i),AS.L(i),AS.T(i));
end

fclose(fid);


%% EMWET

EMWET test

%% reading output

fid     = fopen('test.weight', 'r');
OUT = textscan(fid, '%s'); 
fclose(fid);

out = OUT{1};
Wwing = str2double (out(4));

end