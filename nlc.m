function [c, ceq] = nlc(X) 

scale=[7.00  8.83886  10.27  16.95 0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1  0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1 10 10 10 16 89000 24600];

X=X.*scale;

%% AERODYNAMIC one-g

acc=1;
Cl_Cd=Q3D_Start(X,acc);          % CL/CD of the whole airplane
diff_clcd=Cl_Cd-X(44);

%% WING WEIGHT 

acc=2.5;
Wwing=Start(X,acc);           % Wing weight
Wa_w=89000 - 14467 - 24600;   % MTOW minus calculated wing weight (both included) minus fuel weight
MTOW=Wwing+Wa_w+X(46);        % PROBABLY weight of the wing icludes already both wings so it should be taken 1 time instead of twice
diff_mass = X(45)-MTOW;

%% Volumes and thickness
r=X(1);
k=X(1)-0.37*X(2);
t=X(3)-X(2);

a1=(r+k)*0.5*0.37*X(4);
a2=(k+t)*0.5*0.63*X(4);
area=2*(a1+a2);

xd=linspace(0,1,20);
T=zeros(1,3);
overlap=zeros(1,3);

for i=1:3
    Au=[X(5+(i-1)*12), X(6+(i-1)*12), X(7+(i-1)*12), X(8+(i-1)*12), X(9+(i-1)*12), X(10+(i-1)*12)];
    Al=[X(11+(i-1)*12), X(12+(i-1)*12), X(13+(i-1)*12), X(14+(i-1)*12), X(15+(i-1)*12), X(16+(i-1)*12)];
    [Xtu,Xtl] = D_airfoil2(Au,Al,xd);
    Xtu=(Xtu(:,2))';
    Xtl=(Xtl(:,2))';
    
    Thickness=0;
    
    % Volume
    for k=4:12      % takes into account front and rear spar
        Thickness=Thickness+Xtu(k)-Xtl(k);
    end
    
    T(i)=Thickness/9;
    
    % Not overlapping of the profiles
    diff_overlap=Xtu-Xtl;
    
    switch i
        case 1
            chord=r;
        case 2
            chord=k;
        case 3 
            chord=t;
    end            
    
    for j=1:length(diff_overlap)-2
        c((i-1)*18+j)=-diff_overlap(j+1)*chord;        
    end

end

Thickness_inner_wing = (T(1)*r+T(2)*k)/2;
Thickness_outer_wing = (T(2)*k+T(3)*t)/2;

Vol1=Thickness_inner_wing * ( ( 0.4 * X(1)+0.4*(X(1)-0.37*X(2) ) )*0.37 * X(4) ) / 2;
Vol2= Thickness_outer_wing *( ( 0.4 * (X(1)-0.37*X(2) )+ 0.4 * (X(3)-X(2) ) )*0.48 * X(4) ) / 2;

Vol=2*(Vol1+Vol2); % the volume of the fuel in both wings

Vol_eff=Vol*0.93;
rho_f=817.15;

W_fuel=Vol_eff*rho_f;
diff_fuel=-W_fuel+X(46); % because: -Wfuel +X(46) < 0

% Fuel consinstency

R = 5e6;            % range
C_t = 1.8639e-4;    % spec. fuel consumption
V = 0.76*295.2;     % cruise velocity

CL_CD = X(44);
W_TO_max = X(45);

Wfuel = (1-0.938*(1/(exp(R*C_t/(V*CL_CD)))))*W_TO_max;

%diff_f = Wfuel-X(46); % not necessary now, we use the ratio

% Wing loading

load0=89000/122.4;
load=X(45)/area;

% Constrains
L=length(c);

c(L+1)=diff_fuel;
c(L+2)=load-load0;

ceq = [Cl_Cd/X(44)-1, MTOW/X(45)-1,Wfuel/X(46)-1];
%ceq = [diff_clcd, MTOW/X(45)-1,Wfuel/X(46)-1];

C=[c,ceq];
fileID = fopen('constraint_violation.txt','a+');
fprintf(fileID,'%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n',C);
fclose(fileID);

fileID = fopen('design.txt','a+');
fprintf(fileID,'%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n',X);
fclose(fileID);

end

