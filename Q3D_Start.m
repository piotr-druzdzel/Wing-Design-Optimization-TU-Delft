function Res=Q3D_Start(X,acc)

if acc==1                                 % acceleration 1g - cruise
     vis_inv=1;
     velocity=0.76*295.2;
     W = sqrt(X(45)*(X(45)-X(46)));       % weight given byt he formula in the assignment
     
elseif acc==2.5
     vis_inv=0;                           % pull up manouvre - 2.5 g
     velocity=0.82*295.2;
     W = X(45);                           % half of the MTOW
end

%% Aerodynamic solver setting

kink_perc=0.37;                                           % Calculated from the drawing

% Wing planform geometry 
%                 x                       y                 z    chord(m)               twist angle (deg) 
AC.Wing.Geom = [  0                       0                 0    X(1)                   X(41);  % shouldnt the initial incidence angle be higher than 0? (We have  a swept wing - lets assume 4)
                  kink_perc*X(2)+0.01     kink_perc*X(4)    0    X(1)-kink_perc*X(2)    X(42);  % at the beginning obtained from linear interpolation - then becomes a variable, one centimeter added at the trailing edge to avoid 0 value for the EMWET
                  X(2)                    X(4)              0    X(3)-X(2)              X(43)]; % Lets say 1 degree after torsion in cruise (

% Wing incidence angle (degree)
AC.Wing.inc  = 0;   
                  
% Airfoil coefficients input matrix
%                  | ->                    upper curve coeff.                <-|   | ->                    lower curve coeff.                    <-| 
AC.Wing.Airfoils   = [  X(5)      X(6)      X(7)      X(8)      X(9)      X(10)     X(11)     X(12)     X(13)     X(14)     X(15)     X(16);
                        X(17)     X(18)     X(19)     X(20)     X(21)     X(22)     X(23)     X(24)     X(25)     X(26)     X(27)     X(28);
                        X(29)     X(30)     X(31)     X(32)     X(33)     X(34)     X(35)     X(36)     X(37)     X(38)     X(39)     X(40)];
                       
AC.Wing.eta = [0;kink_perc;1];  % Spanwise location of the airfoil sections

AC.Visc  = vis_inv;             % 1 for viscous analysis and 0 for inviscid


% mean aerodynamic chord (from NASA website: nasascale.org)
r=X(1);
k=X(1)-0.37*X(2);
t=X(3)-X(2);
a1=(r+k)*0.5*0.37*X(4);
a2=(k+t)*0.5*0.63*X(4);

MAC =((r-2/3*((r-k)*(0.5*r+k))/(r+k))*a1+(k-2/3*((k-t)*(0.5*k+t))/(k+t))*a2)/(a1+a2); % mean aerodynamic chord

area=2*(a1+a2);
    
% Flight Condition
AC.Aero.V     = velocity;                                         % flight speed (m/s)
AC.Aero.rho   = 0.364;                                            % air density  (kg/m3)
AC.Aero.alt   = 11000;                                            % flight altitude (m)
AC.Aero.Re    = AC.Aero.rho*AC.Aero.V* MAC/1.44e-5 ;              % reynolds number (bqased on mean aerodynamic chord)
AC.Aero.M     = AC.Aero.V/295.2;                                  % flight Mach number 
AC.Aero.CL    = acc*9.8*W/(0.5*AC.Aero.rho*velocity^2*area);      % lift coefficient - comment this line to run the code for given alpha%

%% 

Res = Q3D_solver(AC);

area0=122.4;               % initial area of the wing from the excel given
CDa_w_in=0.6617/16-0.0248; % from the formula given, after obtaining from the Q3D: Cl wing, Cd wing and assuming L/D = 16

if acc==1
    CDa_w= CDa_w_in * area0/area;
    if Res.CDwing<0 || isnan(Res.CDwing)==1
        Res.CDwing=10000;
    end
    Res=Res.CLwing/(Res.CDwing+CDa_w);
end

end