clear all;

% ----------------- Ustawienia wykresów
LineWidth = 0.7;
gap = [0.12,0.04];
margin_h = [0.13,0.1];
margin_v = [0.05,0.05];
% ----------------- 

alfa = 15;
C = 0.2;

TH = 71;
TC = 22;
TD = 34;

FC = 34;
FH = 22;
FD = 15;

tauH = 130;
tauC = 100;

T = 39.7183;
h = 22.4044;
hL = h;
TL = T;

FH0=FH;
FC0=FC;
FD0=FD;
h0 = h;
T0 = T;

V = C*h^3;
TV = T*V;

% Macierze A,B,C,D - Równań stanu
A_m = [-(4*FC0+4*FD0+4*FH0-3*alfa*sqrt(h0))/(6*C*h0^3),... 
    0;...
    -3*(FH0*TH+FC0*TC+FD0*TD-FH0*T0-FC0*T0-FD0*T0)/(C*h0^4),...
    -(FH0+FC0+FD0)/(C*h0^3)]
B_m = [1/(3*C*h0^2),1/(3*C*h0^2),1/(3*C*h0^2);...
    (TH-T0)/(C*h0^3),(TC-T0)/(C*h0^3),(TD-T0)/(C*h0^3)]
C_m = [1,0;0,1]
D_m = [0,0,0;0,0,0]

% Transmitancja układu
sys = ss(A_m,B_m,C_m,D_m);  
G_s = tf(sys);

% Transmitancje układu ciągłego
G11 = tf(G_s.Numerator(1,1),G_s.Denominator(1,1),'inputdelay',130);
G12 = tf(G_s.Numerator(1,2),G_s.Denominator(1,2),'inputdelay',100);
G21 = tf(G_s.Numerator(2,1),G_s.Denominator(2,1),'inputdelay',130);
G22 = tf(G_s.Numerator(2,2),G_s.Denominator(2,2),'inputdelay',100);
G13 = tf(G_s.Numerator(1,3),G_s.Denominator(1,3));
G23 = tf(G_s.Numerator(2,3),G_s.Denominator(2,3));

Tp = 10;    % Okres próbkowania
N = 60;     % Horyzont predykcji
Nu = 20;    % Horyzont sterowania
ny = 2;     % ilość wyjść
nu = 2;     % ilość wejść
D = 60;     % Horyzont dynamiki
Dz = 40;    % Horyzont dynamiki zakłócenia
psi = diag(repmat([1 1],1,N));
lambda = diag(repmat([1 1],1,Nu));

G_s = [G11 G12; G21 G22];
G_z = c2d(G_s, Tp);

S = DMCstepmatrices(Tp, D, G_z);
[M,Mp] = DMCmatrices(S,N,Nu);
K = inv(M'*psi*M+lambda)*M'*psi;
K1 = K(1:nu,:);

% Odpowiedź skokowa modelu zakłócenia
G_s_dist = [G13; G23];
G_z_dist = c2d(G_s_dist, Tp);
S_dist = DMCstepmatrices(Tp, Dz, G_z_dist);
MPZ = DMCdistMatrix(S_dist, N, Nu);

% Czas próbkowania i czas symulacji
T_max = 10000;
steps = round(T_max/Tp);

% Dyskretne równania stanu A_d, B_d, C_d, D_d
A_d = expm(A_m*Tp)
B_d = (A_m^(-1))*(A_d-eye(2,2))*B_m
C_d = C_m;
D_d = D_m;

du_max = [4; 4];
u_max = [400; 400];
u_min = [-400; -400];

y_zad = repmat([h0; T0],1,steps+N);
y_zad(1,100:300) = h0 + 2;
y_zad(2,200:500) = T0 + 2;
y_zad_v = zeros(N*ny,1);

x = zeros(ny,steps);
y = zeros(ny,steps);

offset = 1000;
Fd(1:steps) = FD0;
Fd(600:800) = FD0+2;
dz = [zeros(1,offset+1), Fd(2:steps)-Fd(1:steps-1),0];

Fc(1:steps) = 0;
Fh(1:steps) = 0;
du = zeros(1,2*steps+offset);

for k = 1:(steps)
    y_zad_v(1:2:N*ny,1) = y_zad(1, k+1:k+N)'-h0;
    y_zad_v(2:2:N*ny,1) = y_zad(2, k+1:k+N)'-T0;

    dUp = flip(du((offset+(k-(D-1))*nu-1):(offset+(k-1)*nu)))';
    dZp = flip(dz((offset+(k-Dz+1)):(offset+k)))';

    du_n = K1*(repmat(y_zad(:,k)-[h0;T0],N,1)-repmat(x(:,k),N,1)-Mp*dUp-MPZ*dZp);
    % du_n = K1*(y_zad_v-repmat(x(:,k),N,1)-Mp*dUp);
    du_n(1) = min(du_max(1),max(-du_max(1),du_n(1))); % Ograniczenie ΔFH
    du_n(2) = min(du_max(2),max(-du_max(2),du_n(2))); % Ograniczenie ΔFC

    % du(offset+k*nu-1:offset+k*nu) = flip(du_n)';

    Fh(k) = Fh(max(1,k-1)) + du_n(1);           % Dodanie przyrostu do sterowania
    Fc(k) = Fc(max(1,k-1)) + du_n(2);           % Dodanie przyrostu do sterowania    
    Fh(k) = min(u_max(1),max(u_min(1),Fh(k)));  % Ograniczenie FH
    Fc(k) = min(u_max(2),max(u_min(2),Fc(k)));  % Ograniczenie FC

    du(offset+k*nu-1:offset+k*nu) = flip([Fh(k)-Fh(max(1,k-1)),Fc(k)-Fc(max(1,k-1))]);
    
    Fh_now = Fh(max(1,k-round(tauH/Tp)));
    Fc_now = Fc(max(1,k-round(tauC/Tp)));
    Fd_now = Fd(k) - FD0;

    % Implementacja równań stanu
    x(:,k+1) = A_d*[x(1,k); x(2,k)]+B_d*[Fh_now;Fc_now;Fd_now];
    y(:,k)   = C_d*[x(1,k); x(2,k)];
    
    % Przesunięcie względem warunków początkowych w punkcie równowagi
    wyh(k)=y(1,k) + h0;
    wyT(k)=y(2,k) + T0;
    czas(k)=(k-1)*Tp;
end
    
    % Wykres wysokości słupa cieczy
    fig_hT = figure();
    %subplot(1,2,1);
    subtightplot(2,2,1,gap,margin_h,margin_v);
    grid on;
    hold on;
    stairs(czas,wyh,'r','LineWidth',LineWidth);
    stairs(czas,y_zad(1,1:steps),'k--');
    title('Wysokość słupa cieczy h');
    %xlabel('t [s]');
    ylabel('h [cm]');
    legend('h(t)','h_{zad}(t)')

    % Wykres temperatury cieczy
    % subplot(1,2,2);
    subtightplot(2,2,2,gap,margin_h,margin_v)
    grid on;
    hold on;
    stairs(czas,wyT,'r','LineWidth',LineWidth);
    stairs(czas,y_zad(2,1:steps),'k--');
    title('Temperatura cieczy T');
    %xlabel('t [s]');
    ylabel('T [°C]');
    legend('T(t)','T_{zad}(t)')

    % Wykres Sterowania FH
    % fig_FcFh = figure();
    subtightplot(2,2,3,gap,margin_h,margin_v)
    grid on;
    hold on;
    stairs(czas,Fh+FH0,'Color',[1, 0.467, 0.467],'LineWidth',LineWidth);
    title('Sterowanie F_H');
    xlabel('t [s]');
    ylabel('F_H [cm^3/s]');
    legend('F_H(t)')

    % Wykres Sterowania FC
    subtightplot(2,2,4,gap,margin_h,margin_v)
    grid on;
    hold on;
    stairs(czas,Fc+FC0,'Color',[1, 0.467, 0.467],'LineWidth',LineWidth);
    title('Sterowanie F_C');
    xlabel('t [s]');
    ylabel('F_C [cm^3/s]');
    legend('F_C(t)')

    set(gcf,'units','points','position',[100,100,900,260])

%---------------------------------------------------- PORÓWNANIE

Tp = 10;    % Okres próbkowania
N = 60;     % Horyzont predykcji
Nu = 20;    % Horyzont sterowania
ny = 2;     % ilość wyjść
nu = 2;     % ilość wejść
D = 60;     % Horyzont dynamiki
Dz = 40;    % Horyzont dynamiki zakłócenia
psi = diag(repmat([1 1],1,N));
lambda = diag(repmat([1 1],1,Nu));

G_s = [G11 G12; G21 G22];
G_z = c2d(0.8*G_s, Tp);

S = DMCstepmatrices(Tp, D, G_z);
[M,Mp] = DMCmatrices(S,N,Nu);
K = inv(M'*psi*M+lambda)*M'*psi;
K1 = K(1:nu,:);

% Odpowiedź skokowa modelu zakłócenia
G_s_dist = [G13; G23];
G_z_dist = c2d(G_s_dist, Tp);
S_dist = DMCstepmatrices(Tp, Dz, G_z_dist);
MPZ = DMCdistMatrix(S_dist, N, Nu);

% Czas próbkowania i czas symulacji
T_max = 10000;
steps = round(T_max/Tp);

% Dyskretne równania stanu A_d, B_d, C_d, D_d
A_d = expm(A_m*Tp)
B_d = (A_m^(-1))*(A_d-eye(2,2))*B_m
C_d = C_m;
D_d = D_m;

du_max = [4; 4];
u_max = [400; 400];
u_min = [-400; -400];

y_zad_v = zeros(N*ny,1);

x = zeros(ny,steps);
y = zeros(ny,steps);

Fc(1:steps) = 0;
Fh(1:steps) = 0;
du = zeros(1,2*steps+offset);

for k = 1:(steps)
    y_zad_v(1:2:N*ny,1) = y_zad(1, k+1:k+N)'-h0;
    y_zad_v(2:2:N*ny,1) = y_zad(2, k+1:k+N)'-T0;

    dUp = flip(du((offset+(k-(D-1))*nu-1):(offset+(k-1)*nu)))';
    dZp = flip(dz((offset+(k-Dz+1)):(offset+k)))';

    % repmat(y_zad(:,k),N,1)
    % return

    du_n = K1*(repmat(y_zad(:,k)-[h0;T0],N,1)-repmat(x(:,k),N,1)-Mp*dUp-MPZ*dZp);
    % du_n = K1*(y_zad_v-repmat(x(:,k),N,1)-Mp*dUp);
    du_n(1) = min(du_max(1),max(-du_max(1),du_n(1))); % Ograniczenie ΔFH
    du_n(2) = min(du_max(2),max(-du_max(2),du_n(2))); % Ograniczenie ΔFC

    % du(offset+k*nu-1:offset+k*nu) = flip(du_n)';

    Fh(k) = Fh(max(1,k-1)) + du_n(1);           % Dodanie przyrostu do sterowania
    Fc(k) = Fc(max(1,k-1)) + du_n(2);           % Dodanie przyrostu do sterowania    
    Fh(k) = min(u_max(1),max(u_min(1),Fh(k)));  % Ograniczenie FH
    Fc(k) = min(u_max(2),max(u_min(2),Fc(k)));  % Ograniczenie FC

    du(offset+k*nu-1:offset+k*nu) = flip([Fh(k)-Fh(max(1,k-1)),Fc(k)-Fc(max(1,k-1))]);
    
    Fh_now = Fh(max(1,k-round(tauH/Tp)));
    Fc_now = Fc(max(1,k-round(tauC/Tp)));
    Fd_now = Fd(k) - FD0;

    % Implementacja równań stanu
    x(:,k+1) = A_d*[x(1,k); x(2,k)]+B_d*[Fh_now;Fc_now;Fd_now];
    y(:,k)   = C_d*[x(1,k); x(2,k)];
    
    % Przesunięcie względem warunków początkowych w punkcie równowagi
    wyh(k)=y(1,k) + h0;
    wyT(k)=y(2,k) + T0;
    czas(k)=(k-1)*Tp;
end

figure(fig_hT);
subtightplot(2,2,1,gap,margin_h,margin_v)
stairs(czas,wyh,'b','LineWidth',LineWidth);
legend('DMC 1: h(t)','h_{zad}(t)', ...
    ' DMC 2: h(t), 20% błąd')
set(gca, 'Children', flipud(get(gca, 'Children')) )
subtightplot(2,2,2,gap,margin_h,margin_v)
stairs(czas,wyT,'b','LineWidth',LineWidth);
legend('DMC 1: T(t)','T_{zad}(t)', ...
    ' DMC 2: T(t), 20% błąd')
set(gca, 'Children', flipud(get(gca, 'Children')) )

% Wykres Sterowania FH
% figure(fig_FcFh);
subtightplot(2,2,3,gap,margin_h,margin_v)
stairs(czas,Fh+FH0,'Color',[0.424, 0.424, 0.961],'LineWidth',LineWidth);
legend('DMC 1: F_h(t)',' DMC 2: F_h(t), 20% błąd')
set(gca, 'Children', flipud(get(gca, 'Children')) )

    % Wykres Sterowania FC
subtightplot(2,2,4,gap,margin_h,margin_v)
stairs(czas,Fc+FC0,'Color',[0.424, 0.424, 0.961],'LineWidth',LineWidth);
legend('DMC 1: F_c(t)',' DMC 2: F_c(t), 20% błąd')
set(gca, 'Children', flipud(get(gca, 'Children')) );





function S=DMCstepmatrices(Tp,D,G_z)
    ny = length(G_z(:,1));
    nu = length(G_z(1,:));
    Y = step(G_z, D*Tp);
    S = zeros(ny,nu,D);
    for i=1:ny
        for j=1:nu
        S(i,j,:)=Y(2:(D+1),i,j);
        end
    end
end

function [M,MP]=DMCmatrices(S,N,Nu)
    ny = length(S(:,1,1));
    nu = length(S(1,:,1));
    D = length(S(1,1,:));
    M = zeros(N*ny,Nu*nu);
    MP = zeros(N*ny,D-1);
    for i=0:(N-1)
        for j=0:(Nu-1)
            t=i+1-j;
            if t > 0
                M(i*ny+1:(i+1)*ny,j*nu+1:(j+1)*nu) = S(:,:,t);
            end
        end
    end
    for i=0:(N-1)
        for j=0:(D-2)
            MP(i*ny+1:(i+1)*ny,j*nu+1:(j+1)*nu) = S(:,:,min(j+i+2,D))-S(:,:,j+1);
        end
    end
end

function MPZ=DMCdistMatrix(S,N,Nu)
    ny = length(S(:,1,1));
    nu = length(S(1,:,1));
    D = length(S(1,1,:));
    MPZ = zeros(N*ny,D-1);
    for i=0:(N-1)
        for j=0:(D-2)
            MPZ(i*ny+1:(i+1)*ny,j*nu+1:(j+1)*nu) = S(:,:,min(j+i+2,D))-S(:,:,j+1);
        end
    end
    MPZ = [zeros(N*ny,nu) MPZ];
    for i=0:(N-1)
        MPZ(i*ny+1:(i+1)*ny,1) = S(:,:,min(i+1,D));
    end
end

function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
% subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input arguments (defaults exist):
%   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis. 
%   marg_h  margins in height in normalized units (0...1)
%            or [lower uppper] for different lower and upper margins 
%   marg_w  margins in width in normalized units (0...1)
%            or [left right] for different left and right margins 
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gap parameter must
%       be defined. For default gap value use empty element- [].      
%
% Usage example: h=subtightplot((2,3,1:2,[0.5,0.2])
if (nargin<4) || isempty(gap),    gap=0.01;  end
if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
if (nargin<5) || isempty(marg_w),  marg_w=marg_h;  end
if isscalar(gap),   gap(2)=gap;  end
if isscalar(marg_h),  marg_h(2)=marg_h;  end
if isscalar(marg_w),  marg_w(2)=marg_w;  end
gap_vert   = gap(1);
gap_horz   = gap(2);
marg_lower = marg_h(1);
marg_upper = marg_h(2);
marg_left  = marg_w(1);
marg_right = marg_w(2);
%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  
% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   
% single subplot dimensions:
%height=(1-(m+1)*gap_vert)/m;
%axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
%width =(1-(n+1)*gap_horz)/n;
%axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;
% merged subplot dimensions:
merged_height=subplot_rows*( height+gap_vert )- gap_vert;
merged_width= subplot_cols*( width +gap_horz )- gap_horz;
% merged subplot position:
merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
pos_vec=[merged_left merged_bottom merged_width merged_height];
% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h=subplot('Position',pos_vec,varargin{:});
if (nargout < 1),  clear h;  end
end
