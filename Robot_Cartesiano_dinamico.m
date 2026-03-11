%Limpieza de pantalla
clear all
close all
clc
tic

disp('====================================================');
disp(' INICIANDO MODELO DINÁMICO (ROBOT CARTESIANO 3P) ');
disp('====================================================');

% 1. Declaración de variables simbólicas
syms l1(t) l2(t) l3(t) t 
syms l1p l2p l3p       % Velocidades lineales
syms l1pp l2pp l3pp    % Aceleraciones lineales
syms m1 m2 m3 lc1 lc2 lc3 g % Masas, centros de gravedad y gravedad

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[1 1 1];

%Creamos el vector de coordenadas (en columna)
Q = [l1; l2; l3];
Qp = [l1p; l2p; l3p];
Qpp = [l1pp; l2pp; l3pp];

%Número de grado de libertad del robot
GDL= size(RP,2);

% --- CINEMÁTICA DIRECTA ---
%Articulación 1 
P(:,:,1)= [0; 0; l1];
R(:,:,1)= [1 0 0; 
           0 0 -1;
           0 -1 0];

%Articulación 2 
P(:,:,2)= [0; 0; l2];
R(:,:,2)= [1  0  0; 
           0  0  1;
           0  -1  0];
    
%Articulación 3 
P(:,:,3)= [0; 0; 0];
R(:,:,3)= [1 0 0;
           0 1 0;     
           0 0 1];

Vector_Zeros= zeros(1, 3);

disp('--- Matrices de Transformación Homogénea (Globales) ---');
for i = 1:GDL
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
    if i == 1
       T(:,:,i)= A(:,:,i);
    else
       T(:,:,i)= simplify(T(:,:,i-1)*A(:,:,i));
    end
    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
    
    disp(['Matriz T', num2str(i), ':']);
    pretty(T(:,:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CÁLCULO DE JACOBIANOS (Analítico)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jv_a = sym(zeros(3, GDL));
Jw_a = sym(zeros(3, GDL));

for k= 1:GDL
    if RP(k)==0 % Junta Rotacional
        if k == 1
            Jv_a(:,k)= cross([0;0;1], PO(:,:,GDL));
            Jw_a(:,k)= [0;0;1];
        else
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        end
     else % Junta Prismática
        if k == 1
            Jv_a(:,k)= [0;0;1];
            Jw_a(:,k)= [0;0;0];
        else
            Jv_a(:,k)= RO(:,3,k-1);
            Jw_a(:,k)= [0;0;0];
        end
     end
end    

Jv_a = simplify (Jv_a);
Jw_a = simplify (Jw_a);

disp('--- Jacobiano Lineal (Jv) ---');
pretty (Jv_a);
disp('--- Jacobiano Angular (Jw) ---');
pretty (Jw_a);

disp('Velocidad Lineal del Efector Final (V):');
V = simplify(Jv_a * diff(Q, t));
pretty(V);

disp('Velocidad Angular del Efector Final (W):');
W = simplify(Jw_a * diff(Q, t));
pretty(W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VECTORES A LOS CENTROS DE MASA (CoM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P01_local = [0; 0; lc1];
P12_local = [0; 0; lc2];
P23_local = [0; 0; lc3];

% Posición global exacta de cada Centro de Masa
Pos_cm1 = PO(:,:,1) + RO(:,:,1) * P01_local;
Pos_cm2 = PO(:,:,2) + RO(:,:,2) * P12_local;
Pos_cm3 = PO(:,:,3) + RO(:,:,3) * P23_local;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGÍA CINÉTICA (K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--- Calculando Energías ---');
% Derivamos las posiciones globales para obtener las velocidades lineales (Vc)
Vc1 = diff(Pos_cm1, t);
Vc2 = diff(Pos_cm2, t);
Vc3 = diff(Pos_cm3, t);

% Sustituimos diff(l(t),t) por nuestras variables l1p, l2p, l3p
derivadas_t = [diff(l1(t),t), diff(l2(t),t), diff(l3(t),t)];
Vc1 = subs(Vc1, derivadas_t, Qp.');
Vc2 = subs(Vc2, derivadas_t, Qp.');
Vc3 = subs(Vc3, derivadas_t, Qp.');

% K = (1/2)*m*v^2 (No hay rotación, por lo que la parte de Inercia es 0)
K1 = simplify(0.5 * m1 * (Vc1.' * Vc1));
K2 = simplify(0.5 * m2 * (Vc2.' * Vc2));
K3 = simplify(0.5 * m3 * (Vc3.' * Vc3));

K_Total = simplify(K1 + K2 + K3);
disp('Energía Cinética Total (K):');
pretty(K_Total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGÍA POTENCIAL (U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U = m*g*h (Usando la coordenada Z global de cada centro de masa)
U1 = m1 * g * Pos_cm1(3);
U2 = m2 * g * Pos_cm2(3);
U3 = m3 * g * Pos_cm3(3);

U_Total = simplify(U1 + U2 + U3);
disp('Energía Potencial Total (U):');
pretty(U_Total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELO DINÁMICO (Lagrangiano y Fuerzas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = simplify(K_Total - U_Total);
disp('====================================================');
disp('Lagrangiano (L = K - U):');
disp('====================================================');
pretty(L);

Fuerzas = sym(zeros(3,1));

for i = 1:GDL
    % Derivada respecto a la velocidad
    dL_dq_punto = diff(L, Qp(i));
    
    % Derivada respecto al tiempo
    dt_dL_dq_punto = diff(dL_dq_punto, t);
    dt_dL_dq_punto = subs(dt_dL_dq_punto, derivadas_t, Qp.');
    dt_dL_dq_punto = subs(dt_dL_dq_punto, [diff(l1p,t), diff(l2p,t), diff(l3p,t)], Qpp.');
    
    % Derivada respecto a la posición (usamos directamente l1, l2, l3)
    dL_dq = diff(L, Q(i));
    
    % Ecuación de Euler-Lagrange
    Fuerzas(i) = simplify(dt_dL_dq_punto - dL_dq);
end

disp('====================================================');
disp('FUERZAS DE LOS ACTUADORES PRISMÁTICOS (F1, F2, F3):');
disp('====================================================');
pretty(Fuerzas);

toc