clear
close all
clc

% User Input
i_x_t = 4;
j_y_t = 24;

H = 1;
k1 = 1;    % for L1
k3 = 9;   % for L2
k4 = 1;    % for H1
k2 = 1;    % for H2

% Inlet width (H1) = H
H1 = k4*H;
% Outlet width (H2) = 2H
H2 = k2*H + k4*H;
% Inlet horz step (L1) = H
L1 = k1*H;
% Total domain length (L2) = 16H
L2 = k1*H + k3*H;

A = 1;
St = 0.5;

rho = 1;
Re = 50;
viscosity = 1/Re;

% CAUTION
% (k1/(k1+k3))*(imax-2) must be integer
% (k2/(k2+k4))*(jmax-2) must be integer

imax = 242;
jmax = 32;

% Tolerance values
epsilon_steady_state = 10^-3;
epsilon_mass_cons = 10^-5;

% Grid size evaluation
delta_x = L2/(imax-2);
delta_y = H2/(jmax-2);

% [Xu,Yu] denote the x and y locations of the points where u velocity is defined
x_locations_u = linspace(0,L2,imax-1);
y_locations_u = zeros(jmax,1);
y_locations_u_counter = 0;
for j=1:jmax
    y_locations_u(j) = y_locations_u_counter;
    if((j==1)||(j==jmax-1))
        y_locations_u_counter = y_locations_u_counter + delta_y/2;
    else
        y_locations_u_counter = y_locations_u_counter + delta_y;
    end    
end

% [Xv,Yv] denote the x and y locations of the points where v velocity is defined
x_locations_v = zeros(imax,1);
y_locations_v = linspace(0,H2,jmax-1);
x_locations_v_counter = 0;
for i=1:imax
    x_locations_v(i) = x_locations_v_counter;
    if((i==1)||(i==imax-1))
        x_locations_v_counter = x_locations_v_counter + delta_x/2;
    else
        x_locations_v_counter = x_locations_v_counter + delta_x;
    end 
end

% [Xp,Yp] denote the x and y locations of the points where pressure is defined
x_locations_p = zeros(imax,1);
y_locations_p = zeros(jmax,1);
x_locations_p_counter = 0;
y_locations_p_counter = 0;
for i=1:imax
    x_locations_p(i) = x_locations_p_counter;
    if((i==1)||(i==imax-1))
        x_locations_p_counter = x_locations_p_counter + delta_x/2;
    else
        x_locations_p_counter = x_locations_p_counter + delta_x;
    end    
end    

for j=1:jmax
    y_locations_p(j) = y_locations_p_counter;
    if((j==1)||(j==jmax-1))
        y_locations_p_counter = y_locations_p_counter + delta_y/2;
    else
        y_locations_p_counter = y_locations_p_counter + delta_y;
    end    
end    

% Creating meshgrid
[Xu,Yu] = meshgrid(x_locations_u,y_locations_u);
[Xv,Yv] = meshgrid(x_locations_v,y_locations_v);
[Xp,Yp] = meshgrid(x_locations_p,y_locations_p);
% BB1
for i=1:k1*(imax-2)/(k1+k3)+1
    Yp(k2*(jmax-2)/(k2+k4)+1,i) = Yp(k2*(jmax-2)/(k2+k4)+1,i) + delta_y/2;
    Yu(k2*(jmax-2)/(k2+k4)+1,i) = Yu(k2*(jmax-2)/(k2+k4)+1,i) + delta_y/2;
end

% LB2
for j=1:k2*(jmax-2)/(k2+k4)+1
    Xp(j,k1*(imax-2)/(k1+k3)+1) = Xp(j,k1*(imax-2)/(k1+k3)+1) + delta_x/2;
    Xv(j,k1*(imax-2)/(k1+k3)+1) = Xv(j,k1*(imax-2)/(k1+k3)+1) + delta_x/2;
end

% Domain data structure creation and Initialization
U = zeros(jmax, imax-1);
V = zeros(jmax-1, imax);
P = zeros(jmax, imax);

U_3d = zeros(jmax, imax-1,455);
V_3d = zeros(jmax-1, imax,455);
P_3d = zeros(jmax, imax,455);


% Initial conditions
U_initial = 0;
V_initial = 0;
P_initial = 0;
P_corr_initial = 0;

% Initial Condition (IC) allotment
U(2:jmax-1, 2:imax-2) = U_initial;
V(2:jmax-2, 2:imax-1) = V_initial;
P(2:jmax-1, 2:imax-1) = P_initial;

% Boundary Parameters
U_top = 0;
V_top = 0;
U_bottom = 0;
V_bottom = 0;
% Define U_left in the outer while loop
%U_left = 0; 
V_left = 0;
% Update in end (Non dirichlet BC)
%U_right = 0;
%V_right = 0;

% Boundary Condition (BC) assignment
% Left
% LB1 - Make this time varying
% U(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = 6*(Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1)-1).*(2-Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1));
V(k2*(jmax-2)/(k2+k4)+2:jmax-2,1) = V_left;
% LB2
U(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = 0;
V(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = 0;
% Bottom
% BB1
U(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = U_bottom;
V(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = V_bottom;
% BB2
U(1,k1*(imax-2)/(k1+k3)+2:imax-2) = U_bottom;
V(1,k1*(imax-2)/(k1+k3)+2:imax-1) = V_bottom;
% Top
U(jmax,:) = U_top;
V(jmax-1,:) = V_top;

% Assigning U_star and V_star values from U, V with IC and BC already implemented
U_star = U; V_star = V;

% Time step calculation
fraction = 0.5;

unsteadiness = 1;

time = 0;
num_time_steps = 0;
num_time_check = 1736;
t_low_lim = num_time_check - 455;

% Outer while loop for steady state convergence
while(num_time_steps <= num_time_check)
    U(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = 6*(Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1)-1).*(2-Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1))*(1 + A*sin(2*pi*St*time));
    U_star(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = 6*(Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1)-1).*(2-Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1))*(1 + A*sin(2*pi*St*time));

    % Time step calculation
    value_1 = max(max(abs(U(:,:))));
    value_2 = max(max(abs(V(:,:))));
    advection_delta_t = 1/(value_1/delta_x + value_2/delta_y);
    diffusion_delta_t = 0.5/(viscosity*(1/delta_x^2 + 1/delta_y^2));
    delta_t = fraction*min(advection_delta_t,diffusion_delta_t);

    U_old = U; V_old = V; P_old = P;                % Copying previous timestep values
    P_corr = zeros(jmax, imax);                     % Defining pressure correction data structure
    P_corr(2:jmax-1, 2:imax-1) = P_corr_initial;    % Initializing Pressure correction to 0 at the beginning of every timestep 
    
    % U velocity prediction
    m_ux = zeros(jmax-1,imax-1);    % mass flux for u in x-dir
    a_ux = zeros(jmax-1,imax-1);    % advection flux for u in x-dir
    d_ux = zeros(jmax-1,imax-1);    % diffusion flux for u in x-dir

    for j=2:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2-1:imax-2
            m_ux(j,i) = rho*(U_old(j,i+1) + U_old(j,i))/2;
            m_ux_plus = max(m_ux(j,i),0); m_ux_minus = min(m_ux(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_ux(j,i) = m_ux_plus*U_old(j,i) + m_ux_minus*U_old(j,i+1);
            d_ux(j,i) = viscosity*(U_old(j,i+1) - U_old(j,i))/delta_x;
        end
    end
    for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
        for i=1:imax-2
            m_ux(j,i) = rho*(U_old(j,i+1) + U_old(j,i))/2;
            m_ux_plus = max(m_ux(j,i),0); m_ux_minus = min(m_ux(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_ux(j,i) = m_ux_plus*U_old(j,i) + m_ux_minus*U_old(j,i+1);
            d_ux(j,i) = viscosity*(U_old(j,i+1) - U_old(j,i))/delta_x;
        end
    end

    m_uy = zeros(jmax-1,imax-1);    % mass flux for u in y-dir
    a_uy = zeros(jmax-1,imax-1);    % advection flux for u in y-dir
    d_uy = zeros(jmax-1,imax-1);    % diffusion flux for u in y-dir

    for j=1:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2:imax-2
            m_uy(j,i) = rho*(V_old(j,i+1) + V_old(j,i))/2;
            m_uy_plus = max(m_uy(j,i),0); m_uy_minus = min(m_uy(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_uy(j,i) = m_uy_plus*U_old(j,i) + m_uy_minus*U_old(j+1,i);
            if(j==1)
                d_uy(j,i) = viscosity*(U_old(j+1,i) - U_old(j,i))/(delta_y/2);
            else
                d_uy(j,i) = viscosity*(U_old(j+1,i) - U_old(j,i))/delta_y;
            end    
        end
    end

    for j=k2*(jmax-2)/(k2+k4)+2-1:jmax-1
        for i=2:imax-2
            m_uy(j,i) = rho*(V_old(j,i+1) + V_old(j,i))/2;
            m_uy_plus = max(m_uy(j,i),0); m_uy_minus = min(m_uy(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_uy(j,i) = m_uy_plus*U_old(j,i) + m_uy_minus*U_old(j+1,i);
            if(((j==k2*(jmax-2)/(k2+k4)+2-1)&&(i>=2)&&(i<=k1*(imax-2)/(k1+k3)+1))||(j==jmax-1))
                d_uy(j,i) = viscosity*(U_old(j+1,i) - U_old(j,i))/(delta_y/2);
            else
                d_uy(j,i) = viscosity*(U_old(j+1,i) - U_old(j,i))/delta_y;
            end    
        end
    end

    for j=2:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2:imax-2
            A_u = (a_ux(j,i) - a_ux(j,i-1))*delta_y + (a_uy(j,i) - a_uy(j-1,i))*delta_x;
            D_u = (d_ux(j,i) - d_ux(j,i-1))*delta_y + (d_uy(j,i) - d_uy(j-1,i))*delta_x;
            S_u = (P_old(j,i) - P_old(j,i+1))*delta_y;
            U_star(j,i) = U_old(j,i) + (delta_t/(rho*delta_x*delta_y))*(D_u - A_u + S_u);    % Predicted u
        end
    end

    for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
        for i=2:imax-2
            A_u = (a_ux(j,i) - a_ux(j,i-1))*delta_y + (a_uy(j,i) - a_uy(j-1,i))*delta_x;
            D_u = (d_ux(j,i) - d_ux(j,i-1))*delta_y + (d_uy(j,i) - d_uy(j-1,i))*delta_x;
            S_u = (P_old(j,i) - P_old(j,i+1))*delta_y;
            U_star(j,i) = U_old(j,i) + (delta_t/(rho*delta_x*delta_y))*(D_u - A_u + S_u);    % Predicted u
        end
    end

    

    % V velocity prediction
    m_vx = zeros(jmax-1,imax-1);    % mass flux for v in x-dir
    a_vx = zeros(jmax-1,imax-1);    % advection flux for v in x-dir
    d_vx = zeros(jmax-1,imax-1);    % diffusion flux for v in x-dir

    for j=2:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2-1:imax-1
            m_vx(j,i) = rho*(U_old(j+1,i) + U_old(j,i))/2;
            m_vx_plus = max(m_vx(j,i),0); m_vx_minus = min(m_vx(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_vx(j,i) = m_vx_plus*V_old(j,i) + m_vx_minus*V_old(j,i+1);
            if(((i==k1*(imax-2)/(k1+k3)+2-1)||(i==imax-1))&&((j>=2)&&(j<=k2*(jmax-2)/(k2+k4)+1)))
                d_vx(j,i) = viscosity*(V_old(j,i+1) - V_old(j,i))/(delta_x/2);
            else
                d_vx(j,i) = viscosity*(V_old(j,i+1) - V_old(j,i))/delta_x;
            end
        end
    end
    for j=k2*(jmax-2)/(k2+k4)+2:jmax-2
        for i=1:imax-1
            m_vx(j,i) = rho*(U_old(j+1,i) + U_old(j,i))/2;
            m_vx_plus = max(m_vx(j,i),0); m_vx_minus = min(m_vx(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_vx(j,i) = m_vx_plus*V_old(j,i) + m_vx_minus*V_old(j,i+1);
            if(((i==1)||(i==imax-1))&&((j>=k2*(jmax-2)/(k2+k4)+2)&&(j<=jmax-2)))
                d_vx(j,i) = viscosity*(V_old(j,i+1) - V_old(j,i))/(delta_x/2);
            else
                d_vx(j,i) = viscosity*(V_old(j,i+1) - V_old(j,i))/delta_x;
            end
        end
    end

    m_vy = zeros(jmax-1,imax-1);    % mass flux for v in y-dir
    a_vy = zeros(jmax-1,imax-1);    % advection flux for v in y-dir
    d_vy = zeros(jmax-1,imax-1);    % diffusion flux for v in y-dir

    for j=1:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2:imax-1
            m_vy(j,i) = rho*(V_old(j+1,i) + V_old(j,i))/2;
            m_vy_plus = max(m_vy(j,i),0); m_vy_minus = min(m_vy(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_vy(j,i) = m_vy_plus*V_old(j,i) + m_vy_minus*V_old(j+1,i);
            d_vy(j,i) = viscosity*(V_old(j+1,i) - V_old(j,i))/delta_y;
        end
    end
    for j=k2*(jmax-2)/(k2+k4)+2-1:jmax-2
        for i=2:imax-1
            m_vy(j,i) = rho*(V_old(j+1,i) + V_old(j,i))/2;
            m_vy_plus = max(m_vy(j,i),0); m_vy_minus = min(m_vy(j,i),0);    % splitting mass flux into mass_flux_plus and mass_flux_minus
            a_vy(j,i) = m_vy_plus*V_old(j,i) + m_vy_minus*V_old(j+1,i);
            d_vy(j,i) = viscosity*(V_old(j+1,i) - V_old(j,i))/delta_y;
        end
    end

    
    for j=2:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2:imax-1
            A_v = (a_vx(j,i) - a_vx(j,i-1))*delta_y + (a_vy(j,i) - a_vy(j-1,i))*delta_x;
            D_v = (d_vx(j,i) - d_vx(j,i-1))*delta_y + (d_vy(j,i) - d_vy(j-1,i))*delta_x;
            S_v = (P_old(j,i) - P_old(j+1,i))*delta_x;
            V_star(j,i) = V_old(j,i) + (delta_t/(rho*delta_x*delta_y))*(D_v - A_v + S_v);    % Predicted v
        end
    end
    for j=k2*(jmax-2)/(k2+k4)+2:jmax-2
        for i=2:imax-1
            A_v = (a_vx(j,i) - a_vx(j,i-1))*delta_y + (a_vy(j,i) - a_vy(j-1,i))*delta_x;
            D_v = (d_vx(j,i) - d_vx(j,i-1))*delta_y + (d_vy(j,i) - d_vy(j-1,i))*delta_x;
            S_v = (P_old(j,i) - P_old(j+1,i))*delta_x;
            V_star(j,i) = V_old(j,i) + (delta_t/(rho*delta_x*delta_y))*(D_v - A_v + S_v);    % Predicted v
        end
    end

    % Non dirichlet Right BC
    U_star(2:jmax-1,imax-1) = U_star(2:jmax-1,imax-2);
    V_star(2:jmax-2,imax) = V_star(2:jmax-2,imax-1);

    m_x_star(:,:) = rho*U_star(:,:);    % Predicted mass flux in x-dir
    m_y_star(:,:) = rho*V_star(:,:);    % Predicted mass flux in x-dir
    
    % Smp_star prediction
    Smp_star = zeros(jmax,imax);
    for j=2:k2*(jmax-2)/(k2+k4)+1
        for i=k1*(imax-2)/(k1+k3)+2:imax-1
            Smp_star(j,i) = (m_x_star(j,i) - m_x_star(j,i-1))*delta_y + (m_y_star(j,i) - m_y_star(j-1,i))*delta_x; % Predicted mass source
        end
    end
    for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
        for i=2:imax-1
            Smp_star(j,i) = (m_x_star(j,i) - m_x_star(j,i-1))*delta_y + (m_y_star(j,i) - m_y_star(j-1,i))*delta_x; % Predicted mass source
        end
    end

    aP = delta_t*(2*delta_y/delta_x + 2*delta_x/delta_y);

    % Inner while loop is for iteratively first correcting pressure
    % correction and then updating predicted mass flux value
    while (max(max(abs(Smp_star))) > epsilon_mass_cons)
        for j=2:k2*(jmax-2)/(k2+k4)+1
            for i=k1*(imax-2)/(k1+k3)+2:imax-1
                % Mass correction expression from existing pressure correction value
                Smp_corr = - delta_t*((P_corr(j,i+1)-P_corr(j,i))/delta_x - (P_corr(j,i)-P_corr(j,i-1))/delta_x)*delta_y - delta_t*((P_corr(j+1,i)-P_corr(j,i))/delta_y - (P_corr(j,i)-P_corr(j-1,i))/delta_y)*delta_x; 
                % Pressure correction expression from updated mass correction value
                P_corr(j,i) = P_corr(j,i) - (Smp_star(j,i) + Smp_corr)/aP;
            end
        end
        for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
            for i=2:imax-1
                % Mass correction expression from existing pressure correction value
                Smp_corr = - delta_t*((P_corr(j,i+1)-P_corr(j,i))/delta_x - (P_corr(j,i)-P_corr(j,i-1))/delta_x)*delta_y - delta_t*((P_corr(j+1,i)-P_corr(j,i))/delta_y - (P_corr(j,i)-P_corr(j-1,i))/delta_y)*delta_x; 
                % Pressure correction expression from updated mass correction value
                P_corr(j,i) = P_corr(j,i) - (Smp_star(j,i) + Smp_corr)/aP;
            end
        end

        % Pressure correction BC implemented 
        % LB1
        P_corr(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = P_corr(k2*(jmax-2)/(k2+k4)+2:jmax-1,2);
        % BB1
        P_corr(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = P_corr(k2*(jmax-2)/(k2+k4)+2,2:k1*(imax-2)/(k1+k3)+1);
        % LB2
        P_corr(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = P_corr(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+2);
        % BB2
        P_corr(1,k1*(imax-2)/(k1+k3)+2:imax-1) = P_corr(2,k1*(imax-2)/(k1+k3)+2:imax-1);
        % Right
        P_corr(:,imax) = 0;
        % Top
        P_corr(jmax,:) = P_corr(jmax-1,:);

        % Mass flux correction
        % x-dir
        for j=2:k2*(jmax-2)/(k2+k4)+1
            for i=k1*(imax-2)/(k1+k3)+2-1:imax-2+1
                m_x_star_corr = - delta_t*(P_corr(j,i+1)-P_corr(j,i))/delta_x;
                m_x_star(j,i) = m_x_star(j,i) + m_x_star_corr;
                U_star(j,i) = m_x_star(j,i)/rho;
            end
        end
        for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
            for i=1:imax-1
                m_x_star_corr = - delta_t*(P_corr(j,i+1)-P_corr(j,i))/delta_x;
                m_x_star(j,i) = m_x_star(j,i) + m_x_star_corr;
                U_star(j,i) = m_x_star(j,i)/rho;
            end
        end
        % y-dir
        for j=1:k2*(jmax-2)/(k2+k4)+2
            for i=k1*(imax-2)/(k1+k3)+2:imax-1
                m_y_star_corr = - delta_t*(P_corr(j+1,i)-P_corr(j,i))/delta_y;
                m_y_star(j,i) = m_y_star(j,i) + m_y_star_corr;
                V_star(j,i) = m_y_star(j,i)/rho;
            end
        end
        for j=k2*(jmax-2)/(k2+k4)+2-1:jmax-2+1
            for i=2:imax-1
                m_y_star_corr = - delta_t*(P_corr(j+1,i)-P_corr(j,i))/delta_y;
                m_y_star(j,i) = m_y_star(j,i) + m_y_star_corr;
                V_star(j,i) = m_y_star(j,i)/rho;
            end
        end
        % Evaluation of updated mass source value 
        for j=2:k2*(jmax-2)/(k2+k4)+1
            for i=k1*(imax-2)/(k1+k3)+2:imax-1
                Smp_star(j,i) = (m_x_star(j,i) - m_x_star(j,i-1))*delta_y + (m_y_star(j,i) - m_y_star(j-1,i))*delta_x;
                % Updating pressure by adding pressure correction term
                P(j,i) = P(j,i) + P_corr(j,i);
            end
        end
        for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
            for i=2:imax-1
                Smp_star(j,i) = (m_x_star(j,i) - m_x_star(j,i-1))*delta_y + (m_y_star(j,i) - m_y_star(j-1,i))*delta_x;
                % Updating pressure by adding pressure correction term
                P(j,i) = P(j,i) + P_corr(j,i);
            end
        end
        % Applying boundary condition to pressure
        % LB1
        P(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = P(k2*(jmax-2)/(k2+k4)+2:jmax-1,2);
        % BB1
        P(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = P(k2*(jmax-2)/(k2+k4)+2,2:k1*(imax-2)/(k1+k3)+1);
        % LB2
        P(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = P(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+2);
        % BB2
        P(1,k1*(imax-2)/(k1+k3)+2:imax-1) = P(2,k1*(imax-2)/(k1+k3)+2:imax-1);
        % Right
        P(:,imax) = 0;
        % Top
        P(jmax,:) = P(jmax-1,:);
    end
    
    % At end of pressure correction loop, U_star = U^(n+1) and V_star = V^(n+1)
    U = U_star;
    V = V_star;
    % Non dirichlet Right BC
    U(2:jmax-1,imax-1) = U(2:jmax-1,imax-2);
    V(2:jmax-2,imax) = V(2:jmax-2,imax-1);
    
    % Unsteadiness calculated from current and previous U,V values
    % unstd_u_val = max(max(abs(U-U_old)))/delta_t;
    % unstd_v_val = max(max(abs(V-V_old)))/delta_t;
    % unsteadiness = max(unstd_u_val,unstd_v_val);
     if num_time_steps > t_low_lim
        U_3d(:,:,num_time_steps - t_low_lim) = U;
        V_3d(:,:,num_time_steps - t_low_lim) = V;
        P_3d(:,:,num_time_steps - t_low_lim) = P;
     end 

    time  = time + delta_t;
    num_time_steps = num_time_steps + 1; 
   
end    

% Part 1 figures

% Only for velocity vector plot
Uc = zeros(jmax, imax);
Vc = zeros(jmax, imax);

% Interpolate U to cell centers
for j=2:k2*(jmax-2)/(k2+k4)+1
    for i=k1*(imax-2)/(k1+k3)+2:imax-1
        Uc(j,i) = 0.5 * (U(j,i) + U(j,i-1));
    end
end
for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
    for i=2:imax-1
        Uc(j,i) = 0.5 * (U(j,i) + U(j,i-1));
    end
end

% Interpolate V to cell centers
for j=2:k2*(jmax-2)/(k2+k4)+1
    for i=k1*(imax-2)/(k1+k3)+2:imax-1
        Vc(j,i) = 0.5 * (V(j,i) + V(j-1,i));
    end
end
for j=k2*(jmax-2)/(k2+k4)+2:jmax-1
    for i=2:imax-1
        Vc(j,i) = 0.5 * (V(j,i) + V(j-1,i));
    end
end
% Apply BC for all 4 walls

% Boundary Condition (BC) assignment
% LB1 - Make this time varying
Uc(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = 6*(Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1)-1).*(2-Yu(k2*(jmax-2)/(k2+k4)+2:jmax-1,1))*(1 + A*sin(2*pi*St*(time-delta_t)));
Vc(k2*(jmax-2)/(k2+k4)+2:jmax-1,1) = V_left;
% LB2
Uc(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = 0;
Vc(2:k2*(jmax-2)/(k2+k4)+1,k1*(imax-2)/(k1+k3)+1) = 0;
% BB1
Uc(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = U_bottom;
Vc(k2*(jmax-2)/(k2+k4)+1,2:k1*(imax-2)/(k1+k3)+1) = V_bottom;
% BB2
Uc(1,k1*(imax-2)/(k1+k3)+2:imax-1) = U_bottom;
Vc(1,k1*(imax-2)/(k1+k3)+2:imax-1) = V_bottom;
% Top
Uc(jmax,:) = U_top;
Vc(jmax,:) = V_top;

figure(1)
plot([0, L2],[H2, H2],'k','LineWidth',2)     % Top boundary
hold on
plot([L1, L1], [0, H2 - H1], 'k','LineWidth',2)    % LB2
plot([0, L1], [H2 - H1, H2 - H1], 'k','LineWidth',2)       % BB1
% plot([0, 0], [H1, H2], 'b','LineWidth',2)            % LB1
% plot([L2, L2], [0, H2], 'b','LineWidth',2)         % Right boundary
plot([L1, L2], [0, 0], 'k','LineWidth',2)       % BB2
xlim([-2, L2+2])
ylim([-1, H2+1])
h = streamslice(Xp, Yp, Uc, Vc, 2); 
set(h, 'Color', 'blue', 'LineWidth', 1.2);
title(sprintf('Streamlines at Re = %0.2f,\n L1 = %0.2f*H, L2 = %0.2f*H,\n H1 = %0.2f*H, H2 = %0.2f*H, Timestep = %0.2f', Re, k1, (k1+k3), k4, (k2+k4),num_time_check))
hold off

% plotting U velocity contour and V velocity contour


% removing the points from plotting which is not in the domain
U(1:k2*(jmax-2)/(k2+k4),1:k1*(imax-2)/(k1+k3)) = NaN; % taking limits from LB2 and BB1
V(1:k2*(jmax-2)/(k2+k4),1:k1*(imax-2)/(k1+k3)) = NaN;

% plotting U velocity contour
figure(2)
plot([0, L2],[H2, H2],'k','LineWidth',2)     % Top boundary
hold on
plot([L1, L1], [0, H2 - H1], 'k','LineWidth',2)    % LB2
plot([0, L1], [H2 - H1, H2 - H1], 'k','LineWidth',2)       % BB1
% plot([0, 0], [H1, H2], 'b','LineWidth',2)           % LB1
% plot([L2, L2], [0, H2], 'b','LineWidth',2)        % Right boundary
plot([L1, L2], [0, 0], 'k','LineWidth',2)      % BB2
contourf(Xu,Yu,U)
colorbar
xlabel('X')
ylabel('Y')
xlim([-2, L2+2])
ylim([-1, H2+1])
title(sprintf('U Velocity Contour at Re = %0.2f,\n L1 = %0.2f*H, L2 = %0.2f*H,\n H1 = %0.2f*H, H2 = %0.2f*H, Timestep = %0.2f', Re, k1, (k1+k3), k4, (k2+k4),num_time_check))
hold off

figure(3)
plot([0, L2],[H2, H2],'k','LineWidth',2)    % Top boundary
hold on
plot([L1, L1], [0, H2 - H1], 'k','LineWidth',2)     % LB2
plot([0, L1], [H2 - H1, H2 - H1], 'k','LineWidth',2)     % BB1
% plot([0, 0], [H1, H2], 'b','LineWidth',2)     % LB1
% plot([L2, L2], [0, H2], 'b','LineWidth',2)      % Right boundary
plot([L1, L2], [0, 0], 'k','LineWidth',2)    % BB2
contourf(Xv,Yv,V) 
colorbar
xlabel('X')
ylabel('Y')
xlim([-2, L2+2])
ylim([-1, H2+1])
title(sprintf('V Velocity Contour at Re = %0.2f,\n L1 = %0.2f*H, L2 = %0.2f*H,\n H1 = %0.2f*H, H2 = %0.2f*H, Timestep = %0.2f', Re, k1, (k1+k3), k4, (k2+k4),num_time_check))
hold off

% 1. Define the start and end indices for your slice
idx_start = num_time_steps - t_low_lim - 300;
idx_end = num_time_steps - t_low_lim - 50;

% 2. Find max and min values and their LOCAL indices within the slice
[val_max, local_index_max] = max(U_3d(j_y_t, i_x_t, idx_start:idx_end));
[val_min, local_index_min] = min(U_3d(j_y_t, i_x_t, idx_start:idx_end));

% 3. Convert local indices to ABSOLUTE indices corresponding to U_3d's 3rd dimension
abs_index_max = idx_start + local_index_max - 1;
abs_index_min = idx_start + local_index_min - 1;
abs_index_mean = floor((abs_index_max + abs_index_min)/2);

% 4. Convert absolute indices to actual TIME STEPS (x-axis coordinates)
t_max_step = t_low_lim + abs_index_max;
t_min_step = t_low_lim + abs_index_min;
t_mean_step = t_low_lim + abs_index_mean;
val_mean = U_3d(j_y_t,i_x_t,abs_index_mean);
figure(4)
plot(t_low_lim+1:num_time_check-1, squeeze(U_3d(j_y_t,i_x_t,1:454)))
hold on
plot([t_max_step, t_mean_step,t_min_step], [val_max,val_mean,val_min],'o', 'MarkerSize', 8, 'LineWidth', 1.5)

% Add text labels for the coordinates
% The spaces before the opening parenthesis provide a small horizontal offset
text(t_max_step, val_max, sprintf('  (%d, %.4f)', t_max_step, val_max), 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(t_mean_step, val_mean, sprintf('  (%d, %.4f)', t_mean_step, val_mean), 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(t_min_step, val_min, sprintf('  (%d, %.4f)', t_min_step, val_min), 'VerticalAlignment', 'top', 'FontWeight', 'bold');

xlabel('time')
ylabel('U-velocity')
xlim([t_low_lim-5, num_time_check+5])
%ylim([-1, H2+1])
title(sprintf('U Velocity at Re = %0.2f,\n L1 = %0.2f*H, L2 = %0.2f*H,\n H1 = %0.2f*H, H2 = %0.2f*H \n x = %0.2f, y = %0.2f', Re, k1, (k1+k3), k4, (k2+k4),x_locations_u(i_x_t),y_locations_u(j_y_t)))
hold off