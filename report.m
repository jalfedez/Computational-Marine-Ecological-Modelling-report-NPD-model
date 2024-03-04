% 1. Parameters struct per day
parameters.n = 200;          % number of cells
parameters.d = 100;          % depth
parameters.Deltaz = parameters.d / parameters.n;  % grid width (volume)
parameters.kw = 0.375;      %1/m
parameters.kc =0.05;        %m2/mm N  absortion    =kp

%all values that have time must be in days
parameters.Io = 350 ;             % light at t=0(W/m2)
parameters.gmax = 1;        % division rate max per day
parameters.u = 1;           % sinking velocity(m/day)
parameters.Av = 5;           % diffusion  m2/d
parameters.Kn = 0.3;         % Half saturation constant of nutrient limited growth for N- and I-species mmol nutrient/m3  mmol N / m3
parameters.e = 0.5;        % natural mortality 1/d
parameters.Nb = 30;               % Nutrient concentration at Z bottom (mmol nutrient/m3)
parameters.y = 0.1;         % grazing   m3/mm*N*d
parameters.tau=0.5;         %remineralization 1/d
parameters.w=5;            %sinking dendrite m/d
parameters.alpha=0.1;       %light sensitivity
    
% 2. Set up grid
parameters.z = parameters.Deltaz / 2 : parameters.Deltaz : (parameters.d - parameters.Deltaz / 2);

% 3. Set initial conditions
% Concentration
P0 = zeros(1, parameters.n);
P0(5) = 1;

% Nutrients
N0 = zeros(1, parameters.n)+ parameters.Nb;

%Detritus   
D0 = zeros(1, parameters.n);
D0(3)=1;


% 4. Solve ODE
[t, Y] = ode45(@(t, Y) derivative(t, Y, parameters), [0 1500], [P0, N0, D0]);



P = Y(:, 1:parameters.n);
N = Y(:, parameters.n+1:2*parameters.n);
D= Y(:, 2*parameters.n+1:3*parameters.n);

% Plot results

%  Show that the solution is converged and independent of the grid size
figure;
subplot(1,3,1);
imagesc(t, parameters.z, P');
xlabel('Time(days)');
ylabel('Depth(m)');
title('Phytoplankton Concentration Profile(mmol N/m^3)');
colorbar;

subplot(1,3,2);
imagesc(t, parameters.z, N');
xlabel('Time(days)');
ylabel('Depth(m)');
title('Nutrients Concentration Profile(mmol N/m^3)');
colorbar;

subplot(1,3,3);
imagesc(t, parameters.z, D');
xlabel('Time(days)');
ylabel('Depth(m)');
title('Detritus Concentration Profile(mmol N/m^3)');
colorbar;

% Plot nutrient profile (N) in the first subplot
figure;
hold on;

% Calculate growth-limiting functions
sigmaN = N(end,:)./(parameters.Kn + N(end,:));
I = calcI(P(end,:), D(end,:), parameters);
sigmaL = parameters.alpha* I ./ sqrt(parameters.gmax^2 + parameters.alpha^2 * I.^2);

% Plot growth-limiting functions
plot(sigmaN, -parameters.z, 'b', 'LineWidth', 1.5);
plot(sigmaL, -parameters.z, 'r', 'LineWidth', 1.5);

xlabel('Growth Limiting Function');
ylabel('Depth (m)');
title('Growth Limiting Functions: Nutrient and Light');
legend('\sigma_N (Nutrient)', '\sigma_L (Light)', 'Location', 'Best');

hold off;

% % Preallocate arrays to store results
% % Define the values of parameters.Io to compare
% Io_values = [50, 100, 200, 300, 400];
% depth = parameters.z;
% figure;
% % Solve ODE for each Io value
% hold on;
% for i = 1:length(Io_values)
%     parameters.Io = Io_values(i);
%     [t, Y] = ode45(@(t, Y) derivative(t, Y, parameters), [0 1500], [P0, N0, D0]);
%     P = Y(:, 1:parameters.n);
%     plot(P(end,:),-depth, 'LineWidth', 1.5);
%     
% end
% hold off;
% 
% xlabel('Phytoplankton Concentration (umol/m^3)');
% ylabel('Depth (m)');
% title('Phytoplankton Concentration Profiles for Different Io Values');
% legend(strcat(cellstr(num2str(Io_values')), " W/m^2"), 'Location', 'Best');


% Derivative function
function dYdt = derivative(t, Y, parameters)
    % Extract parameters from struct
    n = parameters.n;
    Deltaz = parameters.Deltaz;  
    %all values that have time must be in days
    gmax=parameters.gmax;        % division rate max per day
    u=parameters.u;           % sinking velocity
    Av=parameters.Av;           % diffusion  1/d
    kn=parameters.Kn;            % Half saturation constant of nutrient limited growth for N- and I-species mmol nutrient/m3  mmol N / m3
    e=parameters.e;           % natural mortality 1/d
    Nb=parameters.Nb;               % Nutrient concentration at Z bottom (mmol nutrient/m3)
    y=parameters.y;         % grazing   m3/mm*N*d
    tau=parameters.tau;         %remineralization 1/d
    w=parameters.w;   
    alpha=parameters.alpha;

    P = Y(1:parameters.n);
    N = Y(parameters.n+1:2*parameters.n);
    D= Y(2*parameters.n+1:3*parameters.n);


    % initialize fluxes
    Jdp = zeros(1, n+1);
    Jap = zeros(1, n+1);
    % calculate fluxes for P
    for i = 2:n
        Jap(i) = u * P(i-1);                        % Advective flux
        Jdp(i) = -Av * ((P(i) - P(i-1)) / Deltaz); % Diffusive flux
    end

    % calculate overall derivative
    JP = Jap+Jdp;

    % initialize fluxes
    JdN = zeros(1, n+1);

    % calculate fluxes, N doesn't have advective fluxes
    for i = 2:n
        JdN(i) = -Av * ((N(i) - N(i-1)) / Deltaz); % Diffusive flux
    end

    % boundary fluxes
    JdN(1)=0;
    JdN(n+1) = -Av * (Nb - N(n)) / Deltaz;


    %Detritus
    % initialize fluxes
    JdD = zeros(1, n+1);
    JaD = zeros(1, n+1);
    % calculate fluxes, N doesn't have advective fluxes
    for i = 2:n
        JaD(i) = w * D(i-1); 
        JdD(i) = ((D(i) - D(i-1)) / Deltaz); % Diffusive flux
    end
    
    JaD(n+1) = w * D(n);

    JD=JaD+JdD;
   

    % calculate I
    I = calcI(P, D, parameters);

    % calculate dP , dD and dN
    dPdt = zeros(1, n);
    dNdt = zeros(1, n);
    dDdt= zeros(1,n);
    sigmaN=N./(kn+N);
    sigmaL=alpha.*I./(sqrt(gmax^2+alpha^2*I.^2));
    for i = 1:n
        dPdt(i) = gmax*sigmaN(i)*sigmaL(i)*P(i)- e*P(i)-y*P(i)^2 - (JP(i+1) - JP(i)) / Deltaz;
        dNdt(i) = -gmax*sigmaN(i)*sigmaL(i)*P(i)+tau*D(i) - (JdN(i+1) - JdN(i)) / Deltaz;
        dDdt(i)= e*P(i) + y*P(i)^2 - tau*D(i) - JD(i)  -(JD(i+1) - JD(i)) / Deltaz;
    end

    dYdt = [dPdt, dNdt, dDdt]';
end

function I = calcI(P,D, parameters)
    % Extract parameters from struct
    n = parameters.n;
    Deltaz = parameters.Deltaz;
    kw = parameters.kw;
    kc = parameters.kc;
    Io = parameters.Io;
    

    % Calculate cumulative light attenuation
    dI = cumsum((kw + kc * (P+D)) * Deltaz) - (1/2) * Deltaz * (kw + kc * (P+D));

    % Calculate light intensity
    I = Io * exp(-dI);
end
