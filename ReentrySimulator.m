%% ============================================================
%  ATMOSPHERIC RE-ENTRY TRAJECTORY SIMULATOR
%  Simulates a capsule (Orion-like) re-entering Earth's atmosphere
%  Models: drag, lift, gravity, aerodynamic heating
%  Animates: trajectory arc colored by temperature/heating rate
%  Tools: MATLAB | Author: [Your Name]
%% ============================================================

clc; clear; close all;

%% ============================================================
%  SECTION 1: CONSTANTS & VEHICLE PARAMETERS
%% ============================================================

% Earth parameters
R_earth  = 6378.0e3;     % Earth radius (m)
mu_earth = 3.986e14;     % gravitational parameter (m^3/s^2)
g0       = 9.81;         % sea-level gravity (m/s^2)

% Atmospheric model constants (exponential approximation)
rho_0    = 1.225;        % sea-level density (kg/m^3)
H_scale  = 7500;         % scale height (m)

% Vehicle parameters — Orion-like capsule
m        = 9300;         % mass (kg)
D        = 5.0;          % capsule diameter (m)
A_ref    = pi*(D/2)^2;   % reference area (m^2)
CD       = 1.2;          % drag coefficient (blunt body)
CL       = 0.3;          % lift coefficient (trimmed)
beta     = m / (CD*A_ref); % ballistic coefficient (kg/m^2)

% Nose radius for heating calculation (m)
R_nose   = 1.0;

%% ============================================================
%  SECTION 2: INITIAL CONDITIONS
%  Entry interface typically at ~120 km altitude
%  Entry angle: shallow = more heating time, steep = more peak heating
%% ============================================================

h0        = 120e3;               % entry altitude (m)
gamma0    = deg2rad(-6.5);       % entry flight path angle (deg, negative = diving)
V0        = 10800;               % entry velocity (m/s) — ~Mach 32 at Orion return
r0        = R_earth + h0;        % initial radial distance (m)
theta0    = 0;                   % initial downrange angle (rad)

% State vector: [r, theta, V, gamma]
%   r     = radial distance from Earth center (m)
%   theta = downrange angle (rad)
%   V     = velocity magnitude (m/s)
%   gamma = flight path angle (rad, negative = descending)
s0 = [r0; theta0; V0; gamma0];

%% ============================================================
%  SECTION 3: EQUATIONS OF MOTION
%  2D point-mass re-entry equations in polar coordinates
%  Includes gravity, drag, lift, and rotating Earth approximation
%
%  dr/dt     = V * sin(gamma)
%  dtheta/dt = V * cos(gamma) / r
%  dV/dt     = -D/m - g*sin(gamma)
%  dgamma/dt = L/(m*V) + (V/r - g/V)*cos(gamma)
%%  ============================================================

function dsdt = reentry_eom(t, s, m, A_ref, CD, CL, R_earth, mu_earth, rho_0, H_scale)
    r     = s(1);
    V     = s(3);
    gamma = s(4);

    % Altitude and gravity
    h   = r - R_earth;
    g   = mu_earth / r^2;

    % Atmospheric density (exponential model)
    rho = rho_0 * exp(-h / H_scale);

    % Aerodynamic forces per unit mass
    q       = 0.5 * rho * V^2;     % dynamic pressure (Pa)
    D_force = q * CD * A_ref / m;  % drag deceleration (m/s^2)
    L_force = q * CL * A_ref / m;  % lift acceleration (m/s^2)

    % Equations of motion
    drdt     =  V * sin(gamma);
    dthetadt =  V * cos(gamma) / r;
    dVdt     = -D_force - g * sin(gamma);
    dgammadt =  L_force/(V) + (V/r - g/V) * cos(gamma);

    dsdt = [drdt; dthetadt; dVdt; dgammadt];
end

%% ============================================================
%  SECTION 4: INTEGRATE TRAJECTORY
%  Stop when altitude < 10 km (parachute deployment)
%% ============================================================

ode_fn  = @(t,s) reentry_eom(t, s, m, A_ref, CD, CL, R_earth, mu_earth, rho_0, H_scale);
opts    = odeset('Events', @(t,s) chute_event(t,s,R_earth), ...
                 'RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', 1.0);

fprintf('Integrating re-entry trajectory...\n');
[t_traj, S_traj] = ode45(ode_fn, [0 600], s0, opts);

% Extract states
r_traj     = S_traj(:,1);
theta_traj = S_traj(:,2);
V_traj     = S_traj(:,3);
gamma_traj = S_traj(:,4);
h_traj     = r_traj - R_earth;

% Cartesian coordinates for plotting
x_traj = r_traj .* cos(theta_traj);
z_traj = r_traj .* sin(theta_traj);

%% ============================================================
%  SECTION 5: AERODYNAMIC HEATING
%  Sutton-Graves stagnation point heating rate:
%  q_dot = k * sqrt(rho/R_nose) * V^3
%  k = 1.7415e-4 (SI units, W/m^2)
%  Total heat load = integral of q_dot over time
%% ============================================================

k_sg   = 1.7415e-4;    % Sutton-Graves constant
h_traj_m = h_traj;

rho_traj  = rho_0 * exp(-h_traj_m / H_scale);
q_dot     = k_sg * sqrt(rho_traj / R_nose) .* V_traj.^3;   % W/m^2
q_dot_MW  = q_dot / 1e6;                                     % MW/m^2

% Dynamic pressure along trajectory
q_dyn     = 0.5 * rho_traj .* V_traj.^2 / 1e3;  % kPa

% Deceleration in G
g_load    = abs(diff(V_traj) ./ diff(t_traj)) / g0;

%% ============================================================
%  SECTION 6: KEY METRICS
%% ============================================================

[q_dot_peak, idx_peak] = max(q_dot_MW);
[g_peak, idx_g]        = max(g_load);
V_mach_entry           = V0 / 340;

fprintf('\n=====================================================\n');
fprintf('   ATMOSPHERIC RE-ENTRY SIMULATION RESULTS\n');
fprintf('=====================================================\n');
fprintf('Vehicle:         Orion-like capsule\n');
fprintf('Entry velocity:  %.0f m/s  (Mach %.1f)\n', V0, V_mach_entry);
fprintf('Entry angle:     %.1f deg\n', rad2deg(gamma0));
fprintf('Entry altitude:  %.0f km\n', h0/1e3);
fprintf('Ballistic coeff: %.1f kg/m^2\n', beta);
fprintf('-----------------------------------------------------\n');
fprintf('Peak heating:    %.2f MW/m^2  at %.1f km\n', ...
        q_dot_peak, h_traj(idx_peak)/1e3);
fprintf('Peak G-load:     %.1f G  at %.1f km\n', ...
        g_peak, h_traj(min(idx_g,length(h_traj)))/1e3);
fprintf('Total downrange: %.0f km\n', ...
        R_earth * (theta_traj(end) - theta_traj(1)) / 1e3);
fprintf('Flight time:     %.1f s  (%.1f min)\n', ...
        t_traj(end), t_traj(end)/60);
fprintf('=====================================================\n\n');

%% ============================================================
%  SECTION 7: STATIC PLOTS
%% ============================================================

fig1 = figure('Name','Re-entry Analysis','Color',[0.06 0.06 0.10],...
              'Position',[60 60 1200 800]);

ax_bg = [0.10 0.10 0.16];
style_ax = @(ax) set(ax,'Color',ax_bg,...
    'XColor',[0.7 0.7 0.7],'YColor',[0.7 0.7 0.7],...
    'GridColor',[0.22 0.22 0.30],'GridAlpha',1,...
    'XGrid','on','YGrid','on','FontSize',9,'FontName','Consolas');

% Color trajectory by heating rate
n_seg  = length(t_traj)-1;
cmap_h = hot(256);
q_norm = (q_dot_MW - min(q_dot_MW)) / (max(q_dot_MW) - min(q_dot_MW) + 1e-10);

%% --- Plot 1: Altitude vs velocity (entry corridor) ---
ax1 = subplot(2,3,1);
scatter(V_traj/1e3, h_traj/1e3, 8, q_dot_MW, 'filled');
colormap(ax1, hot); colorbar;
style_ax(ax1);
xlabel('Velocity (km/s)','FontName','Consolas');
ylabel('Altitude (km)','FontName','Consolas');
title('Entry corridor','Color','w','FontSize',11);
cb1 = colorbar(ax1); cb1.Color=[0.7 0.7 0.7];
cb1.Label.String='Heating (MW/m^2)';
cb1.Label.Color=[0.7 0.7 0.7];

%% --- Plot 2: Heating rate vs time ---
ax2 = subplot(2,3,2);
plot(t_traj, q_dot_MW,'Color',[1.0 0.45 0.15],'LineWidth',2); hold on;
scatter(t_traj(idx_peak), q_dot_MW(idx_peak), 80,'y','filled','MarkerEdgeColor','w');
text(t_traj(idx_peak)+3, q_dot_MW(idx_peak), ...
     sprintf('%.2f MW/m^2\n@ %.0f km', q_dot_peak, h_traj(idx_peak)/1e3),...
     'Color','y','FontSize',8,'FontName','Consolas');
style_ax(ax2);
xlabel('Time (s)','FontName','Consolas');
ylabel('Heating rate (MW/m^2)','FontName','Consolas');
title('Stagnation point heating','Color','w','FontSize',11);

%% --- Plot 3: G-load vs time ---
ax3 = subplot(2,3,3);
plot(t_traj(1:end-1), g_load,'Color',[0.27 0.82 1.0],'LineWidth',2);
style_ax(ax3);
xlabel('Time (s)','FontName','Consolas');
ylabel('Deceleration (G)','FontName','Consolas');
title('G-load','Color','w','FontSize',11);

%% --- Plot 4: Altitude vs time ---
ax4 = subplot(2,3,4);
plot(t_traj, h_traj/1e3,'Color',[0.4 1.0 0.6],'LineWidth',2);
style_ax(ax4);
xlabel('Time (s)','FontName','Consolas');
ylabel('Altitude (km)','FontName','Consolas');
title('Altitude vs time','Color','w','FontSize',11);

%% --- Plot 5: Dynamic pressure ---
ax5 = subplot(2,3,5);
plot(t_traj, q_dyn,'Color',[1.0 0.85 0.2],'LineWidth',2);
style_ax(ax5);
xlabel('Time (s)','FontName','Consolas');
ylabel('Dynamic pressure (kPa)','FontName','Consolas');
title('Dynamic pressure','Color','w','FontSize',11);

%% --- Plot 6: Velocity vs time ---
ax6 = subplot(2,3,6);
plot(t_traj, V_traj/1e3,'Color',[1.0 0.45 0.2],'LineWidth',2);
style_ax(ax6);
xlabel('Time (s)','FontName','Consolas');
ylabel('Velocity (km/s)','FontName','Consolas');
title('Velocity vs time','Color','w','FontSize',11);

sgtitle(sprintf('Re-Entry Simulation  |  V_{entry}=%.0f m/s  |  \\gamma=%.1f°  |  Peak heat=%.2f MW/m^2',...
        V0, rad2deg(gamma0), q_dot_peak),...
        'Color','w','FontSize',13,'FontWeight','bold','FontName','Consolas');

%% ============================================================
%  SECTION 8: ANIMATED 3D TRAJECTORY
%  Shows capsule arc around Earth, glowing by heating rate
%  Save as MP4 for GitHub/LinkedIn portfolio
%% ============================================================

fprintf('Launching animation...\n');
fprintf('(Close animation window to finish)\n\n');

fig2 = figure('Name','Re-entry Animation','Color',[0.02 0.02 0.06],...
              'Position',[100 50 1100 820]);
ax_anim = axes('Parent',fig2);
set(ax_anim,'Color',[0.03 0.03 0.08],...
    'XColor',[0.3 0.3 0.4],'YColor',[0.3 0.3 0.4],'ZColor',[0.3 0.3 0.4],...
    'GridColor',[0.1 0.1 0.15],'GridAlpha',1,...
    'XGrid','on','YGrid','on','ZGrid','on');
hold on; axis equal;

% Earth sphere
[xe,ye,ze] = sphere(60);
surf(xe*R_earth/1e3, ye*R_earth/1e3, ze*R_earth/1e3,...
     'FaceColor',[0.15 0.35 0.65],'EdgeColor','none','FaceAlpha',0.85);

% Atmosphere shell (thin blue glow)
r_atm = (R_earth + 120e3)/1e3;
surf(xe*r_atm, ye*r_atm, ze*r_atm,...
     'FaceColor',[0.3 0.5 1.0],'EdgeColor','none','FaceAlpha',0.06);

% Full trajectory preview (faint)
x_km = x_traj/1e3;
z_km = z_traj/1e3;
plot3(x_km, zeros(size(x_km)), z_km,...
      '--','Color',[0.4 0.4 0.5],'LineWidth',0.8);

% Axis labels
xlabel('X (km)','FontName','Consolas','Color',[0.5 0.5 0.6]);
ylabel('Y (km)','FontName','Consolas','Color',[0.5 0.5 0.6]);
zlabel('Z (km)','FontName','Consolas','Color',[0.5 0.5 0.6]);
view(-35, 18);

% Pre-build color lookup for heating
q_min = min(q_dot_MW);
q_max = max(q_dot_MW);
cmap_fire = [linspace(0.1,1.0,256)', linspace(0.0,0.4,256)', linspace(0.0,0.0,256)'];

% Animate
skip     = max(1, floor(length(t_traj)/300));  % target ~300 frames
h_trail  = plot3(NaN, NaN, NaN, '-', 'LineWidth', 3, 'Color',[1 0.3 0.0]);
h_cap    = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 14,...
                 'MarkerFaceColor',[1 0.8 0.2],'MarkerEdgeColor','w','LineWidth',1.5);
h_glow   = plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 26,...
                 'MarkerFaceColor',[1 0.4 0.0],'MarkerEdgeColor','none');

title_h  = title('','Color','w','FontSize',12,'FontName','Consolas');

% Trail storage
trail_x = []; trail_z = []; trail_c = [];

for k = 1:skip:length(t_traj)
    % Current heating color
    q_frac = (q_dot_MW(k) - q_min) / (q_max - q_min + 1e-10);
    c_idx  = max(1, round(q_frac * 255) + 1);
    col    = cmap_fire(c_idx, :);

    % Add to trail
    trail_x(end+1) = x_km(k);
    trail_z(end+1) = z_km(k);
    trail_c(end+1) = q_frac;

    % Draw colored trail segments
    if length(trail_x) > 1
        for seg = max(1,length(trail_x)-8):length(trail_x)-1
            q_s   = trail_c(seg);
            ci    = max(1, round(q_s*255)+1);
            cs    = cmap_fire(ci,:);
            plot3([trail_x(seg) trail_x(seg+1)], [0 0], ...
                  [trail_z(seg) trail_z(seg+1)],...
                  '-','Color',cs,'LineWidth',2.5);
        end
    end

    % Capsule marker — size scales with heating
    glow_size = 20 + q_frac * 40;
    set(h_glow,'XData',x_km(k),'YData',0,'ZData',z_km(k),...
               'MarkerSize',glow_size,'MarkerFaceColor',col);
    set(h_cap,'XData',x_km(k),'YData',0,'ZData',z_km(k));

    % Update title
    if k < length(g_load)
        g_now = g_load(k);
    else
        g_now = g_load(end);
    end
    set(title_h,'String',sprintf(...
    't = %.0fs  |  h = %.1f km  |  V = %.2f km/s  |  Heat = %.2f MW/m^2  |  G = %.1f g',...
    t_traj(k), h_traj(k)/1e3, V_traj(k)/1e3, q_dot_MW(k), g_now));

    drawnow limitrate;
    pause(0.01);
end

fprintf('Animation complete.\n');
fprintf('Tip: To save as MP4, uncomment the VideoWriter block below.\n\n');

%% ============================================================
%  OPTIONAL: Save animation as MP4
%  Uncomment to save — takes a few minutes
%% ============================================================

% vid = VideoWriter('reentry_animation.mp4','MPEG-4');
% vid.FrameRate = 30;
% vid.Quality   = 95;
% open(vid);
% % Re-run animation loop here with writeVideo(vid, getframe(fig2)) inside
% close(vid);
% fprintf('Saved reentry_animation.mp4\n');

%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

function [value, isterminal, direction] = chute_event(t, s, R_earth)
    h          = s(1) - R_earth;
    value      = h - 10e3;    % trigger at 10 km altitude
    isterminal = 1;
    direction  = -1;
end