% Full Propeller Analysis at different RPMs, Inflow Velocities, the best
% one so far as of Feb 17, 2026 (Validated)
% Vishal Gautam
% Master's Student
% Virginia Tech, Blacksburg USA

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ── Mode selector ────────────────────────────────────────────────────── %
% test = 1 : fix RPM, sweep velocity  (get CT/CP/CQ vs J)                 %
% test = 0 : fix velocity, sweep RPM  (get Thrust/Torque vs RPM)          %
test = 1;                                                                 %
                                                                          %  
% ── Propeller geometry ───────────────────────────────────────────────── %
% These hard coded values do not work here because the prop geometry      %      
% changes in excel when you change it, so it is done in excel right now   %
%Dp = 24;          % Diameter in inches                                   %
%B  = 12;          % Number of blades                                     %
                                                                          %
% ── Sweep parameters ─────────────────────────────────────────────────── %
fixed_RPM      = 2800;              % RPM used when test=1                % 
fixed_velocity = 00;                % m/s used when test=0                %
                                                                          %
velocity_range = 0:0.5:20;          % m/s  (test=1)                       %
RPM_range      = 100:100:4000;      % RPM  (test=0)                       %
                                                                          %
% ── Atmosphere ───────────────────────────────────────────────────────── %
h      = 634;                       % Blacksburg's height in meters       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ── Atmosphere cont.───────────────────────────────────────────────────────
Tt     = 25;                                                              
rho_0  = 1.225;                                                           
r      = ((273 * (101325 * (1 - (0.0065 * (h / (273 + Tt))))^5.2561)) / ...
          ((101325 * (273 + Tt))) * rho_0);                               

% ── Read blade section data from Excel ────────────────────────────────────
Filename = 'propdata.xlsx';

Dp = xlsread(Filename, 'J2:J2');
P = xlsread(Filename, 'J3:J3');
B = xlsread(Filename, 'J4:J4');
Section_length = xlsread(Filename, 'J9:J9');

Section = xlsread(Filename, 'A3:A88');
Pitch   = xlsread(Filename, 'C3:C88');
Chord   = xlsread(Filename, 'D3:D88');

Cl   = xlsread(Filename, 2, 'B12:B36');
Cd   = xlsread(Filename, 2, 'C12:C36');
Alfa = xlsread(Filename, 2, 'H12:H36');

Di = Dp * 0.0254; % Convert to metres                                     %

%  CORE BEM FUNCTION
function [T, Q] = bem_solve(RPM_rad, velocity, Section, Pitch, Chord, ...
                             Alfa, Cl, Cd, B, r, Section_length)
    T = 0;
    Q = 0;
    for ii = 1:length(Section)
        a  = 0.1;
        ao = 0.01;
        testvalue = 0;
        check     = 1;
        Pitchn    = Pitch(ii);

        while ~testvalue
            Vth = (1 - ao) * RPM_rad * Section(ii);
            Vax = (1 + a)  * velocity;
            v   = sqrt(Vth^2 + Vax^2);

            phi = atan(Vax / Vth);
            ang = Pitchn - phi;

            Clint = interp1(Alfa, Cl, ang, 'spline', 1.5856);
            Cdint = interp1(Alfa, Cd, ang, 'spline', 0.04603);

            if ang < Alfa(1)
                Clint = Cl(1);
                Cdint = Cd(1);
            end

            dT = ((r * v^2) / 2) * (Clint*cos(phi) - Cdint*sin(phi)) * B * Chord(ii);
            dQ = ((r * v^2) / 2) * (Clint*sin(phi) + Cdint*cos(phi)) * B * Chord(ii) * Section(ii);

            anew  = dT  / (4 * pi * Section(ii)   * r * v^2 * (1 + a));
            aonew = dQ  / (4 * pi * Section(ii)^3 * r * v   * (1 + a) * RPM_rad);

            amiddle  = (anew  + a)  / 2;
            aomiddle = (aonew + ao) / 2;

            if (abs(amiddle - a) < 1e-5) && (abs(aomiddle - ao) < 1e-5)
                testvalue = 1;
            end

            check = check + 1;
            a  = amiddle;
            ao = aomiddle;

            if check > 250
                testvalue = 1;
            end
        end

        T = T + dT * Section_length;
        Q = Q + dQ * Section_length;
    end
end

%  TEST = 1 : Fix RPM, sweep velocity >> CT, CP, CQ, eta vs J
if test == 1

    RPM_rad = fixed_RPM * 2 * pi / 60;
    rps     = fixed_RPM / 60;
    nV      = length(velocity_range);

    CT  = zeros(1, nV);
    CP  = zeros(1, nV);
    CQ  = zeros(1, nV);
    eta = zeros(1, nV);
    J   = zeros(1, nV);
    Thr = zeros(1, nV);
    Tor = zeros(1, nV);

    for i = 1:nV
        vel = velocity_range(i);
        [T, Q] = bem_solve(RPM_rad, vel, Section, Pitch, Chord, ...
                           Alfa, Cl, Cd, B, r, Section_length);

        Thr(i) = T;
        Tor(i) = Q;

        CT(i) = T / (r * rps^2 * Di^4);
        CQ(i) = Q / (r * rps^2 * Di^5);
        CP(i) = 2 * pi * CQ(i);
        J(i)  = vel / (rps * Di);

        if (CT(i) > 0) && (CP(i) > 0)
            eta(i) = J(i) * CT(i) / CP(i);
        else
            eta(i) = 0;
        end
    end

    % ── Figure 1: CT vs J ────────────────────────────────────────────────
    figure;
    plot(J, CT, 'b-', 'linewidth', 2);
    title(sprintf('Thrust Coefficient C_T vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('C_T');
    grid on;

    % ── Figure 2: CP vs J ────────────────────────────────────────────────
    figure;
    plot(J, CP, 'r-', 'linewidth', 2);
    title(sprintf('Power Coefficient C_P vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('C_P');
    grid on;

    % ── Figure 3: CQ vs J ────────────────────────────────────────────────
    figure;
    plot(J, CQ, 'm-', 'linewidth', 2);
    title(sprintf('Torque Coefficient C_Q vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('C_Q');
    grid on;

    % ── Figure 4: CT & CP overlaid vs J ──────────────────────────────────
    figure;
    plot(J, CT, 'b-', 'linewidth', 2); hold on;
    plot(J, CP, 'r--', 'linewidth', 2); hold off;
    title(sprintf('C_T and C_P vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('Coefficient');
    legend('C_T', 'C_P');
    grid on;

    % ── Figure 5: Propeller Efficiency vs J ──────────────────────────────
    figure;
    plot(J, eta, 'g-', 'linewidth', 2);
    title(sprintf('Propeller Efficiency \\eta vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('Efficiency \eta');
    ylim([0 1]);
    grid on;

    % ── Figure 6: CT / CP ratio vs J ─────────────────────────────────────
    figure;
    ratio = zeros(size(J));
    for i = 1:nV
        if CP(i) > 0
            ratio(i) = CT(i) / CP(i);
        end
    end
    plot(J, ratio, 'k-', 'linewidth', 2);
    title(sprintf('C_T / C_P vs J  (RPM = %d)', fixed_RPM));
    xlabel('Advance Ratio J');  ylabel('C_T / C_P');
    grid on;

    % ── Figure 7: Thrust vs Velocity ─────────────────────────────────────
    figure;
    plot(velocity_range, Thr, 'b-', 'linewidth', 2);
    title(sprintf('Thrust vs Velocity  (RPM = %d)', fixed_RPM));
    xlabel('Velocity (m/s)');  ylabel('Thrust (N)');
    grid on;

    % ── Figure 8: Torque vs Velocity ─────────────────────────────────────
    figure;
    plot(velocity_range, Tor, 'r-', 'linewidth', 2);
    title(sprintf('Torque vs Velocity  (RPM = %d)', fixed_RPM));
    xlabel('Velocity (m/s)');  ylabel('Torque (Nm)');
    grid on;

%  TEST = 0 : Fix velocity, sweep RPM >> Thrust, Torque, CT, CP vs RPM
else

    nR  = length(RPM_range);
    CT  = zeros(1, nR);
    CP  = zeros(1, nR);
    CQ  = zeros(1, nR);
    eta = zeros(1, nR);
    J   = zeros(1, nR);
    Thr = zeros(1, nR);
    Tor = zeros(1, nR);

    for i = 1:nR
        RPM_rad = RPM_range(i) * 2 * pi / 60;
        rps     = RPM_range(i) / 60;

        [T, Q] = bem_solve(RPM_rad, fixed_velocity, Section, Pitch, Chord, ...
                           Alfa, Cl, Cd, B, r, Section_length);

        Thr(i) = T;
        Tor(i) = Q;

        CT(i) = T / (r * rps^2 * Di^4);
        CQ(i) = Q / (r * rps^2 * Di^5);
        CP(i) = 2 * pi * CQ(i);
        J(i)  = fixed_velocity / (rps * Di);

        if (CT(i) > 0) && (CP(i) > 0)
            eta(i) = J(i) * CT(i) / CP(i);
        else
            eta(i) = 0;
        end
    end

    % ── Figure 1: Thrust vs RPM ───────────────────────────────────────────
    figure;
    plot(RPM_range, Thr, 'b-', 'linewidth', 2);
    title(sprintf('Thrust vs RPM  (V = %g m/s)', fixed_velocity));
    xlabel('RPM');  ylabel('Thrust (N)');
    grid on;

    % ── Figure 2: Torque vs RPM ───────────────────────────────────────────
    figure;
    plot(RPM_range, Tor, 'r-', 'linewidth', 2);
    title(sprintf('Torque vs RPM  (V = %g m/s)', fixed_velocity));
    xlabel('RPM');  ylabel('Torque (Nm)');
    grid on;

    % ── Figure 3: CT vs J ────────────────────────────────────────────────
    figure;
    plot(J, CT, 'b-', 'linewidth', 2);
    title(sprintf('C_T vs Advance Ratio J  (V = %g m/s)', fixed_velocity));
    xlabel('Advance Ratio J');  ylabel('C_T');
    grid on;

    % ── Figure 4: CP vs J ────────────────────────────────────────────────
    figure;
    plot(J, CP, 'r-', 'linewidth', 2);
    title(sprintf('C_P vs Advance Ratio J  (V = %g m/s)', fixed_velocity));
    xlabel('Advance Ratio J');  ylabel('C_P');
    grid on;

    % ── Figure 5: CQ vs J ────────────────────────────────────────────────
    figure;
    plot(J, CQ, 'm-', 'linewidth', 2);
    title(sprintf('C_Q vs Advance Ratio J  (V = %g m/s)', fixed_velocity));
    xlabel('Advance Ratio J');  ylabel('C_Q');
    grid on;

    % ── Figure 6: CT & CP overlaid vs J ──────────────────────────────────
    figure;
    plot(J, CT, 'b-', 'linewidth', 2); hold on;
    plot(J, CP, 'r--', 'linewidth', 2); hold off;
    title(sprintf('C_T and C_P vs J  (V = %g m/s)', fixed_velocity));
    xlabel('Advance Ratio J');  ylabel('Coefficient');
    legend('C_T', 'C_P');
    grid on;

    % ── Figure 7: Efficiency vs J ─────────────────────────────────────────
    figure;
    plot(J, eta, 'g-', 'linewidth', 2);
    title(sprintf('Propeller Efficiency \\eta vs J  (V = %g m/s)', fixed_velocity));
    xlabel('Advance Ratio J');  ylabel('Efficiency \eta');
    ylim([0 1]);
    grid on;

    % ── Figure 8: CT vs RPM ───────────────────────────────────────────────
    figure;
    plot(RPM_range, CT, 'b-', 'linewidth', 2);
    title(sprintf('C_T vs RPM  (V = %g m/s)', fixed_velocity));
    xlabel('RPM');  ylabel('C_T');
    grid on;

    % ── Figure 9: CP vs RPM ───────────────────────────────────────────────
    figure;
    plot(RPM_range, CP, 'r-', 'linewidth', 2);
    title(sprintf('C_P vs RPM  (V = %g m/s)', fixed_velocity));
    xlabel('RPM');  ylabel('C_P');
    grid on;

end