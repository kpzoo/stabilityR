% Control epidemics via targeted thresholds at local scale
function [Ic, Lc, Rc, state] = epiTracking(Pomega, nday, R0, I0, ref)

% Assumptions and notes
% - track a desired reference incidence signal
% - control law computed dynamically via local errors

%% Servo or tracking control

% Controlled incidence, total infectiousness and R
Ic = zeros(1, nday); Lc = Ic; Rc = Ic; state = Ic;
% Initialise deme (no delay), state is whether controlled
Ic(1) = I0; Rc(1) = R0; state(1) = 0;

% Iteratively generate renewal epidemic in deme 
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness (uses i-1 data)
    Lc(i) = Ic(i-1:-1:1)*Pomegat'; 
    
    % Control based on Lam and R and past controls
    [Rc(i), state(i)] = servo(Ic(i-1), ref(i-1), state(i-1), Rc(i-1));
    
    % Renewal incidence
    Ic(i) = poissrnd(Lc(i)*Rc(i));
end


%% Threshold based control and release
function [Rc, state] = servo(c, ref, state, R)

% Assumptons and notes
% - react on max threshold, release on min one
% - c can be I or L or a some measure of incidence

% Control errors for acting and extrema
err = ref - c; Rmax = 8; Rmin = 0.2;
% Control and relaxation gains
ctrl = 0.8; relax = 1.1;

% Scenarios determining law by state
switch(state)
    case 0
        % No control currently acting
        if err > 0
            % Below target so relax 
            Rc = min(R*relax, Rmax);
        else
            % Apply some control
            Rc = max(R*ctrl, Rmin);
            state = 1;
        end
    case 1
        % Controlled state
        if err > 0
            % Stop controlling 
            Rc = min(R*relax, Rmax);
            state = 0;
        else
            % Continue control
             Rc = max(R*ctrl, Rmin);
        end
end

