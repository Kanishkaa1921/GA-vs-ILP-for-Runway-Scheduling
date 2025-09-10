%% ========================= Runway Scheduling GA vs ILP =========================
clc; clear; close all;

% Aircraft sizes to test
Ns = [20, 30, 40, 50];
numTrials = 3;

resultsGA = cell(length(Ns), numTrials);
resultsILP = cell(length(Ns), numTrials);

for i = 1:length(Ns)
    N = Ns(i);
    fprintf("===== Running experiments for N = %d aircraft =====\n", N);

    for t = 1:numTrials
        fprintf(" Trial %d/%d ...\n", t, numTrials);

        % Generate synthetic scheduling problem
        schedTimes = linspace(0, 10*N, N);  % extended horizon for feasibility
        separation = 2 + rand(1,N);        % min separation
        priority   = randi([1 10],1,N);    % priority weights
        isArrival  = rand(1,N) > 0.5;      % arrivals vs departures

        % --- GA Run ---
        tic;
        [gaSol, gaObj] = runGA(N, schedTimes, separation, priority, isArrival);
        gaTime = toc;

        % --- ILP Run ---
        tic;
        [ilpSol, ilpObj] = runILP(N, schedTimes, separation, priority, isArrival);
        ilpTime = toc;

        % Store results
        resultsGA{i,t}  = struct('sol',gaSol,'obj',gaObj,'time',gaTime);
        resultsILP{i,t} = struct('sol',ilpSol,'obj',ilpObj,'time',ilpTime);
    end
end

%% ========================= Analyze Results =========================
gaTimes  = zeros(length(Ns), numTrials);
ilpTimes = zeros(length(Ns), numTrials);
gaObjVal = zeros(length(Ns), numTrials);
ilpObjVal= zeros(length(Ns), numTrials);

for i = 1:length(Ns)
    for t = 1:numTrials
        % --- GA ---
        if ~isempty(resultsGA{i,t}) && ~isempty(resultsGA{i,t}.obj)
            gaObjVal(i,t) = mean(resultsGA{i,t}.obj(:));
        else
            gaObjVal(i,t) = NaN;
        end
        if ~isempty(resultsGA{i,t}) && ~isempty(resultsGA{i,t}.time)
            gaTimes(i,t) = resultsGA{i,t}.time;
        else
            gaTimes(i,t) = NaN;
        end

        % --- ILP ---
        if ~isempty(resultsILP{i,t}) && ~isempty(resultsILP{i,t}.obj)
            ilpObjVal(i,t) = resultsILP{i,t}.obj;
            ilpTimes(i,t)  = resultsILP{i,t}.time;
        else
            ilpObjVal(i,t) = NaN;
            ilpTimes(i,t)  = NaN;
            fprintf('ILP failed for N=%d, Trial=%d\n', Ns(i), t);
        end
    end
end

% Helper function to compute mean ignoring NaNs
meanIgnoringNaN = @(X) mean(X, 2, 'omitnan');

% --- Runtime Plot ---
figure;
plot(Ns, meanIgnoringNaN(gaTimes), '-o','LineWidth',2); hold on;
plot(Ns, meanIgnoringNaN(ilpTimes), '-s','LineWidth',2);
xlabel('Number of Aircraft'); ylabel('Avg Runtime (s)');
legend('GA','ILP'); grid on;
title('GA vs ILP Runtime Scaling');

% --- Objective Plot ---
figure;
plot(Ns, meanIgnoringNaN(gaObjVal), '-o','LineWidth',2); hold on;
plot(Ns, meanIgnoringNaN(ilpObjVal), '-s','LineWidth',2);
xlabel('Number of Aircraft'); ylabel('Avg Objective Value');
legend('GA','ILP'); grid on;
title('GA vs ILP Objective Comparison');

%% ========================= GA Wrapper =========================
function [xPareto,fPareto] = runGA(numAircraft,schedTimes,separation,priority,isArrival)
    nvars = numAircraft;
    lb = schedTimes - 10;
    ub = schedTimes + 10;

    fitnessFcn = @(x) runwayObjectivesStochastic(x,schedTimes,separation,priority,numAircraft,isArrival);
    nonlcon    = @(x) runwayConstraints(x,separation);

    opts = optimoptions('gamultiobj',...
        'MaxGenerations',50,...
        'PopulationSize',100,...
        'CreationFcn',@gacreationsobol,...
        'CrossoverFcn',@crossoverintermediate,...
        'SelectionFcn',@selectiontournament,...
        'MutationFcn',@mutationpower,...
        'Display','off');

    [xPareto,fPareto] = gamultiobj(fitnessFcn,nvars,[],[],[],[],lb,ub,nonlcon,opts);
end

%% ========================= Rectified ILP Wrapper =========================
function [xOpt,fOpt] = runILP(numAircraft,schedTimes,separation,priority,isArrival)
    try
    %%
    % 
    % 
    % 
        % Decision variables: x1..xN + binary y_ij for each aircraft pair
        N = numAircraft;
        pairIdx = nchoosek(1:N,2);
        numPairs = size(pairIdx,1);

        % Bounds for x
        lb = floor(schedTimes - 10);
        ub = ceil(schedTimes + 10);

        % Bounds for binary variables
        lb = [lb zeros(1,numPairs)];
        ub = [ub ones(1,numPairs)];

        % Objective: weighted deviation only for x
        f = [priority(:); zeros(numPairs,1)];

        % Integer constraints for binary variables only
        intcon = N + (1:numPairs);

        % Big-M method
        M = max(ub) - min(lb) + 10;
        A = [];
        b = [];

        for k = 1:numPairs
            i = pairIdx(k,1);
            j = pairIdx(k,2);

            % Constraint: x_i - x_j + M*y_ij >= separation
            row = zeros(1, N + numPairs);
            row(i) = 1; row(j) = -1; 
            row(N + k) = M;
            A = [A; -row]; % intlinprog uses A*x <= b
            b = [b; -separation(i)];

            % Constraint: x_j - x_i + M*(1 - y_ij) >= separation
            row = zeros(1, N + numPairs);
            row(i) = -1; row(j) = 1;
            row(N + k) = -M;
            A = [A; -row];
            b = [b; -separation(i) + M];
        end

        opts = optimoptions('intlinprog','Display','off');
        xFull = intlinprog(f,intcon,A,b,[],[],lb,ub,opts);

        xOpt = xFull(1:N);
        fOpt = sum(priority .* abs(xOpt' - schedTimes));
    catch
        xOpt = [];
        fOpt = [];
    end
end

%% ========================= GA Objective Function =========================
function f = runwayObjectivesStochastic(x,schedTimes,separation,priority,numAircraft,isArrival)
    noise = 2*randn(1,numAircraft);
    xNoisy = x + noise;

    delays = abs(xNoisy - schedTimes);
    f1 = sum(delays);
    f2 = sum(priority .* delays);

    f = [f1, f2];
end

%% ========================= GA Constraints =========================
function [c,ceq] = runwayConstraints(x,separation)
    N = length(x);
    c = [];
    for i = 1:N-1
        for j = i+1:N
            c(end+1) = max(separation(i), separation(j)) - abs(x(j)-x(i));
        end
    end
    ceq = [];
end


