function spine_sweep()
% SPINE_SWEEP
% Run parameter sweeps on the geometric model.
%
% Experiment 1: vary axonDistance, keep irregular angle noise fixed.
% Experiment 2: vary irregular angle noise, keep axonDistance fixed.

%% ----------------- Common base parameters -----------------
baseParams.dendLength = 100;
baseParams.nSpines = 40;

% Regular spines
baseParams.regSpineLength = 1.5;
baseParams.regSpineAngle = 0;

% Irregular spines
baseParams.irPosJitterFrac = 0.25;
baseParams.irLengthMean = 1.8;
baseParams.irLengthStd = 0.7;

% Synapse rule
baseParams.synapseRadius = 0.4;
baseParams.maxAngleDeg = 45;

% Simulation
baseParams.nRuns = 50;

% Precompute alignment threshold
baseParams.maxAngleRad = deg2rad(baseParams.maxAngleDeg);
baseParams.minCosAlignment = cos(baseParams.maxAngleRad);

%% ----------------- Experiment 1: sweep axonDistance -----------------
axonDistances = [1.5 2.0 2.5 3.0 3.5]; 
nAxon = numel(axonDistances);

meanRegExp1 = zeros(nAxon,1);
sdRegExp1 = zeros(nAxon,1);
meanIrExp1 = zeros(nAxon,1);
sdIrExp1 = zeros(nAxon,1);

for i = 1:nAxon
    params = baseParams;
    params.axonDistance = axonDistances(i);
    params.irAngleStd = 20;  % keep irregular angle noise fixed here

    [meanReg, sdReg, meanIr, sdIr] = run_geometric_model(params);

    meanRegExp1(i) = meanReg;
    sdRegExp1(i) = sdReg;
    meanIrExp1(i) = meanIr;
    sdIrExp1(i) = sdIr;
end

% Plot Experiment 1: focus on irregular synapses vs axon distance
figure;
hold on;
errorbar(axonDistances, meanIrExp1, sdIrExp1, '-o', 'LineWidth', 1.5);
errorbar(axonDistances, meanRegExp1, sdRegExp1, '-o', 'LineWidth', 1.5);
xlabel('Axon distance (μm)');
ylabel('Mean potential synapses (irregular)');
title('Experiment 1: synapses vs axon distance');
legend({'Irregular','Regular'}, 'Location', 'best');
grid on;

%% ----------------- Experiment 2: sweep irregular angle noise --------
angleVals = [0 10 20 30 40 50 60]; 
nAng = numel(angleVals);

meanRegExp2 = zeros(nAng,1);
sdRegExp2 = zeros(nAng,1);
meanIrExp2 = zeros(nAng,1);
sdIrExp2 = zeros(nAng,1);

fixedAxonDistance = 2.0;

for i = 1:nAng
    params = baseParams;
    params.axonDistance = fixedAxonDistance;
    params.irAngleStd = angleVals(i);

    [meanReg, sdReg, meanIr, sdIr] = run_geometric_model(params);

    meanRegExp2(i) = meanReg;
    sdRegExp2(i) = sdReg;
    meanIrExp2(i) = meanIr;
    sdIrExp2(i) = sdIr;
end

% Plot Experiment 2: irregular synapses vs angle noise
figure;
hold on;
errorbar(angleVals, meanIrExp2, sdIrExp2, '-o', 'LineWidth', 1.5);
errorbar(angleVals, meanRegExp2, sdRegExp2, '-o', 'LineWidth', 1.5);
xlabel('Irregular spine angle noise (SD, deg)');
ylabel('Mean potential synapses (irregular)');
title(sprintf('Experiment 2: synapses vs angle noise (axon = %.1f μm)', ...
    fixedAxonDistance));
legend({'Irregular','Regular'}, 'Location', 'best');
grid on;

end 


%% ------------------------ HELPER FUNCTIONS ------------------------------

function [meanRegSyn, sdRegSyn, meanIrSyn, sdIrSyn] = run_geometric_model(p)
% RUN_GEOMETRIC_MODEL

nRuns = p.nRuns;

regularSynCount = zeros(nRuns,1);
irregularSynCount = zeros(nRuns,1);

for r = 1:nRuns
    % Regular dendrite
    [regBases, regHeads] = generate_regular_spines( ...
        p.dendLength, p.nSpines, p.regSpineLength, p.regSpineAngle);

    [regCount, ~, ~, ~] = evaluate_synapses( ...
        regBases, regHeads, p.axonDistance, p.synapseRadius, p.minCosAlignment);

    regularSynCount(r) = regCount;

    % Irregular dendrite
    [irBases, irHeads] = generate_irregular_spines( ...
        p.dendLength, p.nSpines, ...
        p.irPosJitterFrac, p.irLengthMean, p.irLengthStd, ...
        p.irAngleStd);

    [irCount, ~, ~, ~] = evaluate_synapses( ...
        irBases, irHeads, p.axonDistance, p.synapseRadius, p.minCosAlignment);

    irregularSynCount(r) = irCount;
end

meanRegSyn = mean(regularSynCount);
sdRegSyn = std(regularSynCount);
meanIrSyn = mean(irregularSynCount);
sdIrSyn = std(irregularSynCount);

end


%% ------------------------ LOCAL COPY OF HELPER FUNCTIONS ------------------------------
function [spineBases, spineHeads] = generate_regular_spines( ...
    dendLength, nSpines, spineLength, spineAngleDeg)

xPositions = linspace(0, dendLength, nSpines)';  
spineBases = [xPositions, zeros(nSpines,1), zeros(nSpines,1)];

theta = deg2rad(spineAngleDeg);
dirXY = repmat([sin(theta), cos(theta), 0], nSpines, 1);

spineVectors = spineLength * dirXY;
spineHeads = spineBases + spineVectors;

end


function [spineBases, spineHeads] = generate_irregular_spines( ...
    dendLength, nSpines, posJitterFrac, lengthMean, lengthStd, angleStdDeg)

basePositions = linspace(0, dendLength, nSpines)';

meanSpacing = dendLength / (nSpines-1);
maxJitter = posJitterFrac * meanSpacing;

jitter = (2*rand(nSpines,1) - 1) * maxJitter;
xPositions = basePositions + jitter;
xPositions = max(min(xPositions, dendLength), 0);

spineBases = [xPositions, zeros(nSpines,1), zeros(nSpines,1)];

lengths = lengthMean + lengthStd*randn(nSpines,1);
lengths(lengths <= 0) = lengthMean * 0.5;

angleNoise = angleStdDeg * randn(nSpines,1);
theta = deg2rad(angleNoise);

dirXY = [sin(theta), cos(theta), zeros(nSpines,1)];
spineVectors = dirXY .* lengths;

spineHeads = spineBases + spineVectors;

end


function [countSyn, meanCosAlign, alignedIdx, withinDist] = evaluate_synapses( ...
    spineBases, spineHeads, axonDistance, synapseRadius, minCosAlignment)

nSpines = size(spineBases,1);

headY= spineHeads(:,2);
distToAxon = abs(headY - axonDistance);
withinDist = distToAxon <= synapseRadius;

spineVecs = spineHeads - spineBases;
baseToAxon = spineBases;
baseToAxon(:,2) = axonDistance - spineBases(:,2);
baseToAxon(:,1) = 0;

spineNorm = vecnorm(spineVecs,2,2);
axonNorm = vecnorm(baseToAxon,2,2);
dotProd = sum(spineVecs .* baseToAxon, 2);

cosTheta = dotProd ./ (spineNorm .* axonNorm + eps);

validIdx = withinDist & ~isnan(cosTheta) & ~isinf(cosTheta);
alignedIdx = validIdx & (cosTheta >= minCosAlignment);

countSyn = sum(alignedIdx);

if any(validIdx)
    meanCosAlign = mean(cosTheta(validIdx));
else
    meanCosAlign = NaN;
end

end
