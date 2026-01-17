function spine_model_stage1_1()

% Stage 1: geometric model of synapse opportunities 
% between a dendrite with spines and a nearby axon.
%
% Dendrite: straight segment along x-axis from 0 to dendLength
% Axon: straight segment parallel to dendrite, offset in y
% Condition for a potential synapse:
%   (1) Spine head is within synapseRadius of the axon
%   (2) Spine is sufficiently aligned toward the axon
%
% Compare:
%   (1) Regular dendrite: evenly spaced spines, fixed length, fixed angle
%   (2) Irregular dendrite: jittered positions, variable lengths & angles

%% ------------------------- PARAMETERS -----------------------------

% Geometry of dendrite and axon
dendLength = 100; % total dendrite length
nSpines = 40; % number of spines
axonDistance = 2.0; % distance in y from dendrite to axon (parallel line), change later if don't want parallel line

% Spine geometry (regular)
regSpineLength = 1.5; % spine length
regSpineAngle = 0; % angle deviation from perfect pointing at axon

% Spine geometry (irregular)
irPosJitterFrac = 0.25; % fraction of mean spacing used as positional jitter
irLengthMean = 1.8; % mean spine length
irLengthStd = 0.7; % SD of spine length (truncated at min>0)
irAngleStd = 20; % SD (deg) around ideal angle: toward axon

% Synapse rule
synapseRadius = 0.4; % max distance from spine head to axon
maxAngleDeg = 45; % max allowed angle (deg) between spine vector and vector to axon

% Simulation 
nRuns = 50; % how many random dendrites to sample for each type

%% --------------- THRESHOLD FOR ALIGNMENT --------------

maxAngleRad = deg2rad(maxAngleDeg);
minCosAlignment = cos(maxAngleRad); % require cos(theta) >= this

%% --------------- STORAGE FOR SUMMARY METRICS ----------------------

regularSynCount = zeros(nRuns,1);
irregularSynCount = zeros(nRuns,1);

regularMeanAlign = nan(nRuns,1);
irregularMeanAlign = nan(nRuns,1);

%% ---------------------- MAIN SIMULATION ---------------------------

for r = 1:nRuns
    % REGULAR DENDRITE
    [regBases, regHeads] = generate_regular_spines( ...
        dendLength, nSpines, regSpineLength, regSpineAngle, axonDistance);

    [regCount, regAlign, regAlignedIdx, regWithinDist] = evaluate_synapses( ...
        regBases, regHeads, axonDistance, synapseRadius, minCosAlignment);

    regularSynCount(r) = regCount;
    regularMeanAlign(r) = regAlign;

    % IRREGULAR DENDRITE
    [irBases, irHeads] = generate_irregular_spines( ...
        dendLength, nSpines, ...
        irPosJitterFrac, irLengthMean, irLengthStd, irAngleStd, axonDistance);

    [irCount, irAlign, irAlignedIdx, irWithinDist] = evaluate_synapses( ...
        irBases, irHeads, axonDistance, synapseRadius, minCosAlignment);

    irregularSynCount(r) = irCount;
    irregularMeanAlign(r) = irAlign;
end

%% ------------------------ PRINT SUMMARY ---------------------------

fprintf('Stage 1 geometric model (nRuns = %d)\n', nRuns);
fprintf('Regular dendrite:\n');
fprintf('  Mean potential synapses: %.2f (SD = %.2f)\n', ...
    mean(regularSynCount), std(regularSynCount));
fprintf('  Mean alignment (cos theta) among candidates: %.3f\n', ...
    mean(regularMeanAlign));

fprintf('\nIrregular dendrite:\n');
fprintf('  Mean potential synapses: %.2f (SD = %.2f)\n', ...
    mean(irregularSynCount), std(irregularSynCount));
fprintf('  Mean alignment (cos theta) among candidates: %.3f\n', ...
    mean(irregularMeanAlign));

%% ------------------------ PLOTS ----------------------------
t = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timeTag = char(t);

% ---------- Figure 1: histograms + bar ----------
fig1 = figure; 
subplot(1,2,1);
histogram(regularSynCount, 'FaceAlpha', 0.7);
hold on;
histogram(irregularSynCount, 'FaceAlpha', 0.7);
xlabel('Potential synapse count');
ylabel('Frequency');
legend({'Regular','Irregular'});
title('Distribution of potential synapses');

subplot(1,2,2);
bar([mean(regularSynCount), mean(irregularSynCount)]);
set(gca,'XTickLabel',{'Regular','Irregular'});
ylabel('Mean potential synapses');
title('Mean across runs');

statsFigName = ['synapse_stats_' timeTag '.png'];
saveas(fig1, statsFigName);

% ---------- Figure 2: geometry examples ----------
allHeadsY = [regHeads(:,2); irHeads(:,2)];
yMaxUnified = max([axonDistance + synapseRadius, max(allHeadsY)])+0.2;
yMaxUnified = ceil(yMaxUnified/0.5)*0.5;
yMinUnified = 0;

fig2 = figure;
subplot(2,1,1);
plot_geometry_example(dendLength, axonDistance, ...
    regBases, regHeads, synapseRadius, regAlignedIdx, regWithinDist);
title('Example: Regular dendrite');

subplot(2,1,2);
plot_geometry_example(dendLength, axonDistance, ...
    irBases, irHeads, synapseRadius, irAlignedIdx, irWithinDist);
title('Example: Irregular dendrite');

subplot(2,1,1);
ylim([yMinUnified yMaxUnified]);
yticks(yMinUnified:0.5:yMaxUnified);

subplot(2,1,2);
ylim([yMinUnified yMaxUnified]);
yticks(yMinUnified:0.5:yMaxUnified);

geomFigName = ['geometry_examples' timeTag '.png'];
saveas(fig2, geomFigName);

%% ------------------------ ANIMATIONS ------------------------------
% Animations showing spines connecting to the axon
regGifName = ['regular_dendrite_' timeTag '.gif'];
irGifName = ['irregular_dendrite_' timeTag '.gif'];

fig3 = figure;
animate_synapse_formation(dendLength, axonDistance, ...
    regBases, regHeads, synapseRadius, regAlignedIdx, regWithinDist, regGifName, yMaxUnified);
title('Regular dendrite: synapse formation');

fig4 = figure;
animate_synapse_formation(dendLength, axonDistance, ...
    irBases, irHeads, synapseRadius, irAlignedIdx, irWithinDist, irGifName, yMaxUnified);
title('Irregular dendrite: synapse formation');

end

%% ------------------------ HELPER FUNCTIONS ------------------------------

function [spineBases, spineHeads] = generate_regular_spines( ...
    dendLength, nSpines, spineLength, spineAngleDeg, axonDistance)
% GENERATE_REGULAR_SPINES

    % Spine base positions along dendrite (x-axis)
    xPositions = linspace(0, dendLength, nSpines)';  
    spineBases = [xPositions, zeros(nSpines,1), zeros(nSpines,1)];

    % Convert uniform angle offset in degrees to radians
    theta = deg2rad(spineAngleDeg);

    % Rotate the baseToAxonDir around z-axis by theta
    dirXY = repmat([sin(theta), cos(theta), 0], nSpines, 1); % unit vector in xy-plane

    spineVectors = spineLength * dirXY; % 1x3 vector
    spineHeads = spineBases + spineVectors; 

end


function [spineBases, spineHeads] = generate_irregular_spines( ...
    dendLength, nSpines, posJitterFrac, lengthMean, lengthStd, ...
    angleStdDeg, axonDistance)
% GENERATE_IRREGULAR_SPINES
% Spines: 
%   base positions jittered around even spacing
%   lengths drawn from N(lengthMean, lengthStd^2), truncated > 0
%   angles drawn from N(0, angleStdDeg^2) around ideal position

    % Base positions
    basePositions = linspace(0, dendLength, nSpines)';
    meanSpacing = dendLength / (nSpines-1);
    maxJitter = posJitterFrac * meanSpacing;

    % Jitter along x (keeping within [0, dendLength])
    jitter = (2*rand(nSpines,1) - 1) * maxJitter; % uniform in [-maxJitter, maxJitter]
    xPositions = basePositions + jitter;
    xPositions = max(min(xPositions, dendLength), 0);

    spineBases = [xPositions, zeros(nSpines,1), zeros(nSpines,1)];

    % Spine lengths from truncated normal
    lengths = lengthMean + lengthStd*randn(nSpines,1);
    lengths(lengths <= 0) = lengthMean * 0.5; % simple truncation for now

    % Angles in degrees around +y direction (ideal toward axon)
    angleNoise = angleStdDeg * randn(nSpines,1);
    theta = deg2rad(angleNoise);

    % For each spine, direction is a unit vector in xy-plane:
    % ideal toward axon = (0,1,0), rotate by theta around z-axis
    % So: dir_i = (sin(theta_i), cos(theta_i), 0)
    dirXY = [sin(theta), cos(theta), zeros(nSpines,1)]; % n x 3
    spineVectors = dirXY .* lengths; % scale each row by length

    spineHeads = spineBases + spineVectors;

end


function [countSyn, meanCosAlign, alignedIdx, withinDist] = evaluate_synapses( ...
    spineBases, spineHeads, axonDistance, synapseRadius, minCosAlignment)
% EVALUATE_SYNAPSES
% return:
%   countSyn: number of spines that satisfy distance + alignment rule
%   meanCosAlign: mean cos(theta) among those that pass the distance rule
%   alignedIdx: logical index of spines that pass both rules
%   withinDist: logical index of spines within synapseRadius

    nSpines = size(spineBases,1);

    % 1) Distance from spine head to nearest point on axon.
    headY = spineHeads(:,2);
    distToAxon = abs(headY - axonDistance);

    % Candidates that are within synapseRadius:
    withinDist = distToAxon <= synapseRadius;

    % 2) Alignment: compare spine vector with vector base -> axon
    spineVecs = spineHeads - spineBases; % n x 3
    baseToAxon = spineBases;
    baseToAxon(:,2) = axonDistance - spineBases(:,2); % adjust y to axon
    baseToAxon(:,1) = 0; % ignore x-direction

    % Compute cos(theta) = (u·v)/(|u||v|)
    spineNorm = vecnorm(spineVecs,2,2); % n x 1
    axonNorm = vecnorm(baseToAxon,2,2); % n x 1
    dotProd = sum(spineVecs .* baseToAxon, 2);

    cosTheta = dotProd ./ (spineNorm .* axonNorm + eps); % avoid /0

    % Apply alignment threshold only to those within distance
    validIdx = withinDist & ~isnan(cosTheta) & ~isinf(cosTheta);
    alignedIdx = validIdx & (cosTheta >= minCosAlignment);

    countSyn = sum(alignedIdx);

    if any(validIdx)
        meanCosAlign = mean(cosTheta(validIdx));
    else
        meanCosAlign = NaN;
    end

end


function plot_geometry_example(dendLength, axonDistance, ...
    spineBases, spineHeads, synapseRadius, alignedIdx, withinDist, yMaxUnified)
% PLOT_GEOMETRY_EXAMPLE
% 2D projection (x vs y) of one dendrite, spines, and axon.

    nSpines = size(spineBases,1);

    if nargin < 6 || isempty(alignedIdx)
        alignedIdx = false(nSpines,1);
    end

    if nargin < 7 || isempty(withinDist)
        withinDist = false(nSpines,1);
    end

    if nargin < 8 || isempty(yMaxUnified)
        allHeadsY   = spineHeads(:,2);
        yMaxUnified = max([axonDistance + synapseRadius, max(allHeadsY)]) + 0.2;
        yMaxUnified = ceil(yMaxUnified/0.5)*0.5;
    end

    hold on; grid on;
    % Dendrite shaft
    plot([0, dendLength], [0, 0], 'k-', 'LineWidth', 2);

    % Axon line
    plot([0, dendLength], [axonDistance, axonDistance], 'r--', 'LineWidth', 2);

    % Spines
    for i = 1:nSpines
        xb = spineBases(i,1); yb = spineBases(i,2);
        xh = spineHeads(i,1); yh = spineHeads(i,2);

        plot([xb, xh], [yb, yh], '-', 'Color','b');

        % Default: blue shaft & head
        headEdge = 'b';
        headFace = 'b';

        % If this spine is a potential synapse, color it green
        if alignedIdx(i)
            headEdge = [0 0.4 0];
            headFace = [0 0.9 0];
        end

        plot(xh, yh, 'o', ...
            'MarkerFaceColor', headFace, ...
            'MarkerEdgeColor', headEdge, ...
            'MarkerSize', 3);

        % Circle showing synapse radius around closest point on axon
        if withinDist(i)
            xClosest = xh;
            yClosest = axonDistance;
            th = linspace(0,2*pi,40);
            xc = xClosest + synapseRadius*cos(th);
            yc = yClosest + synapseRadius*sin(th);
            plot(xc, yc, 'Color',[0.8 0.8 0.8]);
        end

        % If this spine forms a potential synapse: draw green connection
        if alignedIdx(i)
            plot([xh, xh], [yh, axonDistance], ...
                'Color',[0 0.8 0], 'LineWidth', 0.8);
        end
    end

    xlabel('x (along dendrite)');
    ylabel('y');
    axis tight;
    xlim([-5, dendLength+5]);
    ylim([0, yMaxUnified]);
    yticks(0:0.5:yMaxUnified);
    set(gca, 'FontSize', 12);


end


function animate_synapse_formation(dendLength, axonDistance, ...
    spineBases, spineHeads, synapseRadius, alignedIdx, withinDist, gifName, yMaxUnified)
% ANIMATE_SYNAPSE_FORMATION
% show spines one by one, green highlights for those that form potential synapses
%
% alignedIdx: logical vector, true: spine passes distance+alignment rule
% withinDist: logical vector, true: spine is within synapseRadius

    if nargin < 9
        error('animate_synapse_information: not enough input arguments');
    end

    if isempty(gifName)
        saveGIF = false;
    else
        saveGIF = true;
    end

    nSpines = size(spineBases,1);

    clf;
    hold on; grid on;

    % Dendrite shaft
    plot([0, dendLength], [0, 0], 'k-', 'LineWidth', 2);

    % Axon line
    plot([0, dendLength], [axonDistance, axonDistance], 'r--', 'LineWidth', 2);

    xlabel('x (along dendrite)');
    ylabel('y');
    axis tight;
    xlim([-5, dendLength+5]);

    yMin = 0;
    ylim([yMin, yMaxUnified]);
    yticks(yMin:0.5:yMaxUnified);
    set(gca, 'FontSize', 12);
   

    % Animate spine by spine
    for i = 1:nSpines
        xb = spineBases(i,1); yb = spineBases(i,2);
        xh = spineHeads(i,1); yh = spineHeads(i,2);

        % Basic spine shaft + head
        plot([xb, xh], [yb, yh], 'b-');
        plot(xh, yh, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 3);

        % Circle showing synapse radius around closest point on axon
        if withinDist(i)
            xClosest = xh;
            yClosest = axonDistance;
            th = linspace(0,2*pi,40);
            xc = xClosest + synapseRadius*cos(th);
            yc = yClosest + synapseRadius*sin(th);
            plot(xc, yc, 'Color',[0.8 0.8 0.8]);
        end

        % If this spine forms a potential synapse: draw green connection
        if alignedIdx(i)
            plot([xh, xh], [yh, axonDistance], ...
                'Color',[0 0.8 0], 'LineWidth', 0.8);
        
            plot(xh, yh, 'o', ...
                'MarkerSize', 3, ...
                'MarkerEdgeColor',[0 0.4 0], ...
                'MarkerFaceColor',[0 0.9 0]);
        end


        drawnow;
        if saveGIF
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
    
            if i ==1
                imwrite(imind, cm, gifName, 'gif', 'LoopCount', inf, 'DelayTime',0.05);
            else
                imwrite(imind, cm, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
            end
        end
        pause(0.05); % speed
    end

end
