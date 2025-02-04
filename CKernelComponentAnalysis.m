classdef CKernelComponentAnalysis < handle
    properties
        mPrecursor = -1 % mz of the precursor
        
        mMassAccuracy = 0.005; % dalton
        
        mRawMassMap = []; % mRawMassMap[time, massID]
        mResampledMassMap = []; % mRawMassMap[time, massID]
        mFilteredMassMap = []; % mRawMassMap[time, massID]
        
        mMasses = []; % the mz list, corresponds to mRawMassMap and mResampledMassMap
        mFilteredMasses = []; % the filtered mz list, corresponds to mFilteredMassMap
        mFilteredToOriginalMapping = [];

        mRetentionTime = [];
        mResampledRetentionTime = [];

        mMassMapNormarlizationFactor = [];
        mSourceNum = 0; % the number of components (i.e, sources)
        mSources = []; % mSources.mean and mSources.std
        mSourceWeights = []; % mSourceWeights[massID, sourceNum]
        mSourceWeightMask = []; % mSourceWeightMask[massID, sourceNum]
        mSourceWeightThreshold = []; % [1, mSourceNum]
        mSourceStdThreshold = []; % Maximum std;
        mWeightLambda = 1; % L1-regularization
        mEstimation = []; % mEstimation[time, massID]
        mCost = Inf;
        
        mRawData = [];
        mFilteredSpectraRange = []; % [startIdx, endIdx] of the kept spectra
        mNoiseThreshold = 1;
        mSourceData = {};
        
        mComment = '';
        mDefaultPrecursorDuration = 5;
    end
    
    methods
        function clear( obj )
            obj.mPrecursor = -1;
            obj.mComment = '';
            obj.mRawMassMap = [];
            obj.mMasses = [];
            obj.mSources = [];
            obj.mSourceWeights = [];
            obj.mSourceWeightMask = [];
            obj.mSourceWeightThreshold = [];
            obj.mEstimation = [];
            obj.mRawData = [];
            obj.mSourceStdThreshold = [];
            obj.mRetentionTime = [];
            obj.mResampledRetentionTime = [];
            obj.mFilteredMassMap = [];
        end
        
        function e = analyze(obj, numSources, maxIteration, weightLambda, denoise, verbose)
            if nargin < 2 || isempty(numSources), numSources = 2; end
            if nargin < 3 || isempty(maxIteration), maxIteration = 50; end
            if nargin < 4, weightLambda = []; end
            if nargin < 5, denoise = 0; end
            if nargin < 6, verbose = 0; end
            
            obj.mSourceNum = numSources;
            if denoise
                obj.filter_noise(); % This is an important step.
            end
            obj.analyze_init( weightLambda );
            
            disp( 'Fitting the model ...' );
            e = zeros(1, maxIteration);
            for it = 1 : maxIteration
                old.mSources = obj.mSources;
                old.mSourceWeights = obj.mSourceWeights;
                obj.analyze_adjust_weights();
                obj.analyze_adjust_components();
                obj.analyze_cost();
                e(it) = obj.mCost;
                if verbose
                    disp( obj.mCost );
                    disp( [obj.mSources.mean, obj.mSources.var] );
                    disp( sum( obj.mSourceWeights ) );
                end
                if it > 1 && e(it) > e(it-1) % Stop and roll back if not converge (usually learning rate is too large).
                    e = e(1:it);
                    obj.mSources = old.mSources;
                    obj.mSourceWeights = old.mSourceWeights;
                    obj.mEstimation = (obj.mSourceWeights * obj.mSources.profile)';
                    break;
                end
            end
            if verbose
                plot(e);
            end
            
            % Improve fitting
            disp( 'Improving fitting ...' );
            obj.mSourceWeightMask = obj.mSourceWeights > 5;
            obj.analyze_improve_fitting(maxIteration, 0);
            obj.analyze_filter_massmap_outlier(); % Remove bad fits.
            %obj.analyze_improve_component_by_diagnosis_peaks;
            %obj.analyze_adjust_weights;
            disp( 'Done' );
        end
        
        function profile = analyze_init(obj, weightLambda)
            %profile = obj.mRawMassMap(:, end);
            %profile = sgolayfilt( obj.mRawMassMap(:, end), 4, 7);
            %profile = medfilt1( obj.mRawMassMap(:, end), 5);
            %profile = sum( medfilt1(obj.mRawMassMap(:, 1:end-1)), 2 );
            %profile = sum( medfilt1(obj.mRawMassMap(:, 1:end-1) .* (obj.mRawMassMap(:, 1:end-1) > 10^7)), 2 );
            %profile = sum( obj.mRawMassMap, 2);
            profile = mean(obj.mFilteredMassMap(:, 1:end-1), 2);
            %profile = medfilt1( mean(obj.mRawMassMap(:, 1:end-1), 2 ) );
            %profile( bwareaopen(profile < max(profile)/10, 5) ) = 0; % remove noise
            %profile = mean( obj.mRawMassMap, 2);
            % profile = mean( obj.mRawMassMap(:, sum(obj.mRawMassMap)>10^4), 2 ); % Not perfect.
            
            % Estimate GMM. Each mixture corresponds to a glycan.
            a = CKernelComponentAnalysis.resample_profile(profile, obj.mSourceNum*2000);
            warning('off', 'all');
            gm = fitgmdist( a', obj.mSourceNum, 'Options', statset('MaxIter',500), 'Replicates', 10 );
            warning('on');
            [obj.mSources.mean, idx] = sort(gm.mu);
            obj.mSources.var = squeeze(gm.Sigma(idx)); 
            obj.mSources.std = sqrt(obj.mSources.var);
            
            spectrumLength = size(obj.mFilteredMassMap,1);
            x = 1 : spectrumLength;
            obj.mSources.profile = zeros(obj.mSourceNum, spectrumLength);
            for k = 1 : obj.mSourceNum
                obj.mSources.profile(k,:) = exp( -(x - obj.mSources.mean(k)) .^ 2 / (2 * obj.mSources.var(k) ) ) / obj.mSources.std(k);
            end
            if nargin < 2 || isempty( weightLambda )
                obj.mWeightLambda = max( obj.mRawMassMap(:) ) / 1000; % Initialize with heavy penalty to damp noise.
            else
                obj.mWeightLambda = weightLambda;
            end
            obj.mSourceWeights = zeros(length(obj.mFilteredMasses), obj.mSourceNum);
            obj.mSourceWeightMask = [];
            obj.mSources.method = 'KCA';
            
            if obj.mSourceNum > 2
                temp = median( obj.mSources.std );
                r = (obj.mSources.std - temp) / temp;
                flag = r > 0.333;
                if any( flag )
                    temp = obj.mSources.std( ~flag );
                    obj.mSourceStdThreshold = max(temp) * 1.2;
                    obj.mSources.std( flag ) = max(temp) * 1.05;
                    obj.mSources.var( flag ) = obj.mSources.std( flag ) .^ 2;
                    for k = 1 : obj.mSourceNum
                        if ~flag(k), continue; end
                        obj.mSources.profile(k,:) = exp( -(x - obj.mSources.mean(k)) .^ 2 / (2 * obj.mSources.var(k) ) ) / obj.mSources.std(k);
                    end
                end
            end
        end
        
        function analyze_adjust_components(obj)
            diff = obj.mFilteredMassMap - obj.mEstimation; % diff[time, massID]
            t = 1 : size(obj.mFilteredMassMap, 1);
            dt = t - obj.mSources.mean;
            dt2 = dt .^ 2;
            gM = zeros(obj.mSourceNum, 1);
            gV = zeros(obj.mSourceNum, 1);
            for s = 1 : obj.mSourceNum
                gM(s) = -(obj.mSources.profile(s,:) .* dt(s,:) / obj.mSources.var(s)) * ...
                    diff * obj.mSourceWeights(:,s);
                gV(s) = (obj.mSources.profile(s,:) .* (obj.mSources.var(s) - dt2(s,:)) / obj.mSources.std(s)^3) * ...
                    diff * obj.mSourceWeights(:,s);
            end
            scale = 0.01 / max(max(abs(gM)), 0.01);
            obj.mSources.mean = obj.mSources.mean - gM * scale;
            obj.mSources.var = obj.mSources.var - gV * scale;
            obj.mSources.std = sqrt( obj.mSources.var );
            if ~isempty( obj.mSourceStdThreshold )
                flag = obj.mSources.std > obj.mSourceStdThreshold;
                obj.mSources.std(flag) = obj.mSourceStdThreshold;
                obj.mSources.var(flag) = obj.mSourceStdThreshold ^ 2;
            end
            for k = 1 : obj.mSourceNum
                obj.mSources.profile(k,:) = exp( -(t - obj.mSources.mean(k)) .^ 2 / (2 * obj.mSources.var(k) ) ) / obj.mSources.std(k);
            end
        end
        
        function analyze_adjust_weights(obj)
            H = zeros(obj.mSourceNum, obj.mSourceNum);
            for k = 1 : obj.mSourceNum
                H(k,k) = obj.mSources.profile(k,:) * obj.mSources.profile(k,:)';
                H(k,k) = H(k,k) + 1e-6; % avoid singular
                for m = k+1 : obj.mSourceNum
                    H(k,m) = obj.mSources.profile(k,:) * obj.mSources.profile(m,:)';
                    H(m,k) = H(k,m);
                end
            end
            invH = inv(H);
            spectrumLength = size(obj.mFilteredMassMap, 2);
            b = zeros(obj.mSourceNum,1);
            for m = 1 : spectrumLength  % loop through m/z
                if ~isempty( obj.mSourceWeightMask ) && sum(obj.mSourceWeightMask(m,:)) == 0
                % Noise is marked by mSourceWeightMask
                    obj.mSourceWeights(m,:) = 0;
                    continue;
                end
                for a = 1 : obj.mSourceNum
                    b(a) = obj.mSources.profile(a, :) * obj.mFilteredMassMap(:, m);
                end
                b = b - obj.mWeightLambda;
                temp = invH * b;
                temp(temp < 0) = 0;
                if isempty( obj.mSourceWeightMask )
                    obj.mSourceWeights(m,:) = temp;
                else
                    obj.mSourceWeights(m,:) = temp' .* obj.mSourceWeightMask(m,:);
                end
            end
            obj.mEstimation = (obj.mSourceWeights * obj.mSources.profile)';
        end
        
        function analyze_cost(obj)
            obj.mCost = sum(sum( (obj.mFilteredMassMap - obj.mEstimation).^2 )) + ...
                                  obj.mWeightLambda * sum( obj.mSourceWeights(:) );
        end
        
        function analyze_filter_massmap_outlier(obj)
        % Detect more outliers missed by filter_noise().
        % Use the estimated source weights without regularization
        % remove those are not approximated well by the sources
            bg = obj.mFilteredMassMap( obj.mEstimation < 1 ); 
            bg = sort(bg(bg>0));
            if length(bg) > 10
                obj.mSourceWeightThreshold = bg(round(length(bg)*0.5)) ./ max(obj.mSources.profile');
            else
                obj.mSourceWeightThreshold = max( obj.mSourceWeights ) / 1000;
            end
            for k = 1 : length(obj.mFilteredMasses)
                obj.mSourceWeightMask(k, :) = obj.mSourceWeights(k,:) > obj.mSourceWeightThreshold;
                if any(obj.mSourceWeightMask(k, :))
                    obj.mSourceWeights(k, :) = obj.mSourceWeights(k, :) .* obj.mSourceWeightMask(k, :);
                    obj.mEstimation(:,k) = obj.mSourceWeights(k,:) * obj.mSources.profile;
                    r = corr( obj.mEstimation(:, k), obj.mFilteredMassMap(:, k) );
                    if r < 0.85
                        obj.mSourceWeightMask(k,:) = 0;
                        obj.mSourceWeights(k,:) = 0;
                        obj.mEstimation(:,k) = 0;
                    end
                else
                    obj.mSourceWeights(k,:) = 0;
                    obj.mEstimation(:,k) = 0;
                end
            end
        end

        function e = analyze_improve_fitting(obj, maxIteration, weightLambda, verbose)
            if nargin < 2 || isempty(maxIteration), maxIteration = 50; end
            if nargin < 3, weightLambda = 0; end
            if nargin < 4, verbose = 0; end

            oldWeightLambda = obj.mWeightLambda;
            obj.mWeightLambda = weightLambda;
            e = zeros(1, maxIteration);
            for it = 1 : maxIteration
                old.mSources = obj.mSources;
                old.mSourceWeights = obj.mSourceWeights;
                obj.analyze_adjust_weights();
                obj.analyze_adjust_components();
                obj.analyze_cost();
                e(it) = obj.mCost;
                if verbose
                    disp( obj.mCost );
                    disp( [obj.mSources.mean, obj.mSources.var] );
                    disp( sum( obj.mSourceWeights ) );
                end
                if it > 1 && e(it) > e(it-1) % Stop and roll back if not converge (usually learning rate is too large).
                    e = e(1:it);
                    obj.mSources = old.mSources;
                    obj.mSourceWeights = old.mSourceWeights;
                    obj.mEstimation = (obj.mSourceWeights * obj.mSources.profile)';
                    break;
                end
            end
            obj.mWeightLambda = oldWeightLambda;
            if verbose
                plot(e);
            end
        end
        
        function analyze_EM(obj, maxIteration)
            if nargin < 2 || isempty(maxIteration), maxIteration = 100; end
            peakNum = length(obj.mMasses);
            signalLength = size(obj.mRawMassMap, 1);
            signals = zeros(obj.mSourceNum, signalLength);
            psignal = zeros(peakNum, obj.mSourceNum, signalLength);
            threshold = sum(obj.mRawMassMap(:,end))/1000;
            
            for it = 1 : maxIteration
                % estimate p(signal from a source | sources)
                for k = 1 : peakNum-1 % Don't use the precursor profile because its signal is usually 2 order higher and hence affects results decisively.
                    for c = 1 : obj.mSourceNum
                        signals(c,:) = obj.mSourceWeights(k,c) * obj.mSources.profile(c,:);
                    end
                    totalSignal = sum(signals) + eps;
                    for c = 1 : obj.mSourceNum
                        psignal(k,c,:) = signals(c,:) ./ totalSignal;
                    end
                end
                
                % update sources using p(signal from a source | sources) and signals
                x = 1 : signalLength;
                for c = 1 : obj.mSourceNum
                    m = 0;
                    v = 0;
                    totalW = eps;
                    for k = 1 : peakNum-1 % Don't use the precursor profile because its signal is usually 2 order higher and hence affects results decisively.
                        t = sum(obj.mRawMassMap(:, k));
                        if t < threshold, continue; end
                        w = squeeze(psignal(k,c,:)) .* obj.mRawMassMap(:, k) / t;
                        totalW = totalW + sum( w );
                        % temp = squeeze(psignal(k,c,:)) .* obj.mRawMassMap(:, k) / totalS;
                        %m = m + sum( temp .* obj.mRawMassMap(:, k) );
                        m = m + sum( w .* x' );
                        %v = v + sum( temp .* obj.mRawMassMap(:, k).^2 );
                        v = v + sum( w .* x'.^2 );
                    end
                    m = m / totalW;
                    obj.mSources.mean(c) = m;
                    obj.mSources.var(c) = v / totalW - m*m;
                    obj.mSources.std(c) = sqrt( obj.mSources.var(c) );
                    obj.mSources.profile(c,:) = exp( -(x - obj.mSources.mean(c)) .^ 2 / (2 * obj.mSources.var(c) ) ) / obj.mSources.std(c);
                end
                obj.analyze_adjust_weights();
            end
        end
        
        function analyze_improve_KCA_by_diagnosis_peaks(obj)
            dpeaks = obj.detect_diagnosis_peaks;
            
            % update sources using p(signal from a source | sources) and signals
            x = (1 : size(obj.mFilteredMassMap,1))';
            for c = 1 : obj.mSourceNum
                if isempty(dpeaks{c}), disp(['AICD skips ', num2str(c), '-th component']); continue; end
                m = 0;
                totalW = eps;
                profiles = medfilt1(obj.mFilteredMassMap(:, [dpeaks{c}.peak]));
                for k = 1 : length(dpeaks{c})
                    temp = (dpeaks{c}(k).prob > 0.75)' .* profiles(:, k);
                    mv = max(temp);
                    idx = mean(find( temp == mv ));
                    w = dpeaks{c}(k).prob * profiles(:,k);
                    totalW = totalW + w;
                    m = m + idx * w;
                end
                m = m / totalW;
                obj.mSources.mean(c) = m;

                leftV = 0; rightV = 0;
                leftW = eps; rightW = eps;
                leftM = floor(m); rightM = ceil(m);
                for k = 1 : length(dpeaks{c})
                    temp = (dpeaks{c}(k).prob > 0.5)' .* profiles(:, k);
                    leftW = leftW + sum( temp(1:leftM) );
                    rightW = rightW + sum( temp(rightM:end) );
                    leftV = leftV + sum( (x(1:leftM) - m).^2 .* temp(1:leftM));
                    rightV = rightV + sum( (x(rightM:end) - m).^ 2 .* temp(rightM:end));
                end
                obj.mSources.var(c) = min(leftV/leftW, rightV/rightW);
                %obj.mSources.var(c) = (leftV/leftW + rightV/rightW)/2;
                obj.mSources.std(c) = sqrt( obj.mSources.var(c) );
                
                obj.mSources.profile(c,:) = exp( -(x - obj.mSources.mean(c)) .^ 2 / (2 * obj.mSources.var(c) ) ) / obj.mSources.std(c);
            end
            obj.analyze_adjust_weights;
        end

        function analysis_PCA( obj, componentNum )
            % call obj.filter_noise(); before calling this function. 
            % This is an important step.
            [components, scores, variance] = pca( obj.mFilteredMassMap' );
            if nargin < 2
                temp = cumsum(variance) / sum(variance);
                for componentNum = 1 : length(temp)
                    if temp(componentNum) > 0.99
                        break;
                    end
                end
            end
            obj.mSources.method = 'PCA';
            obj.mSources.mean = [];
            obj.mSources.var = [];
            obj.mSources.std = [];
            obj.mSources.profile = components(:, 1:componentNum)';
            obj.mSourceWeights = scores(:, 1:componentNum);
        end
        
        function analysis_NNMF( obj, componentNum )
            [W, H] = nnmf( obj.mFilteredMassMap, componentNum );
            obj.mSources.mean = [];
            obj.mSources.var = [];
            obj.mSources.std = [];
            obj.mSources.profile = W';
            obj.mSourceWeights = H';
            obj.mSources.method = 'NNMF';
        end
        
        function [idxes] = compare_estimated_and_source( obj, idxes, timeRange, filtered )
            if isempty( obj.mEstimation ) || isempty( obj.mFilteredMassMap )
                return;
            end
            if nargin < 3, timeRange = []; end
            if isempty( timeRange ), timeRange = [1, size(obj.mResampledMassMap, 1)]; end
            if nargin < 4 || isempty(filtered)
                filtered = 1;
            end

            scanIdx = max(obj.mFilteredSpectraRange(1), timeRange(1)) : min(obj.mFilteredSpectraRange(2), timeRange(2));
            x_marks = obj.mResampledRetentionTime(scanIdx);
                
            if nargin < 2 || isempty( idxes ) % find the ions to be visualized
                idxes = [];
                for s = 1 : obj.mSourceNum
                    temp = find(obj.mSourceWeights(:,s) > 1);
                    idxes = [idxes; temp];
                end
            end
            idxes = unique(idxes);
            
            colors = 'rgmbcy';
            
            % Use the original signals
            ori = obj.mResampledMassMap(:, obj.mFilteredToOriginalMapping(idxes));
            if filtered
                ori = medfilt1(ori);
            end

            % Use the medfilt1 signals
            % ori = obj.mFilteredMassMap(:, idxes);
            % scanIdx = 1 : size(ori, 1);

            est = obj.mEstimation(:, idxes);
            sourceWeights = obj.mSourceWeights(idxes,:);
            selectedMasses = obj.mFilteredMasses(idxes);
            figure;
            num = length(idxes);
            temp = sqrt(num / 10);
            rows = max(floor( temp * 3 ), 1);
            cols = ceil( num/rows );

            legends = { 'Original signal', 'Approximated by kernels' };
            for c = 1 : obj.mSourceNum
                legends{end+1} = ['Kernel ', num2str(c)];
            end
            for k = 1 : num
                subplot(rows, cols, k);
                plot( x_marks, ori(scanIdx, k), 'k', 'LineWidth', 2 ); hold on;
                plot( x_marks, est(:, k), 'k-.', 'LineWidth', 4 );
                cidx = 1;
                for c = 1 : obj.mSourceNum
                    temp = sourceWeights(k,c) * obj.mSources.profile(c,:);
                    plot( x_marks, temp, colors(cidx), 'LineWidth', 2 );
                    cidx = cidx + 1;
                    if cidx > length(colors), cidx = length(colors); end
                end
                hold off;
                axis tight;
                title( [num2str(idxes(k)), ': ', num2str(selectedMasses(k), '%.3f')] );
                if k == 1
                    legend( legends, 'FontSize', 9 );
                end
            end
        end
        
        function diagnosisPeaks = detect_diagnosis_peaks( obj )
        % Detect strong enough diagnosis peaks of each component.
            diagnosisPeaks = cell(1, obj.mSourceNum);
            numMasses = length( obj.mFilteredMasses );
            signalLength = size(obj.mFilteredMassMap, 1);
            signals = zeros(obj.mSourceNum, signalLength);
            
            widthThreshold = obj.mSources.std * 2;
            startPos1 = obj.mSources.mean - 1.75*obj.mSources.std;
            endPos1 = obj.mSources.mean + 1.75*obj.mSources.std;
            startPos2 = obj.mSources.mean - 2 * obj.mSources.std;
            endPos2 = obj.mSources.mean + 2 * obj.mSources.std;
            centerPos = round(obj.mSources.mean);
            for k = 1 : numMasses-1
                for c = 1 : obj.mSourceNum
                    signals(c,:) = obj.mSourceWeights(k,c) * obj.mSources.profile(c,:);
                end
                totalSignal = sum(signals) + eps;
                threshold = max(totalSignal) / 2;
                fgflag = max(signals,[], 2) > threshold;
                if sum(fgflag) >= 3, continue; end
                for c = 1 : obj.mSourceNum
                    %if max(signals(c,:)) < threshold
                    if fgflag(c) == 0
                        continue;
                    end
                    psource = signals(c,:) ./ totalSignal;
                    mask = (psource .* (obj.mFilteredMassMap(:, k)' > 0) )  > 0.5;
                    idx = centerPos(c);
                    while idx > 0 && mask(idx) > 0; idx = idx - 1; end
                    mask(1:idx) = 0;
                    idx = centerPos(c);
                    while idx < signalLength && mask(idx) > 0; idx = idx + 1; end
                    mask(idx:end) = 0;
                    if sum(mask) >= widthThreshold(c)
                        positions = find(mask);
                        if (positions(1) <= startPos2(c) &&  positions(end) > endPos1(c)) || ...
                                (positions(end) >= endPos2(c) &&  positions(1) < startPos1(c))
                            diagnosisPeaks{c}(end+1) = struct('peak', k, 'prob', psource, ...
                                'reliableStart', positions(1), 'reliableEnd', positions(end) );
                        end
                    end
                end
            end
        end
        
        function decompose_spectrum(obj)
            obj.mSourceData = cell(1, obj.mSourceNum);
            for s = 1 : obj.mSourceNum
                obj.mSourceData{s} = cell(1, size(obj.mRawMassMap,1));
                temp = obj.mSourceWeights(:,s) * obj.mSources.profile(s,:);
                for k = 1 : size(obj.mRawMassMap,1)
                    obj.mSourceData{s}{k} = temp(:,k)';
                end
            end
        end
        
        function filter_noise( obj ) 
            numSpectra = length(obj.mRawData);
            numMasses = length(obj.mMasses);
            
            obj.mFilteredMassMap = medfilt1(obj.mResampledMassMap);
            flag = ones(1, numMasses);
            fgTH = max(max(obj.mFilteredMassMap(:, 1:end-1))) / 100; % Tried and Failed: using mean and std of all signals.
            mask = zeros(numSpectra, 1);
            for m = 1 : numMasses
                temp = obj.mFilteredMassMap(:, m);
                if max(temp) < fgTH
                    flag(m) = 0;
                    continue;
                end
                mask = mask + bwareaopen(temp > max(temp)/10, 5);
            end
            mask = mask > round(numSpectra * 0.05);
            flag = flag > 0;
            obj.mFilteredMasses = obj.mMasses(flag);
            obj.mFilteredMassMap = obj.mFilteredMassMap(:, flag) .* (mask * ones(1, sum(flag)));
            sIdx = max( 1, find(mask, 1, 'first')-10 );
            eIdx = min( find(mask, 1, 'last') + 10, size(obj.mResampledMassMap,1));
            obj.mFilteredSpectraRange = [sIdx, eIdx];
            obj.mFilteredMassMap = obj.mFilteredMassMap(sIdx:eIdx, :);
            
            obj.mFilteredToOriginalMapping = zeros(1, length(obj.mFilteredMasses));
            idx = 1;
            for k = 1 : length(flag)
                if flag(k)
                    obj.mFilteredToOriginalMapping(idx) = k;
                    idx = idx + 1;
                end
            end
        end
        
        function normalize_massmap(obj)
        % Normalize intensity so that the maximum of each mass is 1.
        % But this step may hurt the result (need better way to handle noise).
            massNum = length(obj.mMasses);
            obj.mMassMapNormarlizationFactor = zeros(1, massNum);
            for k = 1 : massNum
                mv = max( obj.mRawMassMap(:,k) );
                if mv < 1, continue; end
                obj.mMassMapNormarlizationFactor(k) = mv;
                obj.mRawMassMap(:,k) = obj.mRawMassMap(:,k) / mv;
            end
        end

        function [weights, estimation, error] = factorize_spectra(obj, spectra)
        % each row of "spectra" is a spectrum.
        % spectra = weights * sources
        % spectra * sources' = weights * sources * sources'
        % spectra * sources' * inv( sources * sources' ) = weights
            H = obj.mSources.profile * obj.mSources.profile';
            for k = 1 : obj.mSourceNum
                H(k,k) = H(k,k) + 1e-6; % avoid singular
            end
            invH = inv(H);
            weights = spectra * obj.mSources.profile' * invH;
            weights(weights<0) = 0;
            estimation = weights * obj.mSources.profile;
            error = spectra - estimation;
        end
        
        function load( obj, filename )
            obj.clear();
            
            % load peaks
            metal = '';
            fid = fopen( [filename, '.txt'], 'r' );
            tline = fgetl(fid);
            while isempty(tline) || tline(1) == '#'
                if contains(tline, '# Metal:')
                    metal = strtrim( tline(9:end) );
                end
                tline = fgetl(fid);
            end
            data = textscan( fid, '%f%s%f%*s' );
            fclose( fid );
            
            % calculate the precuror 
            [~, pidx] = max( data{3} );
            charge = str2double(data{2}{pidx}(1:end-1));
            if charge == 1
                obj.mPrecursor = data{3}(pidx);
            else
                temp = CPeak.protonate_raw_mz(data{1}(pidx), charge, metal);
                tempIdx = find( temp < data{1}(end), 1, 'last' );
                obj.mPrecursor = temp(tempIdx);
            end
            %obj.mPrecursor = obj.mMasses(end); % this will be adjusted later
            
            % remove mz heavier than the precursor
            max_mass = obj.mPrecursor + obj.mMassAccuracy;
            flag = ones(1, length(data{1}));
            for k = 1 : length(data{1})
                if data{1}(k) > max_mass
                    flag(k) = 0;
                    continue;
                end
                charge = str2double(data{2}{k}(1:end-1));
                if charge > 1
                    pm = CPeak.protonate_raw_mz( data{1}(k), charge, metal );
                    if ~any(pm <= max_mass)
                        flag(k) = 0;
                    end
                elseif data{1}(k) > max_mass
                    flag(k) = 0;
                end
            end
            obj.mMasses = data{1}(flag>0);
            
            % load spectra
            disp( ['Loading ', filename, '.mzXML ...'] );
            mzxml_struct = mzxmlread( [filename, '.mzXML'] );
            [spectra, retentionTime] = mzxml2peaks(mzxml_struct, 'Levels', 1 );
            maxLevel = max([mzxml_struct.scan.msLevel]);
            if isempty( spectra )
                [spectra, retentionTime] = mzxml2peaks(mzxml_struct, 'Levels', 2 );
                precursorRetentionTime = retentionTime;
                numSpectra = length(spectra);
                precursorProfile = zeros(numSpectra, 1);
                for s = 1 : numSpectra
                    idxes = find( abs(spectra{s}(:,1) - obj.mPrecursor) < obj.mMassAccuracy );
                    if ~isempty(idxes)
                        precursorProfile(s) = max( spectra{s}(idxes, 2) );
                    end
                end
            elseif maxLevel > 1
                % collect precursor profile
                numSpectra = length(spectra);
                precursorProfile = zeros(numSpectra, 1);
                precursorRetentionTime = retentionTime;
                for s = 1 : numSpectra
                    idxes = find( abs(spectra{s}(:,1) - obj.mPrecursor) < obj.mMassAccuracy );
                    if ~isempty(idxes)
                        precursorProfile(s) = max( spectra{s}(idxes, 2) );
                    end
                end
                
                [spectra, retentionTime] = mzxml2peaks(mzxml_struct, 'Levels', 2 );
                numSpectra = length(spectra);
            end
            obj.mRawData = spectra;
            obj.mRetentionTime = retentionTime;
            numMasses = length(obj.mMasses);
            obj.mRawMassMap = zeros(numSpectra, numMasses);
            obj.mResampledMassMap = zeros(numSpectra, numMasses);

            % Correct baseline & create mRawMassMap
            disp( 'Creating mass map ...' );
            for s = 1 : numSpectra
                flag = (spectra{s}(:,1) < max_mass) & (spectra{s}(:,1) > 100);
                spectra{s} = spectra{s}(flag, :);
                %temp = msbackadj( spectra{s}(:,1), spectra{s}(:,2) );
                %temp(temp < 0)= 0;
                %spectra{s}(:,2) = temp;
                
                for m = 1 : numMasses
                    % [minV, idx] = min( abs(spectra{s}(:,1) - obj.mMasses(m)) );
                    % if minV <= obj.mMassAccuracy
                    %    obj.mRawMassMap(s, m) = spectra{s}(idx, 2);
                    % end
                    intensities = spectra{s}(abs(spectra{s}(:,1) - obj.mMasses(m)) < obj.mMassAccuracy, 2);
                    if ~isempty( intensities )
                        %obj.mRawMassMap(s, m) = max(intensities);
                        obj.mRawMassMap(s, m) = sum(intensities);
                    end
                end
            end
            
            % Resample temporal profiles.
            disp( 'Resampling accordingly to the retention time ...' );
            step = (max(retentionTime) - min(retentionTime)) / (numSpectra-1);
            newTimeStamps = min(retentionTime) : step : max(retentionTime);
            for m = 1 : numMasses
                tsin = timeseries( obj.mRawMassMap(:,m), retentionTime);
                tsout = resample(tsin, newTimeStamps);
                obj.mResampledMassMap(:, m) = tsout.Data(:);
            end
            if maxLevel > 1
                tsin = timeseries(precursorProfile, precursorRetentionTime);
                tsout = resample(tsin, newTimeStamps);
                obj.mResampledMassMap(:, end) = tsout.Data(:);
            end
            
            obj.mResampledRetentionTime = newTimeStamps;
            
            obj.mComment = filename;
        end
        
        function visualize_sources( obj )
            plot( obj.mFilteredSpectraRange(1):obj.mFilteredSpectraRange(2), obj.mSources.profile' );
            ax = gca;
            if ~isempty(num2str(obj.mSources.mean))
                for s = 1 : obj.mSourceNum
                    temp = obj.mSources.mean(s) + obj.mFilteredSpectraRange(1) - 1;
                    label = [num2str(temp, '%4.2f'), ' (', num2str(obj.mSources.std(s), '%4.2f'), ')'];
                    text( temp, max(obj.mSources.profile(s,:)), label );
                end
            end
            title( obj.mSources.method );
            set(ax, 'XLimSpec', 'Tight');
        end % function visualize_sources
        
    end % method
    
    methods (Static)
        function samples = resample_profile( profile, sampleNum )
            p = profile / sum(profile);
            sp = cumsum(p);
            rs = rand(1, sampleNum);
            samples = zeros(1, sampleNum);
            for k = 1 : sampleNum
                samples(k) = find( sp >= rs(k), 1, 'first' );
            end
            
%             num = round( sampleNum * p );
%             diff = sampleNum - sum(num);
%             
%             if diff ~= 0
%                s = sign(diff);
%                diff = abs(diff);
%                k = 1;
%                while k <= diff
%                     a = rand;
%                     idx = find(sp <= a, 1, 'last' );
%                     if s < 0
%                         if num(idx) == 0
%                             continue;
%                         else
%                             num(idx) = num(idx) - 1;
%                         end
%                     else
%                         num(idx) = num(idx) + 1;
%                     end
%                     k = k + 1;
%                 end
%             end
%             
%             samples = [];
%             for k = 1 : length(num)
%                 samples = [samples, ones(1, num(k))*k];
%             end
        end
        
        function kca = process( filename, massAccuracy, sourceNum, weightLambda )
            if nargin < 2 || isempty(massAccuracy)
                massAccuracy = 0.005;
            end
            if nargin < 3 || isempty(sourceNum)
                sourceNum = 2;
            end
            if nargin < 4
                weightLambda = [];
            end
            
            disp( ['Working on ', filename] );
            kca = CKernelComponentAnalysis;
            kca.mMassAccuracy = massAccuracy;
            kca.load( filename );
            fprintf('Detecting peaks ...\n');
            %kca.detect_peaks();
            %kca.detect_precursor();
            %if isempty( kca.mPrecursor )
            %    kca.detect_precursor();
            %    fprintf( 'Detected the precursor as %f\n', kca.mPrecursor );
            %end
            fprintf('Filtering peaks ...\n');
            %kca.filter_noise();
            fprintf('Source decomposition ...\n');
            kca.analyze(sourceNum, 100, weightLambda);
            disp( 'Finished' );
        end
    end
end
