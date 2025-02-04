%%
kca = CKernelComponentAnalysis;
kca.mMasses = IsoMal_Mal_95ms.ions;
kca.mResampledMassMap = IsoMal_Mal_95ms.raw_signals;
kca.mResampledRetentionTime = IsoMal_Mal_95ms.mobility;
kca.mFilteredMassMap = kca.mResampledMassMap;
kca.mFilteredToOriginalMapping = 1 : size(kca.mResampledMassMap, 2);
kca.mFilteredSpectraRange = [1, length( kca.mResampledRetentionTime )];
kca.mFilteredMasses = IsoMal_Mal_95ms.ions;

%%
kca.analyze(4, 500, 10, 0);

%% Compare ions [203.1 365.1 509.2 529.2]
% To use the save kca results, create a kca instance: kca = CKernelComponentAnalysis;
% The KCA results are saved as "/Users/hong/Documents/Projects/Glycomics/data/20240911/kca.mat"
% Run this code segment under "/Users/hong/Documents/Projects/Glycomics/matlab/". Otherwise some functions of kca won't work.
% inds = [2, 10, 15, 16];
ions = [205.1, 245.1, 347.1, 365.1, 367.1, 407.1, 467.2, 185.1, 509.2];
num = length(ions);
mobility = IsoMal_Mal_95ms.mobility;
colors = 'bmrky';
for k = 1 : num
    ionIdx = find( ions(k) == IsoMal_95ms.ions );
    ionname = num2str(IsoMal_95ms.ions(ionIdx));
    subplot(num, 3, (k-1)*3+1);
    plot( mobility, IsoMal_95ms.raw_signals(:,ionIdx), 'k', 'LineWidth', 2 );
    hold on;
    [weights, estimation] = kca.factorize_spectra( IsoMal_95ms.raw_signals(:,ionIdx)' );
    weights(4) = 0;
    plot( mobility, estimation, 'k-.', 'LineWidth', 4 );
    for m = 1 : length( weights )
        plot( mobility, kca.mSources.profile(m, :) * weights(m), colors(m), 'LineWidth', 1.5 );
    end
    hold off;
    % title( ['IsoMal (95ms): ion ', ionname] );
    title( ['IsoMal: ion ', ionname] );

    ionIdx = find( ions(k) == Mal_95ms.ions );
    subplot(num, 3, (k-1)*3+2);
    plot( mobility, Mal_95ms.raw_signals(:, ionIdx), 'k', 'LineWidth', 2 );
    hold on;
    [weights, estimation] = kca.factorize_spectra( Mal_95ms.raw_signals(:, ionIdx)' );
    weights(3) = 0;
    plot( mobility, estimation, 'k-.', 'LineWidth', 4 );
    for m = 1 : length( weights )
        plot( mobility, kca.mSources.profile(m, :) * weights(m), colors(m), 'LineWidth', 1.5 );
    end
    hold off;
    % title( ['Mal (95ms): ion ', ionname] );
    title( ['Mal: ion ', ionname] );

    ionIdx = find( ions(k) == IsoMal_Mal_95ms.ions );
    subplot(num, 3, (k-1)*3+3);
    plot( mobility, IsoMal_Mal_95ms.raw_signals(:, ionIdx), 'k', 'LineWidth', 2 );
    hold on;
    [weights, estimation] = kca.factorize_spectra( IsoMal_Mal_95ms.raw_signals(:, ionIdx)' );
    plot( mobility, estimation, 'k-.', 'LineWidth', 4 );
    for m = 1 : length( weights )
        plot( mobility, kca.mSources.profile(m, :) * weights(m), colors(m), 'LineWidth', 1.5 );
    end
    hold off;
    % title( ['IsoMal & Mal (95ms): ion ', ionname] );
    title( ['IsoMal & Mal: ion ', ionname] );
end

legend( {'Original', 'Approximation', 'Kernel 1', 'Kernel 2', 'Kernel 3', 'Kernel 4'} )

%% Check ions in mixture
% To use the save kca results, create a kca instance: kca = CKernelComponentAnalysis;
% The KCA results are saved as "/Users/hong/Documents/Projects/Glycomics/data/20240911/kca.mat"
% Run this code segment under "/Users/hong/Documents/Projects/Glycomics/matlab/". Otherwise some functions of kca won't work.
% inds = [2, 10, 15, 16];
ions = [205.1, 245.1, 347.1, 365.1, 367.1, 407.1, 467.2, 185.1, 509.2];
num = length(ions);
mobility = IsoMal_Mal_95ms.mobility;
colors = 'bmrky';
for k = 1 : num
    ionIdx = find( ions(k) == IsoMal_Mal_95ms.ions );
    ionname = num2str( ions(k) );

    ax = subplot(3, 3, k);
    plot( mobility, IsoMal_Mal_95ms.raw_signals(:, ionIdx), 'k', 'LineWidth', 2 );
    hold on;
    [weights, estimation] = kca.factorize_spectra( IsoMal_Mal_95ms.raw_signals(:, ionIdx)' );
    plot( mobility, estimation, 'k-.', 'LineWidth', 4 );
    for m = 1 : length( weights )
        plot( mobility, kca.mSources.profile(m, :) * weights(m), colors(m), 'LineWidth', 1.5 );
    end
    hold off;
    % title( ['IsoMal & Mal (95ms): ion ', ionname] );
    title( ['Ion ', ionname] );
    ax.TitleFontSizeMultiplier = 1.1;
end

legend( {'Original', 'Approximation', 'Kernel 1', 'Kernel 2', 'Kernel 3', 'Kernel 4'} )

%% Visualize deconvolution weights
figure;
subplot(4,1,1);
bar( kca.mMasses, kca.mSourceWeights(:, 1), 0.25 );
title("Kernel 1");

subplot(4,1,2);
bar( kca.mMasses, kca.mSourceWeights(:, 2), 0.25 );
title("Kernel 2");

subplot(4,1,3);
bar( kca.mMasses, kca.mSourceWeights(:, 3), 0.25 );
title("Kernel 3")

subplot(4,1,4);
bar( kca.mMasses, kca.mSourceWeights(:, 4), 0.25 );
title("Kernel 4")

