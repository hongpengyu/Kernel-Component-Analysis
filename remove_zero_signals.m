function [results, sIdx, eIdx] = remove_zero_signals( signals )
% Each column of "signals" is a sample

profile = mean( signals, 2 );
N = length(profile);
mask = ones(N, 1);

for k = 1 : N-1
    if profile(k) > 1
        break;
    else
        mask(k) = 0;
    end
end

for k = N : -1 : 2
    if profile(k) > 1
        break;
    else
        mask(k) = 0;
    end
end

sIdx = find( mask, 1, 'first' );
eIdx = find( mask, 1, 'last' );

results = signals(sIdx:eIdx, :);