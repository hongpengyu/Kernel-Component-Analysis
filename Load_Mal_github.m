%% Load data
datapath = './data/';

%% Load IsoMal 95ms
IsoMal_95ms = [];
IsoMal_95ms.name = 'IsoMal (95 ms)';

filename = [datapath, 'all_atds_isomal3_95ms_v3.csv'];

fid = fopen( filename );
header = fgetl(fid);
fclose( fid );

header = strrep( header, 'abund', '');
temp = strsplit( header, ',' );
temp = temp(4:end);

IsoMal_95ms.ions = zeros(1, length(temp));
for k = 1 : length(temp)
    IsoMal_95ms.ions(k) = str2num(temp{k});
end

temp = readmatrix( filename, 'NumHeaderLines', 1 );
IsoMal_95ms.mobility = temp(:, 3);
IsoMal_95ms.raw_signals = temp(:, 4:end);


%% Load Mal 95ms
Mal_95ms = [];
Mal_95ms.name = 'Mal (95 ms)';

filename = [datapath, 'all_atds_mal3_95ms_v3.csv'];

fid = fopen( filename );
header = fgetl(fid);
fclose( fid );

header = strrep( header, 'abund', '');
temp = strsplit( header, ',' );
temp = temp(4:end);

Mal_95ms.ions = zeros(1, length(temp));
for k = 1 : length(temp)
    Mal_95ms.ions(k) = str2num(temp{k});
end

temp = readmatrix( filename, 'NumHeaderLines', 1 );
Mal_95ms.mobility = temp(:, 3);
Mal_95ms.raw_signals = temp(:, 4:end);


%% Load IsoMal & Mal 95ms
IsoMal_Mal_95ms = [];
IsoMal_Mal_95ms.name = 'IsoMal & Mal (95 ms)';

filename = [datapath, 'all_atds_isomal_mal3_95ms_v3.csv'];

fid = fopen( filename );
header = fgetl(fid);
fclose( fid );

header = strrep( header, 'abund', '');
temp = strsplit( header, ',' );
temp = temp(4:end);

IsoMal_Mal_95ms.ions = zeros(1, length(temp));
for k = 1 : length(temp)
    IsoMal_Mal_95ms.ions(k) = str2num(temp{k});
end

temp = readmatrix( filename, 'NumHeaderLines', 1 );
IsoMal_Mal_95ms.mobility = temp(:, 3);
IsoMal_Mal_95ms.raw_signals = temp(:, 4:end);


%% Remove zero signals -- 95 ms
temp = [IsoMal_95ms.raw_signals, Mal_95ms.raw_signals, IsoMal_Mal_95ms.raw_signals];
[~, sIdx, eIdx] = remove_zero_signals( temp );

IsoMal_95ms.mobility = IsoMal_95ms.mobility(sIdx:eIdx);
IsoMal_95ms.raw_signals = IsoMal_95ms.raw_signals(sIdx:eIdx, :);

Mal_95ms.mobility = Mal_95ms.mobility(sIdx:eIdx);
Mal_95ms.raw_signals = Mal_95ms.raw_signals(sIdx:eIdx, :);

IsoMal_Mal_95ms.mobility = IsoMal_Mal_95ms.mobility(sIdx:eIdx);
IsoMal_Mal_95ms.raw_signals = IsoMal_Mal_95ms.raw_signals(sIdx:eIdx, :);
