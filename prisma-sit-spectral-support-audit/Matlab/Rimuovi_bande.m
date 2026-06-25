%% SCRIPT MATLAB: RIMOZIONE BANDE CON ASSORBIMENTO ATMOSFERICO DA HYPERCUBE PRISMA
%
% Questo script:
%   1. Carica l'hypercube PRISMA
%   2. Identifica e rimuove bande con forte assorbimento atmosferico
%   3. Salva l'hypercube pulito
%   4. Mostra un confronto prima/dopo

clear; close all; clc;

%% 1. CONFIGURAZIONE INIZIALE
% MODIFICA QUESTO PERCORSO CON IL TUO FILE .mat
mat_file = 'path/GeoTIFF_L2D_xxxx_xx_xx/hypercube_prisma_ordinato.mat';
output_folder = 'path/GeoTIFF_L2D_xxxx_xx_xx/output/';

% Crea cartella di output se non esiste
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% 2. CARICAMENTO DELL'HYPERCUBE
fprintf('📂 Caricamento hypercube da: %s\n', mat_file);
load(mat_file, 'hypercube', 'wavelengths_sorted', 'R', 'n_rows', 'n_cols','lat','lon');

% Se le variabili hanno nomi diversi, adatta qui
if ~exist('wavelengths_sorted', 'var') && exist('wavelengths', 'var')
    wavelengths_sorted = wavelengths;
end

fprintf('✅ Hypercube caricato: [%d, %d, %d]\n', size(hypercube));
fprintf('   Range spettrale originale: %.1f - %.1f nm\n', ...
    min(wavelengths_sorted), max(wavelengths_sorted));
fprintf('   Numero bande originale: %d\n', length(wavelengths_sorted));

%% 3. IDENTIFICAZIONE BANDE CON ASSORBIMENTO ATMOSFERICO
fprintf('\n🔍 Identificazione bande con assorbimento atmosferico...\n');

% Definizione delle regioni di assorbimento atmosferico (in nm)
% Basato su atmosfera terrestre: vapor d'acqua, ossigeno, CO2, ecc.
absorption_regions = [
    1350, 1450;  % Forte assorbimento vapor d'acqua
    % 1790, 1950;  % Forte assorbimento vapor d'acqua + CO2
    % 1350, 1450;  % Forte assorbimento vapor d'acqua 
    1800, 1950;  % Forte assorbimento vapor d'acqua + CO2
    % 1350, 1480;  % Forte assorbimento H2O
    % 1790, 2040;  % Forte assorbimento H2O + CO2
    % 2500, 2600;  % Oltre il range PRISMA (non serve)
];

% Aggiungi anche le bande estreme (< 400 nm e > 2400 nm)
extreme_regions = [
    0, 420;      % UV/blu estremo (basso SNR, forte scattering)
    2400, 3000;  % Oltre il limite SWIR
];

% Unisci tutte le regioni
all_bad_regions = [absorption_regions; extreme_regions];

% Inizializza maschera: true = banda da tenere
bande_da_tenere = true(size(wavelengths_sorted));

% Per ogni regione di assorbimento, segna le bande da rimuovere
for i = 1:size(all_bad_regions, 1)
    low = all_bad_regions(i, 1);
    high = all_bad_regions(i, 2);
    mask = (wavelengths_sorted >= low) & (wavelengths_sorted <= high);
    bande_da_tenere = bande_da_tenere & ~mask;
    if any(mask)
        fprintf('   Rimossa regione [%.0f-%.0f nm]: %d bande\n', ...
            low, high, sum(mask));
    end
end

% Rimuovi anche bande con lunghezza d'onda zero (spente)
zero_bands = (wavelengths_sorted == 0);
bande_da_tenere = bande_da_tenere & ~zero_bands;
if any(zero_bands)
    fprintf('   Rimosse %d bande con lunghezza d''onda zero\n', sum(zero_bands));
end

% Statistiche
n_bande_rimosse = sum(~bande_da_tenere);
n_bande_tenute = sum(bande_da_tenere);

fprintf('\n📊 Riepilogo:\n');
fprintf('   Bande totali: %d\n', length(wavelengths_sorted));
fprintf('   Bande da rimuovere: %d (%.1f%%)\n', n_bande_rimosse, ...
    100 * n_bande_rimosse / length(wavelengths_sorted));
fprintf('   Bande da tenere: %d (%.1f%%)\n', n_bande_tenute, ...
    100 * n_bande_tenute / length(wavelengths_sorted));

%% 4. VISUALIZZAZIONE DELLE BANDE RIMOSSE
figure('Name', 'Bande con assorbimento atmosferico', ...
       'Position', [100, 100, 1400, 600]);

% Plot 1: Spettro medio con regioni evidenziate
subplot(1, 2, 1);
mean_spectrum = squeeze(mean(mean(hypercube, 1), 2));
plot(wavelengths_sorted, mean_spectrum, 'b-', 'LineWidth', 1.5);
hold on;

% Evidenzia regioni rimosse
for i = 1:size(all_bad_regions, 1)
    low = all_bad_regions(i, 1);
    high = all_bad_regions(i, 2);
    x_fill = [low, high, high, low];
    y_fill = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

xlabel('Lunghezza d''onda (nm)', 'FontSize', 12);
ylabel('Riflettanza media', 'FontSize', 12);
title('Spettro medio con regioni di assorbimento', 'FontSize', 14);
xlim([400, 2500]);
grid on;

% Plot 2: Maschera delle bande tenute
subplot(1, 2, 2);
hold on;
for i = 1:length(wavelengths_sorted)
    if bande_da_tenere(i)
        plot(wavelengths_sorted(i), 1, 'g.', 'MarkerSize', 10);
    else
        plot(wavelengths_sorted(i), 1, 'r.', 'MarkerSize', 10);
    end
end
xlabel('Lunghezza d''onda (nm)', 'FontSize', 12);
ylabel('Bande', 'FontSize', 12);
title('Bande mantenute (verde) / rimosse (rosso)', 'FontSize', 14);
xlim([400, 2500]);
ylim([0.5, 1.5]);
grid on;
legend({'Bande mantenute', 'Bande rimosse'}, 'Location', 'best');

% Salva la figura
saveas(gcf, fullfile(output_folder, 'bande_rimosse.png'));

%% 5. APPLICAZIONE DELLA MASCHERA ALL'HYPERCUBE
fprintf('\n✂️ Applicazione maschera all''hypercube...\n');

hypercube_clean = hypercube(:, :, bande_da_tenere);
wavelengths_clean = wavelengths_sorted(bande_da_tenere);

fprintf('   Nuovo hypercube: [%d, %d, %d]\n', size(hypercube_clean));
fprintf('   Nuovo range spettrale: %.1f - %.1f nm\n', ...
    min(wavelengths_clean), max(wavelengths_clean));

%% 6. VERIFICA QUALITÀ (OPZIONALE)
fprintf('\n🔍 Verifica qualità bande mantenute...\n');

% Calcola SNR approssimativo per ogni banda mantenuta
% (rapporto tra media e deviazione standard in area omogenea)
snr_bands = zeros(length(wavelengths_clean), 1);
for b = 1:length(wavelengths_clean)
    banda = double(hypercube_clean(:, :, b));
    banda_valid = banda(banda > 0);
    if ~isempty(banda_valid)
        snr_bands(b) = mean(banda_valid) / std(banda_valid);
    end
end

fprintf('   SNR medio bande mantenute: %.2f\n', mean(snr_bands));

% Identifica eventuali bande con SNR molto basso (< 5)
low_snr = snr_bands < 5;
if any(low_snr)
    fprintf('   ⚠️ %d bande mantenute hanno SNR < 5:\n', sum(low_snr));
    for i = find(low_snr)'
        fprintf('      %.1f nm: SNR = %.2f\n', wavelengths_clean(i), snr_bands(i));
    end
end

%% 7. VISUALIZZAZIONE RGB PRIMA/DOPO
fprintf('\n🎨 Confronto RGB prima/dopo rimozione...\n');

% Trova bande RGB prima
[~, idx_r_old] = min(abs(wavelengths_sorted - 670));
[~, idx_g_old] = min(abs(wavelengths_sorted - 550));
[~, idx_b_old] = min(abs(wavelengths_sorted - 460));

% Trova bande RGB dopo
[~, idx_r_new] = min(abs(wavelengths_clean - 670));
[~, idx_g_new] = min(abs(wavelengths_clean - 550));
[~, idx_b_new] = min(abs(wavelengths_clean - 460));

% Funzione normalizzazione
normalize_band = @(band) normalize_band_function(band);

% RGB originale (con tutte le bande, ma usa solo VNIR)
r_old = double(hypercube(:, :, idx_r_old));
g_old = double(hypercube(:, :, idx_g_old));
b_old = double(hypercube(:, :, idx_b_old));
rgb_old = cat(3, normalize_band(r_old), normalize_band(g_old), normalize_band(b_old));

% RGB pulito
r_new = double(hypercube_clean(:, :, idx_r_new));
g_new = double(hypercube_clean(:, :, idx_g_new));
b_new = double(hypercube_clean(:, :, idx_b_new));
rgb_new = cat(3, normalize_band(r_new), normalize_band(g_new), normalize_band(b_new));

% Mostra confronto
figure('Name', 'Confronto prima/dopo rimozione bande', ...
       'Position', [100, 100, 1600, 600]);

subplot(1, 2, 1);
imshow(rgb_old);
title('RGB originale (con bande assorbimento)', 'FontSize', 14);
axis off;

subplot(1, 2, 2);
imshow(rgb_new);
title('RGB dopo rimozione bande', 'FontSize', 14);
axis off;

% Salva confronto
saveas(gcf, fullfile(output_folder, 'confronto_rgb.png'));

%% 8. SALVATAGGIO DELL'HYPERCUBE PULITO CON STRUTTURA COMPLETA
fprintf('\n💾 Salvataggio hypercube pulito con struttura completa...\n');

% Prepara tutte le variabili necessarie per ricaricare l'hypercube
% Mantieni la stessa struttura del file originale
n_rows_clean = size(hypercube_clean, 1);
n_cols_clean = size(hypercube_clean, 2);
n_bands_clean = size(hypercube_clean, 3);

% Salva come file .mat con tutte le variabili
mat_clean = fullfile(output_folder, 'hypercube_prisma_clean.mat');
save(mat_clean, ...
    'hypercube_clean', ...        % I dati puliti
    'wavelengths_clean', ...      % Le lunghezze d'onda pulite
    'R', ...                       % La georeferenziazione (dal file originale)
    'n_rows_clean', ...            % Numero di righe (esplicito)
    'n_cols_clean', ...            % Numero di colonne (esplicito)
    'n_bands_clean', ...           % Numero di bande (esplicito)
    'bande_da_tenere', ...         % Maschera delle bande mantenute
    'wavelengths_sorted', ...      % Lunghezze d'onda originali (per confronto)
    'lat', ...
    'lon', ...
    '-v7.3');                       % Supporto per file grandi

fprintf('   ✅ Hypercube pulito salvato in: %s\n', mat_clean);
fprintf('   Dimensioni salvate: %d righe, %d colonne, %d bande\n', ...
    n_rows_clean, n_cols_clean, n_bands_clean);

% Salva anche le lunghezze d'onda in CSV (utile per altri software)
csv_clean = fullfile(output_folder, 'wavelengths_clean.csv');
wl_table = table((1:n_bands_clean)', wavelengths_clean, ...
    'VariableNames', {'Band', 'Wavelength_nm'});
writetable(wl_table, csv_clean);
fprintf('   ✅ Lunghezze d''onda salvate in: %s\n', csv_clean);

% Opzionale: salva anche la maschera delle bande per riferimento futuro
mask_file = fullfile(output_folder, 'bande_rimosse_mask.mat');
save(mask_file, 'bande_da_tenere', 'wavelengths_sorted');
fprintf('   ✅ Maschera bande salvata in: %s\n', mask_file);

% % % Opzionale: salva anche in formato GeoTIFF se necessario
% % try
% %     tif_clean = fullfile(output_folder, 'hypercube_prisma_clean.tif');
% % 
% %     % Nota: geotiffwrite potrebbe non supportare molte bande
% %     if n_bands_clean <= 100  % Solo se non sono troppe bande
% %         geotiffwrite(tif_clean, hypercube_clean, R);
% %         fprintf('   ✅ Hypercube salvato anche come GeoTIFF (prime %d bande)\n', ...
% %             min(n_bands_clean, 100));
% %     else
% %         fprintf('   ⚠️ Troppe bande per GeoTIFF (%d), salvo solo le prime 100\n', n_bands_clean);
% %         geotiffwrite(tif_clean, hypercube_clean(:, :, 1:100), R);
% %     end
% % catch ME
% %     fprintf('   ⚠️ Impossibile salvare GeoTIFF: %s\n', ME.message);
% % end
%% 9. STATISTICHE RIASSUNTIVE
% fprintf('\n' + repmat('=', 1, 60) + '\n');
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('             PULIZIA COMPLETATA\n');
%fprintf(repmat('=', 1, 60) + '\n');
fprintf('%s\n', repmat('=', 1, 60)); 
fprintf('📊 RIEPILOGO FINALE:\n');
fprintf('   Bande originali: %d\n', size(hypercube, 3));
fprintf('   Bande rimosse: %d\n', n_bande_rimosse);
fprintf('   Bande mantenute: %d\n', size(hypercube_clean, 3));
fprintf('   Range originale: %.1f - %.1f nm\n', ...
    min(wavelengths_sorted), max(wavelengths_sorted));
fprintf('   Range pulito: %.1f - %.1f nm\n', ...
    min(wavelengths_clean), max(wavelengths_clean));
fprintf('   File salvati in: %s\n', output_folder);
% fprintf(repmat('=', 1, 60) + '\n');
fprintf('%s\n', repmat('=', 1, 60));

%% FUNZIONE DI SUPPORTO
function band_norm = normalize_band_function(band)
    % Normalizza una banda per visualizzazione usando percentili
    band = double(band);
    pixel_validi = band(band > 0);
    
    if isempty(pixel_validi)
        band_norm = zeros(size(band));
        return;
    end
    
    vmin = prctile(pixel_validi, 2);
    vmax = prctile(pixel_validi, 98);
    
    if vmax > vmin
        band_norm = (band - vmin) / (vmax - vmin);
    else
        band_norm = zeros(size(band));
    end
    
    band_norm = max(0, min(1, band_norm));
end