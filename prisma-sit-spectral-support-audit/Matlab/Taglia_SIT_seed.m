%% SCRIPT MATLAB DEFINITIVO: INTEGRAZIONE IPERCUBO PRISMA + USO SUOLO PUGLIA
% Versione finale con:
% - Caricamento ipercubo e rigenerazione coordinate UTM da R
% - Lettura CSV con poligoni (generato da GDAL)
% - Estrazione WKT, filtraggio spaziale, rasterizzazione
% - ESCLUSIONE DEL MARE tramite soglia NIR
% - Preparazione dati per contrastive learning
% 
% Data: 4 Giugno 2026

% attenzione modificare scane_id

clear; close all; clc;

%% 1. CONFIGURAZIONE INIZIALE
hypercube_file = '/Users/andreaguerriero/Documents/Matlab/PRISMA/GeoTIFF_L2D_2022_03_14/output/hypercube_prisma_clean.mat';
csv_path = '/Users/andreaguerriero/Documents/Matlab/PRISMA/SIT/Shape/Completo/puglia_completo_UTM.csv';
% output_folder = '/Users/andreaguerriero/Documents/Matlab/PRISMA/SIT/Shape/Completo/Tagli/';
output_folder = '/Users/andreaguerriero/Documents/Matlab/PRISMA/GeoTIFF_L2D_2022_03_14/output/';

livello_da_usare = 'LIVELLO_4';
soglia_area_minima = 0.5;  % Ettari

scene_id = "Scena_2";

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

random_seed = 49;

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf(' INTEGRAZIONE IPERCUBO PRISMA + USO SUOLO PUGLIA\n');
fprintf('%s\n', repmat('=', 1, 70));

%% 2. CARICA IPERCUBO E RIGENERA COORDINATE UTM
fprintf('\n📂 Caricamento ipercubo PRISMA...\n');
%load(hypercube_file, 'hypercube_clean', 'wavelengths_sorted', 'R');
load(hypercube_file, "hypercube_clean", "wavelengths_clean", "R");

wl = wavelengths_clean;
wavelengths_sorted = wavelengths_clean;

fprintf("Bands in cube: %d\n", size(hypercube_clean, 3));
fprintf("Wavelengths:   %d\n", numel(wl));


n_rows = size(hypercube_clean, 1);
n_cols = size(hypercube_clean, 2);
n_bands = size(hypercube_clean, 3);

fprintf('   Dimensioni: %d x %d x %d\n', n_rows, n_cols, n_bands);

%% 2b. RIGENERAZIONE COORDINATE E NORMALIZZAZIONE CRS
fprintf('\n🔄 Rigenerazione coordinate UTM e normalizzazione CRS...\n');

if isa(R, 'map.rasterref.MapCellsReference')
    [x_intrinsic, y_intrinsic] = meshgrid(1:n_cols, 1:n_rows);
    [x_world, y_world] = intrinsicToWorld(R, x_intrinsic, y_intrinsic);
else
    error('Tipo di R non supportato: %s', class(R));
end

% CRS del SIT/poligoni CSV.
sit_epsg = 32633;

% -------------------------------------------------------------------------
% Stima automatica del CRS PRISMA dai valori di Easting.
%
% Regola pratica per la Puglia:
% - EPSG:32633 -> Easting tipicamente ~600000-800000
% - EPSG:32634 -> Easting tipicamente ~200000-350000
%
% Se Easting è basso, assumiamo UTM 34N.
% Se Easting è alto, assumiamo UTM 33N.
% -------------------------------------------------------------------------

x_min_raw = min(x_world(:));
x_max_raw = max(x_world(:));
x_med_raw = median(x_world(:), 'omitnan');

y_min_raw = min(y_world(:));
y_max_raw = max(y_world(:));

fprintf('Coordinate raw da R:\n');
fprintf('   Easting raw:  [%.1f, %.1f], median %.1f\n', x_min_raw, x_max_raw, x_med_raw);
fprintf('   Northing raw: [%.1f, %.1f]\n', y_min_raw, y_max_raw);

if x_med_raw < 400000
    prisma_epsg = 32634;
elseif x_med_raw > 500000
    prisma_epsg = 32633;
else
    warning('Easting mediano ambiguo: %.1f. Uso EPSG:32633 come default.', x_med_raw);
    prisma_epsg = 32633;
end

fprintf('   CRS PRISMA stimato: EPSG:%d\n', prisma_epsg);
fprintf('   CRS SIT target:     EPSG:%d\n', sit_epsg);

% -------------------------------------------------------------------------
% Se la scena è in UTM 34N, la riproiettiamo in UTM 33N.
% -------------------------------------------------------------------------

if prisma_epsg ~= sit_epsg
    fprintf('   Riproiezione coordinate PRISMA da EPSG:%d a EPSG:%d...\n', ...
        prisma_epsg, sit_epsg);

    crs_prisma = projcrs(prisma_epsg);
    crs_sit = projcrs(sit_epsg);

    % UTM PRISMA -> lat/lon
    [lat_geo, lon_geo] = projinv(crs_prisma, x_world, y_world);

    % lat/lon -> UTM SIT
    [x_sit, y_sit] = projfwd(crs_sit, lat_geo, lon_geo);

    lon = x_sit;
    lat = y_sit;
else
    lon = x_world;
    lat = y_world;
end

fprintf('✅ Coordinate finali nel CRS SIT EPSG:%d:\n', sit_epsg);
fprintf('   Easting:  [%.1f, %.1f]\n', min(lon(:)), max(lon(:)));
fprintf('   Northing: [%.1f, %.1f]\n', min(lat(:)), max(lat(:)));

lat_min = min(lat(:));
lat_max = max(lat(:));
lon_min = min(lon(:));
lon_max = max(lon(:));
%% 3. CREA MASCHERA MARE/TERRA USANDO BANDA NIR
%% CREAZIONE MASCHERA MARE ROBUSTA


% Trova bande per NDVI (Red ≈ 670nm, NIR ≈ 800nm)
[~, idx_red] = min(abs(wavelengths_clean - 670));
[~, idx_nir] = min(abs(wavelengths_clean - 800));

red_band = double(hypercube_clean(:, :, idx_red));
nir_band = double(hypercube_clean(:, :, idx_nir));

% Calcola NDVI
ndvi = (nir_band - red_band) ./ (nir_band + red_band + eps);
% Soglia per il mare (valori tipici NDVI per acqua < -0.1)
soglia_mare_ndvi = -0.01;
fprintf('   Soglia mare NDVI: %.2f\n', soglia_mare_ndvi);

% Maschera mare basata su NDVI
mask_mare = ndvi < soglia_mare_ndvi;


percentuale_mare = 100 * sum(mask_mare(:)) / numel(mask_mare);
fprintf('   Pixel classificati come mare: %.1f%%\n', percentuale_mare);
%% 4. LEGGI CSV CON I POLIGONI
fprintf('\n📦 Lettura file CSV...\n');
T = readtable(csv_path);
fprintf('   Caricati %d poligoni\n', height(T));


poligoni_file = '/Users/andreaguerriero/Documents/Matlab/PRISMA/SIT/Shape/Completo/poligoni_estratti.mat';
load(poligoni_file, 'S');

%% 6. FILTRAGGIO SPAZIALE
fprintf('\n🔍 Filtraggio poligoni nell''area PRISMA...\n');

S_filtered = struct([]);
filt_count = 0;

for i = 1:length(S)
    if mod(i, 10000) == 0
        fprintf('   Filtrati %d/%d...\n', i, length(S));
    end
    
    xmin = min(S(i).X); xmax = max(S(i).X);
    ymin = min(S(i).Y); ymax = max(S(i).Y);
    
    if xmax < lon_min || xmin > lon_max || ymax < lat_min || ymin > lat_max
        continue;
    end
    
    area_ha = 1;
    if isfield(S(i), 'SHAPE_AREA')
        area_ha = S(i).SHAPE_AREA / 10000;
    end
    if area_ha < soglia_area_minima, continue; end
    
    filt_count = filt_count + 1;
    if isempty(S_filtered)
        S_filtered = S(i);
    else
        S_filtered(filt_count) = S(i);
    end
end

fprintf('✅ Trovati %d poligoni nell''area\n', filt_count);
fprintf('\n--- DEBUG FILTRAGGIO SPAZIALE ---\n');
fprintf('Scene bounds:\n');
fprintf('  lon/easting:  %.2f - %.2f\n', lon_min, lon_max);
fprintf('  lat/northing: %.2f - %.2f\n', lat_min, lat_max);
fprintf('S_filtered length: %d\n', length(S_filtered));

if isempty(S_filtered)
    warning('S_filtered è vuoto: nessun poligono SIT interseca questa scena.');
else
    fprintf('Primi poligoni filtrati:\n');
    for kk = 1:min(5, length(S_filtered))
        fprintf('  kk=%d, LIVELLO_4=%d, area=%.2f ha, bbox=[%.1f %.1f %.1f %.1f]\n', ...
            kk, ...
            double(S_filtered(kk).LIVELLO_4), ...
            double(S_filtered(kk).SHAPE_AREA)/10000, ...
            min(S_filtered(kk).X), max(S_filtered(kk).X), ...
            min(S_filtered(kk).Y), max(S_filtered(kk).Y));
    end
end
fprintf('---------------------------------\n');

%% 7. RASTERIZZAZIONE CON ESCLUSIONE DEL MARE
fprintf('\n🎨 Rasterizzazione uso suolo (escludendo il mare)...\n');

mask = zeros(n_rows, n_cols, 'uint16');
id_map = zeros(n_rows, n_cols, 'uint32');
desc_map = strings(n_rows, n_cols);

lon_vec = lon(:);
lat_vec = lat(:);
[x_idx, y_idx] = meshgrid(1:n_cols, 1:n_rows);
x_idx = x_idx(:); y_idx = y_idx(:);

for i = 1:length(S_filtered)
    if mod(i, 500) == 0
        fprintf('   Rasterizzati %d/%d poligoni...\n', i, length(S_filtered));
    end
    
    px = S_filtered(i).X;
    py = S_filtered(i).Y;
    
    in = inpolygon(lon_vec, lat_vec, px, py);
    idx = find(in);
    
    codice = uint16(S_filtered(i).(livello_da_usare));
    desc = S_filtered(i).DESC_;
    if iscell(desc), desc = desc{1}; end
    
    for j = 1:length(idx)
        r = y_idx(idx(j));
        c = x_idx(idx(j));
        % Assegna solo se non già assegnato E non è mare
        if mask(r,c) == 0 && ~mask_mare(r,c)
            mask(r,c) = codice;
            id_map(r,c) = i;
            desc_map(r,c) = desc;
        end
    end
end


fprintf('\n--- DEBUG RASTERIZZAZIONE ---\n');
fprintf('Pixel mask > 0 prima del controllo finale: %d\n', sum(mask(:) > 0));
fprintf('Pixel mare: %d / %d (%.2f%%)\n', ...
    sum(mask_mare(:)), numel(mask_mare), 100*sum(mask_mare(:))/numel(mask_mare));

if ~isempty(S_filtered)
    test_hits = 0;

    lon_vec = lon(:);
    lat_vec = lat(:);

    for kk = 1:min(20, length(S_filtered))
        in_test = inpolygon(lon_vec, lat_vec, S_filtered(kk).X, S_filtered(kk).Y);
        n_in = sum(in_test);
        n_in_land = sum(in_test & ~mask_mare(:));

        fprintf('  poly %d: pixels inside=%d, inside non-mare=%d, LIVELLO_4=%d\n', ...
            kk, n_in, n_in_land, double(S_filtered(kk).LIVELLO_4));

        test_hits = test_hits + n_in_land;
    end

    fprintf('Totale non-mare nei primi poligoni testati: %d\n', test_hits);
end
fprintf('-----------------------------\n');


pix_valid = sum(mask(:) > 0);
fprintf('✅ Pixel con uso suolo: %d (%.1f%% dell''area totale)\n', ...
    pix_valid, 100 * pix_valid / numel(mask));

%% 8. TABELLA CLASSI ORIGINALI SIT
fprintf('\n📋 Tabella classi uso suolo originali SIT:\n');

classi = unique(mask(mask > 0));

tab = table();
tab.Codice = classi;
tab.Descrizione = strings(length(classi), 1);
tab.NumPixel = zeros(length(classi), 1);
tab.Area_ha = zeros(length(classi), 1);
tab.Percentuale = zeros(length(classi), 1);

area_totale_ha = sum(mask(:) > 0) * 900 / 10000;

for i = 1:length(classi)
    c = classi(i);
    mask_c = (mask == c);

    tab.NumPixel(i) = sum(mask_c(:));
    tab.Area_ha(i) = tab.NumPixel(i) * 900 / 10000;
    tab.Percentuale(i) = 100 * tab.Area_ha(i) / area_totale_ha;

    [r0, c0] = find(mask_c, 1);

    if ~isempty(r0)
        tab.Descrizione(i) = desc_map(r0, c0);
    else
        tab.Descrizione(i) = getClassDescription(c);
    end

    fprintf('   %5d: %-40s %8d pixel %10.1f ha (%5.2f%%)\n', ...
        c, tab.Descrizione(i), tab.NumPixel(i), tab.Area_ha(i), tab.Percentuale(i));
end

writetable(tab, fullfile(output_folder, 'classi_uso_suolo.csv'));

fprintf('\n✅ Nessuna unificazione applicata: uso dei codici SIT originali LIVELLO_4.\n');
%% 9. VISUALIZZAZIONE CON E SENZA MARE
fprintf('\n🎨 Creazione visualizzazioni...\n');

figure('Name', 'Uso suolo con/senza mare', 'Position', [100, 100, 1800, 800]);

subplot(1,3,1);
imagesc(lon(1,:), lat(:,1), mask);
title('Uso suolo (senza mare)');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar;
axis xy equal tight;

subplot(1,3,2);
imagesc(lon(1,:), lat(:,1), mask_mare);
title('Maschera mare');
xlabel('Easting (m)');
ylabel('Northing (m)');
colormap(gca, 'gray');
axis xy equal tight;

subplot(1,3,3);
[~, idx_r] = min(abs(wavelengths_clean - 670));
[~, idx_g] = min(abs(wavelengths_clean - 550));
[~, idx_b] = min(abs(wavelengths_clean - 460));

r_band = double(hypercube_clean(:, :, idx_r));
g_band = double(hypercube_clean(:, :, idx_g));
b_band = double(hypercube_clean(:, :, idx_b));

r_norm = r_band / prctile(r_band(r_band>0), 98);
g_norm = g_band / prctile(g_band(g_band>0), 98);
b_norm = b_band / prctile(b_band(b_band>0), 98);
rgb = cat(3, r_norm, g_norm, b_norm);
rgb = min(rgb, 1);
imshow(rgb);
title('PRISMA RGB');
axis on;
xlabel('Easting (m)');
ylabel('Northing (m)');
% Aggiungi i tick delle coordinate
xticks = linspace(1, n_cols, 5);
yticks = linspace(1, n_rows, 5);
set(gca, 'XTick', xticks, 'XTickLabel', arrayfun(@(x) sprintf('%.2f', x), interp1(1:n_cols, lon(1,:), xticks), 'UniformOutput', false));
set(gca, 'YTick', yticks, 'YTickLabel', arrayfun(@(y) sprintf('%.2f', y), interp1(1:n_rows, lat(:,1), yticks), 'UniformOutput', false));

saveas(gcf, fullfile(output_folder, 'diagnostica_mare.png'));


%% 10. VISUALIZZAZIONE COLORI DISTINTI (solo terra) - CON ORIENTAMENTO CORRETTO
figure('Name', 'Uso suolo - colori distinti', 'Position', [100, 100, 1200, 800]);

% Crea mappa indici
classi_uniques = classi;
n_classi = length(classi_uniques);
mask_index = zeros(size(mask));

% Crea un array di descrizioni per ogni classe
descrizioni_classi = strings(n_classi, 1);

for i = 1:n_classi
    codice = classi_uniques(i);
    mask_index(mask == codice) = i;
    
    % Cerca la descrizione nella tabella
    idx_tab = find(tab.Codice == codice, 1);
    if ~isempty(idx_tab)
        descrizioni_classi(i) = tab.Descrizione(idx_tab);
    else
        descrizioni_classi(i) = getClassDescription(codice);
    end
end

% APPLICA LA STESSA ROTAZIONE/USO CHE HAI USATO PER LE ALTRE FIGURE
% (supponendo che tu abbia usato flipud o rot90 in precedenza)

% Verifica l'orientamento attuale
% Se stavi usando flipud nelle altre figure, applicalo anche qui:
mask_index_display = mask_index;  % o flipud(mask_index) se necessario
lon_display = lon(1,:);  % mantieni le coordinate originali
lat_display = lat(:,1);  % mantieni le coordinate originali

% Se necessario, applica la correzione (adatta in base al tuo caso)
% mask_index_display = flipud(mask_index);  % decommenta se serve

% Mostra l'immagine con axis xy per origine in alto a sinistra
imagesc(lon_display, lat_display, mask_index_display);
axis xy;  % FORZA l'orientamento corretto (origine in alto a sinistra)

cmap = colorcube(n_classi);
colormap(gca, cmap);

% Configura la colorbar
c = colorbar;
c.Ticks = 1:n_classi;
c.TickLabels = descrizioni_classi;

% Tronca etichette troppo lunghe
if n_classi > 5
    for i = 1:n_classi
        if strlength(descrizioni_classi(i)) > 20
            descrizioni_classi(i) = extractBefore(descrizioni_classi(i), 20) + "...";
        end
    end
    c.TickLabels = descrizioni_classi;
end

c.Label.String = 'Classe di uso suolo';
c.Label.FontSize = 10;

title('Uso suolo - colori distinti');
xlabel('Easting (m)');
ylabel('Northing (m)');
axis equal tight;  % Mantiene le proporzioni

saveas(gcf, fullfile(output_folder, 'uso_suolo_colori_distinti_corretto.png'));
fprintf('✅ Visualizzazione con nomi classi salvata (orientamento corretto)\n');
%% 11. SALVATAGGIO GEOTIFF
fprintf('\n💾 Salvataggio risultati...\n');

try
    % PRISMA/SIT data are in UTM zone 33N
    epsg_code = 32633;

    if prisma_epsg == sit_epsg
    geotiffwrite(fullfile(output_folder, 'uso_suolo_prisma.tif'), ...
        mask, R, ...
        'CoordRefSysCode', sit_epsg, ...
        'TiffType', 'bigtiff');
else
    warning('GeoTIFF non salvato: la scena è stata riproiettata da EPSG:%d a EPSG:%d ma R non è stato aggiornato.', prisma_epsg, sit_epsg);
end

    fprintf('   ✅ GeoTIFF salvato (EPSG:%d)\n', epsg_code);

catch ME
    fprintf('   ⚠️ Errore GeoTIFF: %s\n', ME.message);
end
%% 12. PREPARAZIONE DATI PER CONTRASTIVE LEARNING (VERSIONE OTTIMIZZATA)
fprintf('\n🧠 Preparazione dati per contrastive learning...\n');

% 1. Crea una maschera dell'immagine PRISMA (pixel con dati validi)
banda_test = double(hypercube_clean(:, :, 1));
mask_prisma = banda_test > 0;
pixel_prisma = sum(mask_prisma(:));
fprintf('   Passo 1: Immagine PRISMA ha %d pixel validi\n', pixel_prisma);

% 2. Trova TUTTI i pixel con uso suolo (valori originali, prima dell'unificazione)
mask_uso_suolo = mask > 0;
pixel_uso = sum(mask_uso_suolo(:));
fprintf('   Passo 2: Uso suolo ha %d pixel\n', pixel_uso);

% 3. Intersezione: pixel validi in entrambe
mask_valida = mask_prisma & mask_uso_suolo;
pixel_validi_tot = sum(mask_valida(:));
fprintf('   Passo 3: Intersezione: %d pixel validi\n', pixel_validi_tot);

if pixel_validi_tot == 0
    error('❌ Nessun pixel valido trovato! Verifica le maschere.');
end

% 4. Estrai coordinate e codici ORIGINALI (versione ottimizzata)
fprintf('   Passo 4: Estrazione coordinate e codici...\n');

% Usa find con output separati per righe e colonne
[pixel_r, pixel_c] = find(mask_valida);
n_pixel_totali = length(pixel_r);
fprintf('      Trovate %d coordinate\n', n_pixel_totali);

% Prealloca il vettore dei codici
codici_originali = zeros(n_pixel_totali, 1, 'uint16');

% Estrai i codici in blocchi per evitare problemi di memoria
blocco_size = 100000;
fprintf('      Estrazione codici in blocchi da %d...\n', blocco_size);

for start_idx = 1:blocco_size:n_pixel_totali
    end_idx = min(start_idx + blocco_size - 1, n_pixel_totali);
    for i = start_idx:end_idx
        codici_originali(i) = mask(pixel_r(i), pixel_c(i));
    end
    fprintf('         Estratti %d/%d codici\n', end_idx, n_pixel_totali);
end

fprintf('   Passo 4 completato: %d pixel estratti\n', n_pixel_totali);

% 5. Statistiche rapide sui codici originali
classi_originali = unique(codici_originali);
fprintf('   Classi originali trovate: %d\n', length(classi_originali));

% Mostra le prime classi per numerosità
conteggi_originali = zeros(length(classi_originali), 1);
for i = 1:length(classi_originali)
    conteggi_originali(i) = sum(codici_originali == classi_originali(i));
end
[~, ordine] = sort(conteggi_originali, 'descend');
fprintf('   Top 5 classi originali:\n');
for i = 1:min(5, length(classi_originali))
    idx = ordine(i);
    fprintf('      Codice %d: %d pixel\n', classi_originali(idx), conteggi_originali(idx));
end
% 6. Nessuna unificazione: le label sono i codici SIT originali
fprintf('\n   Passo 5: Nessuna unificazione codici.\n');
fprintf('   Le label y useranno direttamente i codici SIT originali %s.\n', livello_da_usare);

codici_finali = codici_originali;

% Per tracciabilità esplicita
original_sit_codes_all = codici_originali;
expected_y_all = codici_originali;

fprintf('   Classi originali mantenute: %d\n', length(unique(codici_finali)));



% % Verifica presenza seminativi
% presente_2111 = any(codici_unificati == 2111);
% fprintf('\n   Codice 2111 presente dopo unificazione? %s\n', string(presente_2111));
% if presente_2111
%     n_2111 = sum(codici_unificati == 2111);
%     fprintf('      → %d pixel di seminativi (%.1f%% del totale)\n', ...
%         n_2111, 100 * n_2111 / n_pixel_totali);
% end

% 7. Campiona casualmente
fprintf('\n   Passo 6: Campionamento casuale...\n');

max_pix = 100000;

if n_pixel_totali > max_pix
rng(random_seed, 'twister');
campione_idx = randperm(n_pixel_totali, max_pix);

    r_idx = pixel_r(campione_idx);
    c_idx = pixel_c(campione_idx);

    y = codici_finali(campione_idx);
    original_sit_codes = original_sit_codes_all(campione_idx);
    expected_y = expected_y_all(campione_idx);

    n_valid = max_pix;
else
    r_idx = pixel_r;
    c_idx = pixel_c;

    y = codici_finali;
    original_sit_codes = original_sit_codes_all;
    expected_y = expected_y_all;

    n_valid = n_pixel_totali;
end

fprintf('   Selezionati %d pixel per il training\n', n_valid);


%% SAVE PER-SCENE CAP SUMMARY

cap_summary = table();

cap_summary.scene_id = string(scene_id);
cap_summary.n_valid_pre_cap = n_pixel_totali;
cap_summary.n_exported_post_cap = n_valid;
cap_summary.retained_fraction = n_valid / n_pixel_totali;
cap_summary.cap_applied = n_pixel_totali > max_pix;
cap_summary.max_pix = max_pix;
cap_summary.random_seed = random_seed;

cap_summary.n_classes_pre_cap = numel(unique(codici_finali));
cap_summary.n_classes_post_cap = numel(unique(y));

cap_summary_path = fullfile( ...
    output_folder, ...
    sprintf('%s_seed_%d_cap_summary.csv', scene_id, random_seed) ...
);
writetable(cap_summary, cap_summary_path);

fprintf('✅ Cap summary salvato in:\n%s\n', cap_summary_path);




% Verifica classi nel campione
classi_campione = unique(y);
fprintf('   Classi nel campione: %d\n', length(classi_campione));

% 8. Prealloca e estrai firme spettrali
fprintf('\n   Passo 7: Estrazione firme spettrali...\n');

X = zeros(n_valid, n_bands, 'single');
fid = zeros(n_valid, 1, 'uint16');

for i = 1:n_valid
    X(i,:) = single(hypercube_clean(r_idx(i), c_idx(i), :));
    fid(i) = id_map(r_idx(i), c_idx(i));
    
    if mod(i, 20000) == 0
        fprintf('      Estratti %d/%d pixel\n', i, n_valid);
    end
end


%% 12b. CHECK LIVE COERENZA y / fid / S_filtered SENZA UNIFICAZIONE
fprintf('\n🔎 Controllo LIVE coerenza y/fid con S_filtered in memoria...\n');

expected_y_from_field = nan(length(fid), 1);
original_code_from_field = nan(length(fid), 1);

for ii = 1:length(fid)
    f = double(fid(ii));

    if f < 1 || f > length(S_filtered)
        expected_y_from_field(ii) = NaN;
        original_code_from_field(ii) = NaN;
    else
        original_code = double(S_filtered(f).(livello_da_usare));

        original_code_from_field(ii) = original_code;
        expected_y_from_field(ii) = original_code;
    end
end

y_double = double(y(:));
fid_double = double(fid(:));

mismatch_pixel = y_double ~= expected_y_from_field;
mismatch_pixel(isnan(expected_y_from_field)) = true;

n_mismatch = sum(mismatch_pixel);

fprintf('   Pixel totali: %d\n', length(y_double));
fprintf('   Pixel mismatch y vs original SIT code from field_id: %d (%.4f%%)\n', ...
    n_mismatch, 100 * n_mismatch / length(y_double));

if n_mismatch > 0
    warning('Sono presenti mismatch tra y e S_filtered(field_id). Controllare prima di esportare.');
end

unique_f = unique(fid_double);
check_rows = {};

for kk = 1:length(unique_f)
    f = unique_f(kk);
    idx = fid_double == f;

    yy = y_double(idx);
    ee = expected_y_from_field(idx);

    y_mode = mode(yy);
    e_mode = mode(ee);

    original_code = NaN;
    original_desc = "";

    if f >= 1 && f <= length(S_filtered)
        original_code = double(S_filtered(f).(livello_da_usare));

        if isfield(S_filtered(f), 'DESC_')
            d = S_filtered(f).DESC_;
            if iscell(d), d = d{1}; end
            original_desc = string(d);
        end
    end

    status = "ok";
    if y_mode ~= e_mode
        status = "mismatch";
    end

    check_rows = [check_rows; { ...
        f, sum(idx), y_mode, original_code, e_mode, original_desc, status ...
    }];
end

live_check_table = cell2table(check_rows, ...
    'VariableNames', { ...
        'field_id', ...
        'n_pixels', ...
        'dataset_y_mode', ...
        'original_sit_code', ...
        'expected_original_code', ...
        'original_sit_description', ...
        'status' ...
    });

live_check_path = fullfile(output_folder, ...
    sprintf('LIVE_field_id_label_consistency_check_seed_%d.csv', random_seed));
sfiltered_path = fullfile(output_folder, ...
    sprintf('%s_seed_%d_S_filtered_used_for_dataset.mat', scene_id, random_seed));

writetable(live_check_table, live_check_path);

fprintf('   Tabella LIVE salvata in:\n   %s\n', live_check_path);

live_mismatch_table = live_check_table(live_check_table.status == "mismatch", :);

fprintf('   Field_id mismatch: %d\n', height(live_mismatch_table));

if height(live_mismatch_table) > 0
    fprintf('\n   Primi mismatch:\n');
    disp(live_mismatch_table(1:min(20, height(live_mismatch_table)), :));

    error('STOP: y e field_id non sono coerenti. Dataset non esportabile.');
end




% 9. Salva
fprintf('\n   Passo 8: Salvataggio...\n');
save(fullfile(output_folder, sprintf('dati_contrastive_seed_%d.mat', random_seed)), ...
    'X', ...
    'y', ...
    'fid', ...
    'tab', ...
    'original_sit_codes', ...
    'expected_y', ...
    'r_idx', ...
    'c_idx', ...
    'S_filtered', ...
    '-v7.3');
fprintf('   ✅ Dati salvati: %d campioni, %d feature\n', n_valid, n_bands);

% 10. Statistiche finali sul campione
fprintf('\n📊 Statistiche campione finale:\n');
[classi_finali, ~, ic] = unique(y);
conteggi_finali = histcounts(ic, length(classi_finali));
[conteggi_finali, ordine] = sort(conteggi_finali, 'descend');
classi_finali = classi_finali(ordine);

for i = 1:min(10, length(classi_finali))
    fprintf('   Classe %d: %d pixel (%.1f%%)\n', ...
        classi_finali(i), conteggi_finali(i), 100*conteggi_finali(i)/n_valid);
end

%% 13. VISUALIZZAZIONE DEI PIXEL SELEZIONATI
fprintf('\n🎯 Visualizzazione pixel selezionati...\n');

% Crea una maschera binaria dei pixel selezionati
mask_training = zeros(n_rows, n_cols, 'uint8');
for i = 1:n_valid
    mask_training(r_idx(i), c_idx(i)) = 1;
end

% Crea RGB per visualizzazione
[~, idx_r] = min(abs(wavelengths_sorted - 670));
[~, idx_g] = min(abs(wavelengths_sorted - 550));
[~, idx_b] = min(abs(wavelengths_sorted - 460));

r_band = double(hypercube_clean(:, :, idx_r));
g_band = double(hypercube_clean(:, :, idx_g));
b_band = double(hypercube_clean(:, :, idx_b));

% Normalizzazione robusta
p2_r = prctile(r_band(r_band>0), 2); p98_r = prctile(r_band(r_band>0), 98);
p2_g = prctile(g_band(g_band>0), 2); p98_g = prctile(g_band(g_band>0), 98);
p2_b = prctile(b_band(b_band>0), 2); p98_b = prctile(b_band(b_band>0), 98);

r_norm = max(0, min(1, (r_band - p2_r) / (p98_r - p2_r + eps)));
g_norm = max(0, min(1, (g_band - p2_g) / (p98_g - p2_g + eps)));
b_norm = max(0, min(1, (b_band - p2_b) / (p98_b - p2_b + eps)));

rgb = cat(3, r_norm, g_norm, b_norm);

% APPLICA ROTAZIONE 180 GRADI
rgb = flipud(rgb);  
mask_training = flipud(mask_training);

% Trova i pixel per il plot sulla versione ruotata
[rows_plot, cols_plot] = find(mask_training);


% Trova i pixel per il plot
[rows_plot, cols_plot] = find(mask_training);

% Crea figura
figure('Name', 'Pixel selezionati per training', 'Position', [100, 100, 1600, 700]);

subplot(1,2,1);
imagesc(rgb);
hold on;
plot(cols_plot, rows_plot, 'r.', 'MarkerSize', 1);
hold off;
title(sprintf('RGB + pixel selezionati (%d pixel)', length(rows_plot)));
xlabel('Colonna (pixel)');
ylabel('Riga (pixel)');
axis xy equal tight;

subplot(1,2,2);
imagesc(mask_training);
colormap(gca, 'gray');
title('Maschera pixel selezionati');
xlabel('Colonna (pixel)');
ylabel('Riga (pixel)');
axis xy equal tight;
colorbar;

saveas(gcf, fullfile(output_folder, 'pixel_training_definitivo.png'));
fprintf('✅ Visualizzazione pixel training salvata\n');

%% Verifica bilanciamento classi nel training set
fprintf('\n📊 Verifica bilanciamento classi nel training set:\n');

% Pulisci eventuali variabili che potrebbero interferire
clear xticks yticks;

% Usa i dati REALI del campionamento (y)
classi_campione_uniques = unique(y);
n_classi_campione = length(classi_campione_uniques);

% Calcola i pixel totali per classe NEL CAMPIONE
pixel_per_classe_campione = zeros(n_classi_campione, 1);
for i = 1:n_classi_campione
    pixel_per_classe_campione(i) = sum(y == classi_campione_uniques(i));
end

% Ordina per numero di pixel (decrescente)
[pixel_ordinati, ordine] = sort(pixel_per_classe_campione, 'descend');
classi_ordinate = classi_campione_uniques(ordine);

% CREA UNA MAPPA DI DESCRIZIONI PER TUTTE LE CLASSI
% Unisci le descrizioni dalla tabella tab con quelle mancanti
descrizioni_complete = strings(n_classi_campione, 1);

for i = 1:n_classi_campione
    codice = classi_ordinate(i);
    
    % Cerca nella tabella tab (che ha già le descrizioni unificate)
    idx_tab = find(tab.Codice == codice, 1);
    
    if ~isempty(idx_tab)
        % Se trovato nella tabella, usa quella descrizione
        descrizioni_complete(i) = tab.Descrizione(idx_tab);
    else
        % Se non trovato, usa la funzione di fallback
        descrizioni_complete(i) = getClassDescription(codice);
    end
end

% Crea figura
figure('Name', 'Bilanciamento classi (aggiornato)', 'Position', [100, 100, 1400, 700]);

% Subplot 1: Pixel totali nel campione
subplot(1, 2, 1);
bar(pixel_ordinati);
title(sprintf('Pixel nel campione (%d totali)', n_valid));
xlabel('Classe');
ylabel('Numero pixel');

% Imposta gli xtick con le descrizioni
set(gca, 'XTick', 1:n_classi_campione);
set(gca, 'XTickLabel', descrizioni_complete);
xtickangle(45);
grid on;

% Aggiungi etichette con i valori (solo se non ci sono troppe classi)
if n_classi_campione <= 15
    for i = 1:n_classi_campione
        text(i, pixel_ordinati(i), num2str(pixel_ordinati(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end

% Subplot 2: Percentuali
subplot(1, 2, 2);
percentuali = 100 * pixel_ordinati / n_valid;
bar(percentuali);
title('Percentuale sul campione');
xlabel('Classe');
ylabel('Percentuale (%)');

set(gca, 'XTick', 1:n_classi_campione);
set(gca, 'XTickLabel', descrizioni_complete);
xtickangle(45);
grid on;

if n_classi_campione <= 15
    for i = 1:n_classi_campione
        text(i, percentuali(i), sprintf('%.1f%%', percentuali(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end

saveas(gcf, fullfile(output_folder, 'bilanciamento_classi_aggiornato.png'));
fprintf('✅ Grafico bilanciamento classi aggiornato salvato\n');

% Sovrascrivi il vecchio grafico
saveas(gcf, fullfile(output_folder, 'bilanciamento_classi.png'));




%% ESPORTAZIONE DATI PER PYTHON (CSV)
fprintf('\n🐍 Esportazione dati per Python...\n');

% Crea una cartella per i file Python
python_folder = fullfile(output_folder, sprintf('python_data_seed_%d', random_seed));
if ~exist(python_folder, 'dir')
    mkdir(python_folder);
end
sfiltered_path_python = fullfile(python_folder, sprintf('%s_S_filtered_used_for_dataset.mat', scene_id));
save(sfiltered_path_python, 'S_filtered', '-v7.3');

fprintf('💾 S_filtered salvato anche in python_data:\n%s\n', sfiltered_path_python);

% 1. SALVA LE FIRME SPETTRALI (X)
% X è già una matrice di dimensioni [n_valid, n_bands]
csv_x = fullfile(python_folder, 'X.csv');
writematrix(X, csv_x);
fprintf('   ✅ X.csv salvato (%d x %d)\n', size(X,1), size(X,2));

% 2. SALVA LE ETICHETTE (y)
csv_y = fullfile(python_folder, 'y.csv');
writematrix(y, csv_y);
fprintf('   ✅ y.csv salvato (%d x 1)\n', length(y));

% 3. SALVA GLI ID DEI CAMPI (field_ids)
csv_fid = fullfile(python_folder, 'field_ids.csv');
writematrix(fid, csv_fid);
fprintf('   ✅ field_ids.csv salvato (%d x 1)\n', length(fid));

% 3b. SALVA I CODICI SIT ORIGINALI
csv_original_codes = fullfile(python_folder, 'original_sit_codes.csv');
writematrix(original_sit_codes, csv_original_codes);
fprintf('   ✅ original_sit_codes.csv salvato (%d x 1)\n', length(original_sit_codes));

% 3c. SALVA LE LABEL ATTESE DA field_id
csv_expected_y = fullfile(python_folder, 'expected_y_from_field_id.csv');
writematrix(expected_y, csv_expected_y);
fprintf('   ✅ expected_y_from_field_id.csv salvato (%d x 1)\n', length(expected_y));

% 3d. SALVA ROW/COL DEI PIXEL CAMPIONATI
csv_rows = fullfile(python_folder, 'rows.csv');
csv_cols = fullfile(python_folder, 'cols.csv');

writematrix(r_idx, csv_rows);
writematrix(c_idx, csv_cols);

fprintf('   ✅ rows.csv salvato (%d x 1)\n', length(r_idx));
fprintf('   ✅ cols.csv salvato (%d x 1)\n', length(c_idx));



%% 4. SALVA LE LUNGHEZZE D'ONDA (CORRETTO)
csv_wavelengths = fullfile(python_folder, 'wavelengths.csv');
wavelengths_table = table((1:length(wavelengths_sorted))', wavelengths_sorted, ...
    'VariableNames', {'band_index', 'wavelength_nm'});
writetable(wavelengths_table, csv_wavelengths);
fprintf('   ✅ wavelengths.csv salvato (%d bande)\n', length(wavelengths_sorted));

% 5. SALVA LA TABELLA DELLE CLASSI
csv_classes = fullfile(python_folder, 'classi.csv');
writetable(tab, csv_classes);
fprintf('   ✅ classi.csv salvato (%d classi)\n', height(tab));

% 6. SALVA UN CAMPIONE DELLE COORDINATE (opzionale)
% Utile per visualizzare i dati in Python
r_idx_sample = r_idx(1:min(10000, length(r_idx)));
c_idx_sample = c_idx(1:min(10000, length(c_idx)));
lon_sample = lon(sub2ind(size(lon), r_idx_sample, c_idx_sample));
lat_sample = lat(sub2ind(size(lat), r_idx_sample, c_idx_sample));

coord_table = table(r_idx_sample, c_idx_sample, lon_sample, lat_sample, ...
    'VariableNames', {'row', 'col', 'easting_m', 'northing_m'});
csv_coords = fullfile(python_folder, 'coordinate_campione.csv');
writetable(coord_table, csv_coords);
fprintf('   ✅ coordinate_campione.csv salvato (%d punti)\n', height(coord_table));

% 7. CREA UN FILE README CON LE ISTRUZIONI
readme_file = fullfile(python_folder, 'README.txt');
fid_readme = fopen(readme_file, 'w');
fprintf(fid_readme, 'DATI PRISMA PER PYTHON\n');
fprintf(fid_readme, '======================\n\n');
fprintf(fid_readme, 'File generati il: %s\n\n', datestr(now));
fprintf(fid_readme, 'DESCRIZIONE FILE:\n');
fprintf(fid_readme, '----------------\n');
fprintf(fid_readme, 'X.csv                : matrice [n_campioni x n_bande] - firme spettrali\n');
fprintf(fid_readme, 'y.csv                : vettore [n_campioni x 1] - classe di uso suolo\n');
fprintf(fid_readme, 'field_ids.csv        : vettore [n_campioni x 1] - ID del poligono\n');
fprintf(fid_readme, 'wavelengths.csv      : tabella [n_bande x 2] - lunghezze d''onda (nm)\n');
fprintf(fid_readme, 'classi.csv           : tabella classi SIT originali LIVELLO_4 con pixel, area e percentuale\n');
fprintf(fid_readme, 'original_sit_codes.csv : vettore [n_campioni x 1] - codice SIT originale LIVELLO_4\n');
fprintf(fid_readme, 'expected_y_from_field_id.csv : label attesa ricostruita da S_filtered(field_id)\n');
fprintf(fid_readme, 'rows.csv              : riga pixel PRISMA campionata\n');
fprintf(fid_readme, 'cols.csv              : colonna pixel PRISMA campionata\n');fprintf(fid_readme, 'coordinate_campione.csv: coordinate UTM di un campione di pixel\n\n');
fprintf(fid_readme, 'STATISTICHE:\n');
fprintf(fid_readme, '------------\n');
fprintf(fid_readme, 'n_campioni: %d\n', n_valid);
fprintf(fid_readme, 'n_bande: %d\n', n_bands);
fprintf(fid_readme, 'n_classi: %d\n', height(tab));
fprintf(fid_readme, 'range spettrale: %.1f - %.1f nm\n', min(wavelengths_sorted), max(wavelengths_sorted));
fprintf(fid_readme, 'sistema coordinate: UTM 33N (EPSG:32633)\n\n');
fprintf(fid_readme, 'Per caricare in Python:\n');
fprintf(fid_readme, '----------------------\n');
fprintf(fid_readme, 'import pandas as pd\n');
fprintf(fid_readme, 'import numpy as np\n\n');
fprintf(fid_readme, 'X = pd.read_csv("X.csv", header=None).values\n');
fprintf(fid_readme, 'y = pd.read_csv("y.csv", header=None).values.ravel()\n');
fprintf(fid_readme, 'field_ids = pd.read_csv("field_ids.csv", header=None).values.ravel()\n');
fprintf(fid_readme, 'wavelengths = pd.read_csv("wavelengths.csv")["wavelength_nm"].values\n');
fprintf(fid_readme, 'classi = pd.read_csv("classi.csv")\n');
fprintf(fid_readme, 'original_sit_codes = pd.read_csv("original_sit_codes.csv", header=None).values.ravel()\n');
fprintf(fid_readme, 'expected_y = pd.read_csv("expected_y_from_field_id.csv", header=None).values.ravel()\n');
fprintf(fid_readme, 'rows = pd.read_csv("rows.csv", header=None).values.ravel()\n');
fprintf(fid_readme, 'cols = pd.read_csv("cols.csv", header=None).values.ravel()\n');
fclose(fid_readme);
fprintf('   ✅ README.txt creato\n');

fprintf('\n✅ Dati per Python salvati in: %s\n', python_folder);



%% 13. STATISTICHE FINALI
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('             ELABORAZIONE COMPLETATA\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('📊 RIEPILOGO:\n');
fprintf('   Ipercubo: %d x %d x %d\n', n_rows, n_cols, n_bands);
fprintf('   Poligoni CSV: %d\n', height(T));
fprintf('   Poligoni estratti: %d\n', length(S));
fprintf('   Poligoni in area: %d\n', length(S_filtered));
fprintf('   Pixel mare: %.1f%%\n', percentuale_mare);
fprintf('   Pixel con classe: %d (%.1f%% terra)\n', pix_valid, 100 * pix_valid / (numel(mask) - sum(mask_mare(:))));
fprintf('   Classi distinte: %d\n', length(classi));
fprintf('   Campioni contrastive: %d\n', n_valid);
fprintf('   Output in: %s\n', output_folder);
fprintf('%s\n', repmat('=', 1, 70));

% ============================================================
% FUNZIONE DI SUPPORTO PER DESCRIZIONI MANCANTI
% ============================================================

function desc = getClassDescription(codice)
    % Mappa di fallback per codici che potrebbero mancare nella tabella
    switch codice
        case 2111
            desc = "seminativi semplici in aree non irrigue";
        case 2112
            desc = "colture orticole in aree non irrigue";
        case 2121
            desc = "seminativi semplici in aree irrigue";
        case 2123
            desc = "colture orticole in aree irrigue";
        case 221
            desc = "vigneti";
        case 222
            desc = "frutteti e frutti minori";
        case 223
            desc = "uliveti";
        case 231
            desc = "superfici a copertura erbacea densa";
        case 241
            desc = "colture temporanee associate a colture permanenti";
        case 242
            desc = "sistemi colturali e particellari complessi";
        case 243
            desc = "aree prevalentemente occupate da colture agrarie con presenza di spazi naturali";
        case 311
            desc = "boschi di latifoglie";
        case 312
            desc = "boschi di conifere";
        case 313
            desc = "boschi misti di conifere e latifoglie";
        case 321
            desc = "aree a pascolo naturale e praterie";
        case 322
            desc = "brughiere e cespuglieti";
        case 323
            desc = "aree a vegetazione sclerofilla";
        case 324
            desc = "aree a ricolonizzazione naturale";
        case 3241
            desc = "aree a ricolonizzazione naturale (tipo 1)";
        case 3242
            desc = "aree a ricolonizzazione naturale (tipo 2)";
        case 331
            desc = "spiagge, dune, sabbie";
        case 332
            desc = "rocce nude, falseie, rupi, affioramenti";
        case 333
            desc = "aree con vegetazione rada";
        otherwise
            desc = sprintf("Classe %d", codice);
    end
end
%class_map = codici_originali;

class_map = mask;

save(fullfile(output_folder, "Scena_18_dataset_maps.mat"), ...
    "id_map", ...
    "class_map", ...
    "hypercube_clean", ...
    "wavelengths_sorted", ...
    "-v7.3");