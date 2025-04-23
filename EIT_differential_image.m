% === Load EIDORS ===
addpath('C:\Users\loren\OneDrive\Desktop\EIDORS Folder\Senior Project\eidors-v3.12\eidors-v3.12\eidors');
run('C:\Users\loren\OneDrive\Desktop\EIDORS Folder\Senior Project\eidors-v3.12\eidors-v3.12\eidors\startup.m');

% === Load Gmsh mesh ===
mesh_file = 'C:\Users\loren\OneDrive\Desktop\EIDORS Folder\Senior Project\matlab_senior_project\gmsh_files\test_15.msh';
fmdl = gmsh_mk_fwd_model(mesh_file);

% === Electrode layout ===
width = 0.13; depth = 0.13; offset = 0.02;
z_levels = [0.1125, 0.075, 0.0375];
electrode_centers = [];

% FRONT
for z = z_levels
    electrode_centers = [electrode_centers; width/2 - offset, 0, z; width/2 + offset, 0, z];
end
% RIGHT
for z = z_levels
    electrode_centers = [electrode_centers; width, depth/2 - offset, z; width, depth/2 + offset, z];
end
% BACK
for z = z_levels
    electrode_centers = [electrode_centers; width/2 + offset, depth, z; width/2 - offset, depth, z];
end
% LEFT
for z = z_levels
    electrode_centers = [electrode_centers; 0, depth/2 + offset, z; 0, depth/2 - offset, z];
end

fmdl.electrode = struct('nodes', {}, 'z_contact', {});
for k = 1:24
    dists = vecnorm(fmdl.nodes - electrode_centers(k,:), 2, 2);
    [~, closest] = min(dists);
    fmdl.electrode(k).nodes = closest;
    fmdl.electrode(k).z_contact = 0.01;
end

% === Load baseline and inclusion data ===
raw_base = readmatrix('C:\Users\loren\OneDrive\Desktop\EIT DATA\eit_baserun_noinclusion.csv');
raw_incl = readmatrix('C:\Users\loren\OneDrive\Desktop\EIT DATA\eit_inclusionrun2.csv');

Ei_base = round(raw_base(:,1)); Ej_base = round(raw_base(:,2)); Vdiff_base = double(raw_base(:,3));
Ei_incl = round(raw_incl(:,1)); Ej_incl = round(raw_incl(:,2)); Vdiff_incl = double(raw_incl(:,3));

valid_base = ~isnan(Ei_base) & ~isnan(Ej_base) & Ei_base >= 1 & Ej_base <= 24;
valid_incl = ~isnan(Ei_incl) & ~isnan(Ej_incl) & Ei_incl >= 1 & Ej_incl <= 24;

Ei = Ei_base(valid_base); Ej = Ej_base(valid_base);
Vdiff_base = Vdiff_base(valid_base); Vdiff_incl = Vdiff_incl(valid_incl);

if length(Vdiff_base) ~= length(Vdiff_incl)
    error('Baseline and inclusion measurements do not match.');
end

% === Define stimulation (E20 → E10)
n_elec = 24;
stim = struct;
stim.stim_pattern = zeros(n_elec,1);
stim.stim_pattern(20) = +1; stim.stim_pattern(10) = -1;

meas_pattern = zeros(length(Ei), n_elec);
for i = 1:length(Ei)
    meas_pattern(i, Ei(i)) = 1;
    meas_pattern(i, Ej(i)) = -1;
end
stim.meas_pattern = meas_pattern;
fmdl.stimulation = stim;

% === Scale and reshape
scale = 4.096 / 32768;
vh_meas = Vdiff_base * scale;
vi_meas = Vdiff_incl * scale;
vh_meas = reshape(vh_meas, [], 1);
vi_meas = reshape(vi_meas, [], 1);

% === Inverse model setup
inv_model = eidors_obj('inv_model', 'eit_diff_model');
inv_model.fwd_model = fmdl;
inv_model.reconst_type = 'difference';
inv_model.jacobian_bkgnd.value = 1;
inv_model.hyperparameter.value = 1e-2;
inv_model.solve = @inv_solve_gn;
inv_model.RtR_prior = @prior_laplace;

% === Reconstruct image
img_rec = inv_safe(inv_model, vh_meas, vi_meas);
elem_centers = interp_mesh(fmdl);

% === Full field: Resistive + Conductive
img_rec.calc_colours.map = jet(256);
img_rec.calc_colours.clim = min(0.05, max(abs(img_rec.elem_data)));
img_rec.calc_colours.ref_level = 0;
figure;
show_fem_enhanced(img_rec, 1);
title('Full Field: Conductive & Resistive Anomalies');
colorbar;

% === Resistive-only view
resistive_img = img_rec;
resistive_img.elem_data(img_rec.elem_data >= -0.0025) = 0;
resistive_img.calc_colours.map = hot(256);
resistive_img.calc_colours.clim = 0.01;
figure;
show_fem_enhanced(resistive_img, 1);
title('Resistive Anomalies Only (Δσ < -0.0025)');
colorbar;

% === Conductive-only view
conductive_img = img_rec;
conductive_img.elem_data(img_rec.elem_data <= 0.0025) = 0;
conductive_img.calc_colours.map = parula(256);
conductive_img.calc_colours.clim = 0.01;
figure;
show_fem_enhanced(conductive_img, 1);
title('Conductive Anomalies Only (Δσ > +0.0025)');
colorbar;

% === Inclusion (near injection path: E20 → E10)
e20_pos = fmdl.nodes(fmdl.electrode(20).nodes, :);
e10_pos = fmdl.nodes(fmdl.electrode(10).nodes, :);
midline = (e20_pos + e10_pos) / 2;
dists = vecnorm(elem_centers - midline, 2, 2);
anomaly_mask = dists < 0.03 & (img_rec.elem_data > 0.0025 | img_rec.elem_data < -0.0025);

highlight_img = img_rec;
highlight_img.elem_data(~anomaly_mask) = 0;
highlight_img.calc_colours.map = jet(256);
highlight_img.calc_colours.clim = 0.01;

figure;
show_fem_enhanced(highlight_img, 1);
title('Labeled Inclusion Near Injection Path (E20 → E10)');
colorbar;
hold on;

% Label Δσ values on inclusion
idx_inclusion = find(anomaly_mask);
for i = 1:length(idx_inclusion)
    idx = idx_inclusion(i);
    pos = elem_centers(idx,:);
    text(pos(1), pos(2), pos(3)+0.003, sprintf('%.4f', img_rec.elem_data(idx)), ...
        'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'center');
end

% Correct E1–E24 labeling (face-by-face)
electrode_num = 1;
for face = 1:4
    for row = 1:3
        for col = 1:2
            idx = electrode_num;
            pos = fmdl.nodes(fmdl.electrode(idx).nodes, :);
            text(pos(1), pos(2), pos(3)+0.004, ['E', num2str(idx)], ...
                'FontSize', 8, 'Color', 'm', 'HorizontalAlignment', 'center');
            plot3(pos(1), pos(2), pos(3), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
            electrode_num = electrode_num + 1;
        end
    end
end
hold off;

% === Voltage comparison plot
figure;
plot(vh_meas, 'b', 'LineWidth', 1.3); hold on;
plot(vi_meas, 'r', 'LineWidth', 1.3);
legend('Baseline (vh)', 'Inclusion (vi)', 'Location', 'Best');
title('Voltage Comparison: Baseline vs. Inclusion');
xlabel('Measurement Index'); ylabel('Voltage (V)');
grid on;
