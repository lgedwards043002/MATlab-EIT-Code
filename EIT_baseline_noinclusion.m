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

% === Load baseline data only ===
raw = readmatrix('C:\Users\loren\OneDrive\Desktop\EIT DATA\eit_baserun_noinclusion.csv');
Ei = round(raw(:,1)); Ej = round(raw(:,2)); Vdiff_raw = double(raw(:,3));

valid = ~isnan(Ei) & ~isnan(Ej) & Ei >= 1 & Ej >= 1 & Ei <= 24 & Ej <= 24;
Ei = Ei(valid); Ej = Ej(valid); Vdiff_raw = Vdiff_raw(valid);

% === Define stimulation (E20 → E10)
n_elec = 24;
stim = struct;
stim.stim_pattern = zeros(n_elec,1);
stim.stim_pattern(20) = +1;
stim.stim_pattern(10) = -1;

meas_pattern = zeros(length(Ei), n_elec);
for i = 1:length(Ei)
    meas_pattern(i, Ei(i)) = 1;
    meas_pattern(i, Ej(i)) = -1;
end
stim.meas_pattern = meas_pattern;
fmdl.stimulation = stim;

% === Scale Arduino baseline data
scale = 4.096 / 32768;
vi_meas = Vdiff_raw * scale;

% === Simulate homogeneous reference (for comparison)
img_homg = mk_image(fmdl, 1);
vh_struct = fwd_solve(img_homg);
vh_meas = double(vh_struct.meas(:));

% === Absolute inverse model setup
inv_model = eidors_obj('inv_model', 'eit_abs_model');
inv_model.fwd_model = fmdl;
inv_model.reconst_type = 'absolute';
inv_model.jacobian_bkgnd.value = 1;
inv_model.hyperparameter.value = 1e-2;
inv_model.solve = @inv_solve_gn;
inv_model.RtR_prior = @prior_laplace;

% === Reconstruct absolute image
img_rec = inv_solve(inv_model, vi_meas);
img_rec.calc_colours.clim = 5;

% === Display image
figure;
show_fem_enhanced(img_rec, 1);
title('Baseline-Only Absolute Reconstruction (E20 → E10)');
colorbar;

% === Electrode labels
hold on;
for k = 1:24
    pos = fmdl.nodes(fmdl.electrode(k).nodes, :);
    plot3(pos(1), pos(2), pos(3), 'go', 'MarkerSize', 6, 'LineWidth', 1.5);
    text(pos(1), pos(2), pos(3)+0.005, ['E', num2str(k)], ...
        'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'center');
end
hold off;

% === 3D Slice
figure;
show_3d_slices(img_rec, 0.065, 0.065, 0.05);
title('Baseline-Only 3D Slices (E20 → E10)');

% === Voltage plot (blue = simulated, red = measured baseline)
figure;
plot(vh_meas, 'b', 'LineWidth', 1.3); hold on;
plot(vi_meas, 'r', 'LineWidth', 1.3);
legend('Simulated (vh)', 'Measured (vi)', 'Location', 'Best');
title('Voltage Comparison (Baseline Only, E20 → E10)');
xlabel('Measurement Index'); ylabel('Voltage (V)');
grid on;
