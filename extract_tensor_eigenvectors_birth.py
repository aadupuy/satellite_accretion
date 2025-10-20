import numpy as np
import idlsave as idl
import os

# ----------------------------
# Configuration
# ----------------------------
SIMS = ['09_18', '17_11', '37_11']  # Simulations to process
Lgrid = 100.  # Mpc
Ngrid = 256   # grid resolution
tensors = ['Shear', 'Tidal']
smoothing_scales = ['1', '2', '5']
rvirs = ['_0.5rvir', '_1.0rvir', '_1.5rvir', '_2.0rvir']
snapshot_file = 'redshift_snap.txt'
output_dir = 'new'

# ----------------------------
# Load snapshot-redshift table
# ----------------------------
snapshot_l = np.genfromtxt(snapshot_file, dtype=[('snapshot', int), ('z', float)])

# ----------------------------
# Main loop over simulations
# ----------------------------
for SIM in SIMS:
    print(f"Processing simulation: {SIM}")
    prefix = f'CIC_8192_GAL_FOR_{SIM}_'
    field_path = f'/store/clues/HESTIA/RE_SIMS/8192/GAL_FOR/{SIM}/LSS/FIELDS/'

    for tensor in tensors:
        print(f"  Tensor type: {tensor}")

        for smo in smoothing_scales:
            suffix = f'_256_{smo}.0_{tensor}Tensor_Eigen_cell_All'

            for rvir in rvirs:
                print(f"    rvir selection: {rvir}")

                # Load satellite infall/birth data
                infile = f'out_infall_{SIM}{rvir}.txt'
                dtype = [('id', int), ('sat', int), ('zinf', float),
                         ('xcinf','f8'), ('ycinf','f8'), ('zcinf','f8'),
                         ('vxinf','f8'), ('vyinf','f8'), ('vzinf','f8'),
                         ('Mvir_inf','f8'), ('Mgas_inf','f8'), ('Mstar_inf','f8'),
                         ('xc0','f8'), ('yc0','f8'), ('zc0','f8'),
                         ('vx0','f8'), ('vy0','f8'), ('vz0','f8'),
                         ('Mvir0','f8'), ('Mgas0','f8'), ('Mstar0','f8'),
                         ('xclg','f8'), ('yclg','f8'), ('zclg','f8'),
                         ('vxlg','f8'), ('vylg','f8'), ('vzlg','f8'),
                         ('zbirth','f8'), ('xcbirth','f8'), ('ycbirth','f8'), ('zcbirth','f8'),
                         ('rvir','f8')]
                names = [n for n,_ in dtype]
                out = np.genfromtxt(infile, dtype=dtype, names=names)

                save = np.zeros((out.shape[0], 19))
                save[:,0] = out['id']

                # --- Convert coordinates to grid indices ---
                def coords_to_idx(x, y, z):
                    return (x*Ngrid/Lgrid).astype(int), (y*Ngrid/Lgrid).astype(int), (z*Ngrid/Lgrid).astype(int)

                icinf, jcinf, kcinf = coords_to_idx(out['xcinf']/1000, out['ycinf']/1000, out['zcinf']/1000)
                icbirth, jcbirth, kcbirth = coords_to_idx(out['xcbirth']/1000, out['ycbirth']/1000, out['zcbirth']/1000)

                z_inf = out['zinf']
                z_birth = out['zbirth']

                # --- Loop over snapshots ---
                for snap, z_snap in snapshot_l:
                    snapstr = f"{snap:03d}" if snap < 100 else str(snap)

                    # --- Extract at infall ---
                    mask_inf = z_inf == z_snap
                    if np.any(mask_inf):
                        file_path = os.path.join(field_path, f"{prefix}{snapstr}{suffix}.sav")
                        tensor_data = idl.read(file_path, verbose=False)
                        eig = tensor_data.sheartensor_eigen_cell.eigenvectors[0] \
                              if 'shear' in suffix.lower() else tensor_data.tidaltensor_eigen_cell.eigenvectors[0]
                        save[mask_inf,1:4] = eig[0,:,kcinf[mask_inf],jcinf[mask_inf],icinf[mask_inf]]
                        save[mask_inf,4:7] = eig[1,:,kcinf[mask_inf],jcinf[mask_inf],icinf[mask_inf]]
                        save[mask_inf,7:10] = eig[2,:,kcinf[mask_inf],jcinf[mask_inf],icinf[mask_inf]]

                    # --- Extract at birth ---
                    mask_birth = z_birth == z_snap
                    if np.any(mask_birth):
                        file_path = os.path.join(field_path, f"{prefix}{snapstr}{suffix}.sav")
                        tensor_data = idl.read(file_path, verbose=False)
                        eig = tensor_data.sheartensor_eigen_cell.eigenvectors[0] \
                              if 'shear' in suffix.lower() else tensor_data.tidaltensor_eigen_cell.eigenvectors[0]
                        save[mask_birth,10:13] = eig[0,:,kcbirth[mask_birth],jcbirth[mask_birth],icbirth[mask_birth]]
                        save[mask_birth,13:16] = eig[1,:,kcbirth[mask_birth],jcbirth[mask_birth],icbirth[mask_birth]]
                        save[mask_birth,16:] = eig[2,:,kcbirth[mask_birth],jcbirth[mask_birth],icbirth[mask_birth]]

                # --- Save results ---
                outfile = os.path.join(output_dir, f"{prefix}_{'_'.join(suffix.split('_')[1:5])}_vecs_LGs{rvir}.npy")
                np.save(outfile, save)
