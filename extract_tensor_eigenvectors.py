import numpy as np
import idlsave as idl
import os

# ----------------------------
# Configuration
# ----------------------------
SIMS = ['09_18', '17_11', '37_11']  # Simulations to process
Lgrid = 127.0   # Box size in Mpc
Ngrid = 256     # Grid resolution per dimension
snapshot_file = 'redshift_snap.txt'  # Snapshot-redshift mapping file

tensors = ['Shear', 'Tidal']  # Tensor types to process
smoothings = ['1', '2', '5']  # Smoothing scales

# Host halo IDs per simulation (M31, MW)
host_ids = {
    '09_18': ('127000000000002', '127000000000003'),
    '17_11': ('127000000000002', '127000000000003'),
    '37_11': ('127000000000001', '127000000000002')
}

# ----------------------------
# Load snapshot-redshift table
# ----------------------------
snapshot_l = np.genfromtxt(snapshot_file, dtype=[('snapshot', int), ('z', float)])

# ----------------------------
# Main loop over simulations
# ----------------------------
for SIM in SIMS:
    print(f"Processing simulation: {SIM}")
    id_m31, id_mw = host_ids[SIM]

    # Base path for tensor grids
    pwd = f'/store/clues/HESTIA/RE_SIMS/8192/GAL_FOR/{SIM}/LSS/FIELDS/'
    prefix = f'CIC_8192_GAL_FOR_{SIM}_'

    # Loop over host halos
    for id_halo in [id_m31, id_mw]:
        print(f"  Host halo: {id_halo}")

        # Load satellite infall positions and redshifts
        halo_file = f'/store/clues/HESTIA/RE_SIMS/8192/GAL_FOR/{SIM}/AHF_output_2x2.5Mpc/HESTIA_100Mpc_8192_{SIM}.127_halo_{id_halo}.dat'
        infall_data = np.genfromtxt(halo_file, usecols=(0, 6, 7, 8), names=['zinf', 'xcinf', 'ycinf', 'zcinf'])

        # Convert physical coordinates to Mpc and then to grid indices
        x = infall_data['xcinf'] / 1000
        y = infall_data['ycinf'] / 1000
        z = infall_data['zcinf'] / 1000
        i = (x * Ngrid / Lgrid).astype(int)
        j = (y * Ngrid / Lgrid).astype(int)
        k = (z * Ngrid / Lgrid).astype(int)

        # Loop over tensor types and smoothing scales
        for tensor in tensors:
            print(f"    Processing tensor: {tensor}")
            for smo in smoothings:
                suffix = f'_256_{smo}.0_{tensor}Tensor_Eigen_cell_All'

                # Initialize array to store results
                save = np.zeros((infall_data.shape[0], 13))
                save[:, 0:4] = np.column_stack((infall_data['zinf'], infall_data['xcinf'], infall_data['ycinf'], infall_data['zcinf']))

                # Loop over snapshots
                for snap, z_snap in snapshot_l:
                    sat_snap = infall_data['zinf'] == z_snap
                    if not np.any(sat_snap):
                        continue  # Skip if no satellites at this snapshot

                    # Format snapshot number for file name
                    snap_str = f"{snap:03d}" if snap < 100 else str(snap)
                    tensor_file = os.path.join(pwd, f"{prefix}{snap_str}{suffix}.sav")

                    # Load tensor from IDL file
                    tensor_data = idl.read(tensor_file, verbose=False)
                    if 'tidal' in suffix.lower():
                        eigenvectors = tensor_data.tidaltensor_eigen_cell.eigenvectors[0]
                    else:
                        eigenvectors = tensor_data.sheartensor_eigen_cell.eigenvectors[0]

                    # Extract eigenvectors at satellite positions
                    save[sat_snap, 4:7] = eigenvectors[0, :, k[sat_snap], j[sat_snap], i[sat_snap]]
                    save[sat_snap, 7:10] = eigenvectors[1, :, k[sat_snap], j[sat_snap], i[sat_snap]]
                    save[sat_snap, 10:] = eigenvectors[2, :, k[sat_snap], j[sat_snap], i[sat_snap]]

                # Save results
                out_file = f"{prefix}_{'_'.join(suffix.split('_')[1:5])}_vecs_LG_{id_halo}.npy"
                np.save(out_file, save)
                print(f"      Saved: {out_file}")