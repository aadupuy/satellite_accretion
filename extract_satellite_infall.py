import numpy as np
import os

# -----------------------------
# Parameters
# -----------------------------
rvir_factor = 2.0
list_SIMS = ['09_18','17_11','37_11']

# Column names for input files
LGs_names = [
    'SIM0','SIM1','id_m31','id_mw','id_sim_virgo','M_m31','M_mw','Mratio',
    'Md_star_m31','Md_star_mw','Md_gas_m31','Md_gas_mw','c_gas_m31','c_gas_mw',
    'c_star_m31','c_star_mw','sep','vrad2','vtan','m_sim_virgo','sim_d_virgo',
    'delta_virgo','fhr_m31','fhr_mw','xlg','ylg','zlg','x_m31','y_m31','z_m31',
    'x_mw','y_mw','z_mw','vx_m31','vy_m31','vz_m31','vx_mw','vy_mw'
]
halos_names = [
    'ID','hostHalo','numSubStruct','Mvir','npart','Xc','Yc','Zc','VXc','VYc','VZc',
    'Rvir','Rmax','r2','mbp_offset','com_offset','Vmax','v_esc','sigV','lambda',
    'lambdaE','Lx','Ly','Lz','b','c','Eax','Eay','Eaz','Ebx','Eby','Ebz','Ecx','Ecy','Ecz',
    'ovdens','nbins','fMhires','Ekin','Epot','SurfP','Phi0','cNFW','n_gas','M_gas','lambda_gas',
    'lambdaE_gas','Lx_gas','Ly_gas','Lz_gas','b_gas','c_gas','Eax_gas','Eay_gas','Eaz_gas',
    'Ebx_gas','Eby_gas','Ebz_gas','Ecx_gas','Ecy_gas','Ecz_gas','Ekin_gas','Epot_gas',
    'n_star','M_star','lambda_star','lambdaE_star','Lx_star','Ly_star','Lz_star','b_star','c_star',
    'Eax_star','Eay_star','Eaz_star','Ebx_star','Eby_star','Ebz_star','Ecx_star','Ecy_star','Ecz_star',
    'Ekin_star','Epot_star','mean_z_gas','mean_z_star','n_star_excised','M_star_excised','mean_z_star_excised'
]
tree_names = ['z'] + halos_names  # redshift first, then halo columns

# -----------------------------
# Load Local Group table
# -----------------------------
tab_LGs = np.genfromtxt('LGs_8192_GAL_FOR.txt', names=LGs_names)

# -----------------------------
# Function to process one simulation
# -----------------------------
def process_simulation(SIM):
    print(f"Processing simulation: {SIM}")

    # Load satellites list
    tab_dwarfs = np.genfromtxt(f'dwarfs_{SIM}.txt', dtype=[int,float,int], names=['id_halo','M_halo','label'])
    id_M31 = tab_dwarfs['id_halo'][tab_dwarfs['label']==0][0]
    id_MW = tab_dwarfs['id_halo'][tab_dwarfs['label']==1][0]
    id_M31sat = tab_dwarfs['id_halo'][tab_dwarfs['label']==2]
    id_MWsat = tab_dwarfs['id_halo'][tab_dwarfs['label']==3]

    # Load merger trees
    halos_path = f'/store/clues/HESTIA/RE_SIMS/8192/GAL_FOR/{SIM}/AHF_output_2x2.5Mpc'
    tree_M31 = np.genfromtxt(os.path.join(halos_path, f'HESTIA_100Mpc_8192_{SIM}.127_halo_{id_M31}.dat'), names=tree_names)
    tree_MW = np.genfromtxt(os.path.join(halos_path, f'HESTIA_100Mpc_8192_{SIM}.127_halo_{id_MW}.dat'), names=tree_names)

    all_sats = [(id_M31sat, tree_M31, 1), (id_MWsat, tree_MW, 2)]
    n_sats_total = len(id_M31sat) + len(id_MWsat)
    out = np.zeros((n_sats_total, 35))

    idx = 0
    for sat_ids, lg_tree, lg_label in all_sats:
        for sat_id in sat_ids:
            # Load satellite tree
            tree_sat = np.genfromtxt(os.path.join(halos_path, f'HESTIA_100Mpc_8192_{SIM}.127_halo_{sat_id}.dat'), names=tree_names)
            
            # Skip empty satellites
            if tree_sat.shape == ():
                out[idx, :] = np.nan
                idx += 1
                continue

            # Match lengths
            nzmax = min(lg_tree.shape[0], tree_sat.shape[0])
            tree_lg = lg_tree[:nzmax]
            tree_sat = tree_sat[:nzmax]

            # Distance between LG and satellite at each redshift
            d = np.sqrt((tree_lg['Xc'] - tree_sat['Xc'])**2 +
                        (tree_lg['Yc'] - tree_sat['Yc'])**2 +
                        (tree_lg['Zc'] - tree_sat['Zc'])**2)
            cond = d > rvir_factor * tree_lg['Rvir']

            if np.any(cond) and not np.all(cond):
                first_idx = np.where(cond)[0][0]
                # Extract infall properties
                z_infall = tree_sat['z'][first_idx]
                xc_infall, yc_infall, zc_infall = tree_sat['Xc'][first_idx], tree_sat['Yc'][first_idx], tree_sat['Zc'][first_idx]
                vx_infall, vy_infall, vz_infall = tree_sat['VXc'][first_idx], tree_sat['VYc'][first_idx], tree_sat['VZc'][first_idx]
                Mvir_infall, Mgas_infall, Mstar_infall = tree_sat['Mvir'][first_idx], tree_sat['M_gas'][first_idx], tree_sat['M_star'][first_idx]

                # Current properties at z=0
                z0_idx = np.where(tree_sat['z']==0.)[0][0]
                xc0, yc0, zc0 = tree_sat['Xc'][z0_idx], tree_sat['Yc'][z0_idx], tree_sat['Zc'][z0_idx]
                vx0, vy0, vz0 = tree_sat['VXc'][z0_idx], tree_sat['VYc'][z0_idx], tree_sat['VZc'][z0_idx]
                Mvir0, Mgas0, Mstar0 = tree_sat['Mvir'][z0_idx], tree_sat['M_gas'][z0_idx], tree_sat['M_star'][z0_idx]

                # LG properties at infall
                xclg, yclg, zclg = tree_lg['Xc'][first_idx], tree_lg['Yc'][first_idx], tree_lg['Zc'][first_idx]
                vxlg, vylg, vzlg = tree_lg['VXc'][first_idx], tree_lg['VYc'][first_idx], tree_lg['VZc'][first_idx]
                rvir = tree_lg['Rvir'][first_idx]

                # Birth properties
                zbirth = tree_sat['z'][-1]
                xc_birth, yc_birth, zc_birth = tree_sat['Xc'][-1], tree_sat['Yc'][-1], tree_sat['Zc'][-1]
                xblg, yblg, zblg = tree_lg['Xc'][-1], tree_lg['Yc'][-1], tree_lg['Zc'][-1]
            else:
                # Fill with NaNs if no infall detected
                vals = [np.nan]*35
                out[idx, :] = vals
                idx += 1
                continue

            # Fill output row
            out[idx, :] = [
                sat_id, lg_label, z_infall, xc_infall, yc_infall, zc_infall,
                vx_infall, vy_infall, vz_infall, Mvir_infall, Mgas_infall, Mstar_infall,
                xc0, yc0, zc0, vx0, vy0, vz0, Mvir0, Mgas0, Mstar0,
                xclg, yclg, zclg, vxlg, vylg, vzlg,
                zbirth, xc_birth, yc_birth, zc_birth, xblg, yblg, zblg, rvir
            ]
            idx += 1

    # Remove rows with NaNs
    out = out[~np.isnan(out[:,2:]).any(axis=1)]

    # Save
    out_file = f'out_infall_{SIM}_{rvir_factor:.1f}rvir.txt'
    header = "halo_id 1:M31/2:MW z_inf Xc_inf Yc_inf Zc_inf Vx_inf Vy_inf Vz_inf Mvir_inf Mgas_inf Mstar_inf Xc_0 Yc_0 Zc_0 Vx_0 Vy_0 Vz_0 Mvir_0 Mgas_0 Mstar_0 XcLG YcLG ZcLG VxLG VyLG VzLG zbirth xcbirth ycbirth zcbirth xlg_birth ylg_birth zlg_birth rvir"
    np.savetxt(out_file, out, fmt='%s', header=header)
    print(f"Saved: {out_file}")


# -----------------------------
# Main loop
# -----------------------------
for SIM in list_SIMS:
    process_simulation(SIM)