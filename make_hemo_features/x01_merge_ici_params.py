""" Merge WSS param csv files.
"""

import pandas as pd
import sys
from pathlib import Path 

if __name__ == "__main__":
    # search_dir = Path(sys.argv[1])
    search_dir = Path('/scratch/s/steinman/macdo708/surge_cfd_analysis/data/ici_data')
    csv_files = sorted(search_dir.glob('art*.csv')) # Bug - csv files are named as npz, but can open this way.
    df_out_file = search_dir / 'ici_params.csv'

    df_out = pd.DataFrame()

    # Split list into pairs
    chunk_size = 2
    paired_list = [csv_files[i:i+chunk_size] for i in range(0, len(csv_files), chunk_size)]

    # Double check all pairs match
    for pair in paired_list:
        assert pair[0].stem.split('_')[1] == pair[1].stem.split('_')[1]

    case_names = sorted(list(set([x.stem.split('_')[1] for x in csv_files])))
    
    for pair in paired_list:
        case = pair[0].stem.split('_')[1]
        df_std = pd.read_csv(pair[0], index_col=0)
        df_srg = pd.read_csv(pair[1], index_col=0)

        df_out.at[case, 'ICI_std'] = df_std['ICI'].values[0]
        df_out.at[case, 'ICI_srg'] = df_srg['ICI'].values[0]

        df_out.at[case, 'u_in_mean_std'] = df_std['u_in_mean'].values[0]
        df_out.at[case, 'u_in_mean_srg'] = df_srg['u_in_mean'].values[0]

        df_out.at[case, 'u_in_mean_std'] = df_std['u_in_mean'].values[0]
        df_out.at[case, 'u_in_mean_srg'] = df_srg['u_in_mean'].values[0]

        df_out.at[case, 'u_n_in_mean_std'] = df_std['u_n_in_mean'].values[0]
        df_out.at[case, 'u_n_in_mean_srg'] = df_srg['u_n_in_mean'].values[0]

        df_out.at[case, 'u_n_out_mean_std'] = df_std['u_n_out_mean'].values[0]
        df_out.at[case, 'u_n_out_mean_srg'] = df_srg['u_n_out_mean'].values[0]

        df_out.at[case, 'a_in_mean_std'] = df_std['a_in_mean'].values[0]
        df_out.at[case, 'a_in_mean_srg'] = df_srg['a_in_mean'].values[0]

        df_out.at[case, 'a_out_mean_std'] = df_std['a_out_mean'].values[0]
        df_out.at[case, 'a_out_mean_srg'] = df_srg['a_out_mean'].values[0]

        df_out.at[case, 'Q_i_mean_std'] = df_std['Q_i_mean'].values[0]
        df_out.at[case, 'Q_i_mean_srg'] = df_srg['Q_i_mean'].values[0]

        df_out.at[case, 'Q_v_mean_std'] = df_std['Q_v_mean'].values[0]
        df_out.at[case, 'Q_v_mean_srg'] = df_srg['Q_v_mean'].values[0]

        df_out.at[case, 'A_i_mean_std'] = df_std['A_i_mean'].values[0]
        df_out.at[case, 'A_i_mean_srg'] = df_srg['A_i_mean'].values[0]

        df_out.at[case, 'A_o_std'] = df_std['A_o'].values[0]
        df_out.at[case, 'A_o_srg'] = df_srg['A_o'].values[0]

    df_out.to_csv(df_out_file)