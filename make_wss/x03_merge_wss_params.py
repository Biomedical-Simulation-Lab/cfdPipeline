""" Merge WSS param csv files.
"""

import pandas as pd
import sys
from pathlib import Path 

if __name__ == "__main__":
    # search_dir = Path(sys.argv[1])
    search_dir = Path('/scratch/s/steinman/macdo708/surge_cfd_analysis/data/wss_params_per_case')
    csv_files = sorted(search_dir.glob('art*.csv'))
    df_out_file = search_dir / 'wss_params.csv'

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

        df_out.at[case, 'TAWSS_sac_std'] = df_std['TAWSS_sac'].values[0]
        df_out.at[case, 'TAWSS_sac_srg'] = df_srg['TAWSS_sac'].values[0]
        df_out.at[case, 'TAWSS_parent_std'] = df_std['TAWSS_parent'].values[0]
        df_out.at[case, 'TAWSS_parent_srg'] = df_srg['TAWSS_parent'].values[0]
        df_out.at[case, 'SCI_std'] = df_std['SCI'].values[0]
        df_out.at[case, 'SCI_srg'] = df_srg['SCI'].values[0]
        df_out.at[case, 'LSA_std'] = df_std['LSA'].values[0]
        df_out.at[case, 'LSA_srg'] = df_srg['LSA'].values[0]

    df_out.to_csv(df_out_file)