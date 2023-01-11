""" Merge neck features
"""

import pandas as pd
import sys
from pathlib import Path 

if __name__ == "__main__":
    search_dir = Path(sys.argv[1])

    csv_files = sorted(search_dir.glob('*.csv'))
    case_names = [x.stem.split('_')[0] for x in csv_files]
    df_out_file = search_dir / (search_dir.parents[1].stem + '_morph_params.csv')

    dfs = [pd.read_csv(x, index_col=0) for x in csv_files]

    for cdx, c in enumerate(case_names):
        dfs[cdx]['Case'] = c
        dfs[cdx].set_index('Case', inplace=True)

    df = pd.concat(dfs)
    df.to_csv(df_out_file)