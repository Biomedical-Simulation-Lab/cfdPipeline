""" Calc bandedness (and rough SPI) from spectrograms.
"""

from pathlib import Path
import numpy as np 
import pandas as pd 
from bsl.spectral_librosa import Bandedness 
import sys 

if __name__ == "__main__":
    # First arg contains the spectrogram npz files
    # Second is where the output csv is
    # Third is where you want the SBI timeseries to go
    proj_dir = Path(sys.argv[1]) 
    hemo_d = Path(sys.argv[2]) 
    sbi_out_d = Path(sys.argv[2]) 
    spec_data_d = proj_dir / 'spectrogram_data'

    seq = proj_dir.stem

    df = pd.DataFrame()

    for ff in [sbi_out_d]:
        if not ff.exists():
            ff.mkdir(exist_ok=True, parents=True)

    spec_files = sorted(spec_data_d.glob('art*.npz'))

    hemo_summary_file = hemo_d / (seq + '_hemo_params.csv')

    for f in spec_files:
        case = f.stem.split('_')[1]

        bandedness_timeseries = {}
        spec_data = np.load(f)
        S = spec_data['S']

        bins = spec_data['bins']
        freqs = spec_data['freqs']
        sr = spec_data['sr'], 
        n_fft = spec_data['n_fft']
        S[S < -20] = -20 

        b = Bandedness(S=np.exp(S), freqs=freqs, bins=bins, n_fft=n_fft, sr=sr[0], n_chroma=24)
        b.make_chroma_filterbank()
        b.make_chromagram()
        b.calc_chroma_entropy()

        df.at[case, 'SBI'] = b.chroma_entropy.mean()

        bandedness_timeseries['SBI'] = b.chroma_entropy

        np.savez(sbi_out_d / (case + '.npz'), **bandedness_timeseries)

    df_summary = pd.read_csv(hemo_summary_file, dtype={'Case':str})
    df_summary = df_summary.set_index('Case')
    
    df_summary = df_summary.join(df)
    df_summary.to_csv(hemo_summary_file)
