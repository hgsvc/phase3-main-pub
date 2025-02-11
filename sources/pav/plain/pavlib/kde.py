"""
K-mer density estimation.
"""

import importlib.util
import inspect
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.stats

import pavlib

# K-mer orientation matrix: Tig k-mer against forward (vertical axis)
# and reverse (horizontal axis) k-mers in the reference
#
#       F         T
#   ---------------------------------
# F | No data (-1)  | Rev (2)       |
# T | Reference (0) | Fwd & Rev (1) |
#   ---------------------------------
KMER_ORIENTATION_STATE = np.asarray(
    [
        [-1, 2],
        [0, 1]
    ]
)

SAMPLE_INDEX_CHUNK_SIZE = 400


class KdeTruncNorm(object):
    """
    Kernel Density Estimation (KDE) with a truncated-normal distribution. Uses FFT convolution to solve density.
    """

    def __init__(self,
                 bandwidth=pavlib.const.INV_KDE_BANDWIDTH,
                 trunc_z=pavlib.const.INV_KDE_TRUNC_Z,
                 conv=pavlib.const.INV_KDE_FUNC):

        # Check parameters
        if bandwidth <= 0:
            raise ValueError(f'Bandwidth must be > 0: {bandwidth}')

        if trunc_z <= 0:
            raise ValueError(f'Truncation SD must be > 0: {trunc_z}')

        self.bandwidth = float(bandwidth)
        self.band_sd = float(trunc_z)

        # Get convolution method
        if isinstance(conv, str):
            conv_lower = conv.strip().lower()

            is_auto = conv_lower == 'auto'

            self.conv_method = None

            if is_auto:
                conv_lower = 'fft'

            if conv_lower == 'fft' and self.conv_method is None:

                spec = importlib.util.find_spec('scipy.signal')

                if spec is None:
                    if not is_auto:
                        raise RuntimeError(f'Error initializing KdeTruncNorm: Missing package for KDE convolution method {conv}: scipy.signal')
                    else:
                        conv_upper = 'conv'  # Try next

                import scipy.signal
                self.conv_method = scipy.signal.fftconvolve

            if conv_lower == 'conv' and self.conv_method is None:
                self.conv_method = np.convolve

            if self.conv_method is None:
                if is_auto:
                    raise RuntimeError(f'Error initializing KdeTruncNorm: Could not automatically resolve convolution method')

                raise RuntimeError(f'Error initializing KdeTruncNorm: Unknown convolution method: {conv}')

        elif hasattr(conv, '__call__'):
            self.conv_method = conv

        # Inspect convolution method
        n_arg = len(inspect.getfullargspec(self.conv_method).args)
        n_default = len(inspect.getfullargspec(self.conv_method).defaults)

        if n_arg < 2:
            raise RuntimeError(f'Convolution method does not take at least 2 arguments (n = {n_arg}): {self.conv_method}')

        if n_arg - n_default > 2:
            raise RuntimeError(f'Convolution method requires more than 2 arguments (n = {n_arg - n_default} without default values): {self.conv_method}')

        # Set normal vector
        tnorm = scipy.stats.truncnorm(-trunc_z, trunc_z)  # Truncate at Z = band_sd

        self.band_bound = int(np.ceil(trunc_z * bandwidth))  # number of positions in the truncated normal vector

        self.v_tnorm = tnorm.pdf(
            np.arange(-self.band_bound, self.band_bound + 1) / self.bandwidth  # Range from -band_sd to +band_sd after accounting for bandwith
        )  # Pre-compute PDF at positions

        self.v_tnorm = self.v_tnorm / self.v_tnorm.sum()  # Vector sums to 1 (density 1 is a position with ones from -band_sd to band_sd)

    def __call__(self, x, n=None):
        """
        Get the density of an n-length vector with ones at x positions and zeros elsewhere.

        :param x: Location of 1's in the array (if `n` is defined) or an array to directly conolve with the density
            kernel and `n` becomes the length of this array.
        :param n: Length of the array. If `None`, `x` is the full-length array to convolve.

        :return: A density vector of length `n`.
        """

        if n is not None:
            v_state = np.zeros(n)

            for v in x:
                v_state[v] = 1
        else:
            if not isinstance(x, np.ndarray):
                v_state = np.array(x)
            else:
                v_state = x

        y = self.conv_method(v_state, self.v_tnorm)


        return y[self.band_bound:-self.band_bound]


def rl_encoder(df):
    """
    Take a full table of states and KDE estimates per site and generate a table of contiguous states.

    Recall that the state table has columns removed.

    The table returned has these fields:
    * STATE: State of a contiguous region.
    * POS_KDE: Position in the k-mer table.
    * LEN_KDE: Length in the KDE table.
    * POS_QRY: Relative position in the query.
    * LEN_QRY: Length including excluded k-mers
    * MAX_GAIN: Maximum difference between the highest value and the middle (next highest) value in the KDE values.

    :param df: Dataframe of states.
    :param pos_qry: Query position at the start of the table (first query base is 'INDEX' + pos_qry).

    :return: DataFrame with STATE, POS_KDE, LEN_KDE, POS_QRY, LEN_QRY, and MAX_GAIN columns.
    """

    # Stop if dataframe is empty
    if df.shape[0] == 0:
        return pd.DataFrame(
            [], index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
        )

    # Get a table of run-length encoded states (collapse consecutive states)
    df_list = list()

    state = int(df.loc[0, 'STATE'])
    pos = 0

    index = int(df.loc[0, 'INDEX'])

    for i in np.where(df['STATE'][:-1].array != df['STATE'][1:].array)[0] + 1:

        row = df.iloc[i]
        new_index = int(row['INDEX'])

        df_list.append(
            pd.Series(
                [state, pos, i - pos, index, new_index - index],
                index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
            )
        )

        index = new_index
        pos = i
        state = int(row['STATE'])

    df_list.append(
        pd.Series(
            [state, pos, df.shape[0] - pos, index, df.iloc[-1]['INDEX'] - index],
            index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
        )
    )

    df_rl = pd.concat(df_list, axis=1).T.astype(int)

    # Score states
    df_rl['MAX_GAIN'] = 0.0

    for index in range(df_rl.shape[0]):
        row = df_rl.loc[index]

        kde_matrix = np.asarray(
            df.loc[
                row['POS_KDE'] :
                row['POS_KDE'] + row['LEN_KDE'],
                ['KDE_FWD', 'KDE_FWDREV', 'KDE_REV']
            ]
        )

        kde_matrix[:, int(row['STATE'])] - (kde_matrix.sum(axis=1) - kde_matrix.min(axis=1))

        df_rl.loc[index, 'MAX_GAIN'] = np.max(
            2 * kde_matrix[:, int(row['STATE'])] - (kde_matrix.sum(axis=1) - kde_matrix.min(axis=1))
        )

        # If A, B, and C are the KDE values where A <= B <= C (A is the maximum state), then MAX_GAIN is the maximum of A - B:
        # A - B = A - B + (A - A) + (C - C)
        #       = A + A - (A + B + C) + C
        #       = 2A - (A + B + C - C)
        #       = 2A - (S - C), where S is the sum across all state (approx. 1) and C is the minimum value

    return df_rl
