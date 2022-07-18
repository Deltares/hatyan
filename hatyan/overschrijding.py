"""
All kind of functions to support calculation of probabilities

Author:         Boyan Domhof (Deltares)
Last update:    17/12/2020
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy import optimize, signal
from typing import Union, List
import datetime as dt
import os


def delete_values_between_peak_trough(times_to_delete, times, values):
    mask = np.in1d(times, times_to_delete)
    return times[~mask], values[~mask]


def go_through_peak(values, _i, _i_extra, multiplier):
    while values[_i + _i_extra + multiplier] < values[_i + _i_extra]:
        _i_extra += multiplier
    return _i_extra


def filter_identified_peaks(values, times, times_to_delete):
    _values = values.copy()
    _values[~np.in1d(times, times_to_delete)] = -9999.0
    return _values


def check_peakside(values, _i, multiplier, window, threshold):
    _i_extra = 0
    check = True
    while check:
        try:
            while (values[_i + _i_extra + multiplier] <= values[_i + _i_extra]) and (values[_i + _i_extra + multiplier] != -9999.0):
                _i_extra += multiplier
        except IndexError:
            pass
        try:
            _i1, _i2 = (_i + _i_extra + (multiplier * window)), (_i + _i_extra)
            if any(values[np.min([_i1, _i2]):np.max([_i1, _i2])] > threshold):
                new_peak = values[np.min([_i1, _i2]):np.max([_i1, _i2])].max()
                while values[_i + _i_extra] != new_peak:
                    _i_extra += multiplier
            else:
                check = False
        except IndexError:
            pass
    return _i_extra


def detect_peaks_hkv(df: pd.DataFrame, window: int, inverse: bool = False, threshold: float = None) -> pd.DataFrame:
    _df = df.copy()
    if inverse:
        _df['values'] = -_df['values']
    times, values = _df.index.values, _df['values'].values
    indices = np.arange(len(times))

    __df = _df.sort_values(by=['values'], axis=0, ascending=False)
    times_sorted, values_sorted = __df.index.values, __df['values'].values

    if threshold is None:
        threshold = df['values'].mean() + 2*df['values'].std()
    else:
        if inverse:
            threshold = -threshold

    peaks = np.ones(len(values)) * np.nan
    t_peaks, dt_left_peaks, dt_right_peaks = peaks.copy(), peaks.copy(), peaks.copy()

    if inverse:
        print('Determining peaks (inverse)')
    else:
        print('Determining peaks')
    peak_count = 0
    while len(values_sorted) != 0:
        _t, _p = times_sorted[0], values_sorted[0]
        t_peaks[peak_count], peaks[peak_count] = _t, _p

        _i = indices[_t == times][0]

        # first left peak
        _i_extra = check_peakside(filter_identified_peaks(values, times, times_sorted),
                                  _i, -1, window, threshold)
        dt_left_peaks[peak_count] = (times[_i] - times[_i + _i_extra]) / np.timedelta64(1, 's')
        times_sorted, values_sorted = delete_values_between_peak_trough(times[(_i + _i_extra):(_i+1)],
                                                                        times_sorted, values_sorted)

        # right peak
        _i_extra = check_peakside(filter_identified_peaks(values, times, times_sorted),
                                  _i, 1, window, threshold)
        dt_right_peaks[peak_count] = (times[_i + _i_extra] - times[_i]) / np.timedelta64(1, 's')
        times_sorted, values_sorted = delete_values_between_peak_trough(times[(_i-1):(_i + _i_extra)],
                                                                        times_sorted, values_sorted)

        peak_count += 1

    t_peaks = t_peaks[~np.isnan(t_peaks)]
    peaks = peaks[~np.isnan(peaks)]
    dt_left_peaks = dt_left_peaks[~np.isnan(dt_left_peaks)]
    dt_right_peaks = dt_right_peaks[~np.isnan(dt_right_peaks)]

    dt_total_peaks = dt_left_peaks + dt_right_peaks

    df_peaks = pd.DataFrame(index=pd.to_datetime(t_peaks),
                            data={'values': peaks, 'dt_left': dt_left_peaks,
                                  'dt_right': dt_right_peaks, 'dt_total': dt_total_peaks})
    df_peaks = df_peaks.loc[df_peaks['dt_total'] > 0]

    if inverse:
        df_peaks['values'] = -df_peaks['values']

    return df_peaks


def distribution(df: pd.DataFrame, col: str = None,
                 c: float = -0.3, d: float = 0.4, inverse: bool = False) -> pd.DataFrame:
    col = df.columns[0] if col is None else col
    years = get_total_years(df)
    if inverse:
        df = df.sort_values(by=col, ascending=False)
    else:
        df = df.sort_values(by=col)
    rank = np.array(range(len(df[col]))) + 1
    df['{}_Tfreq'.format(col)] = (1 - (rank + c) / (len(rank) + d)) * (len(rank) / years)
    return df.sort_values(by='{}_Tfreq'.format(col), ascending=False)


def get_weibull(df: pd.DataFrame, threshold: float, Tfreqs: np.ndarray, col: str = None,
                inverse: bool = False) -> pd.DataFrame:
    col = df.columns[0] if col is None else col

    values = df[col].values
    if inverse:
        values = -values
        threshold = -threshold
    p_val_gt_threshold = df[col+'_Tfreq'].loc[values > threshold].iloc[0]

    def pfunc(x, p_val_gt_threshold, threshold, sigma, alpha):
        return p_val_gt_threshold * np.exp(-((x/sigma)**alpha) + ((threshold/sigma)**alpha))

    def pfunc_inverse(p_X_gt_x, p_val_gt_threshold, threshold, sigma, alpha):
        return sigma * (((threshold/sigma)**alpha) - np.log(p_X_gt_x / p_val_gt_threshold))**(1/alpha)

    def der_pfunc(x, p_val_gt_threshold, threshold, alpha, sigma):
        return -p_val_gt_threshold * (alpha * x**(alpha - 1)) * (sigma**(-alpha)) * np.exp(-((x/sigma)**alpha) + ((threshold/sigma)**alpha))

    def cost_func(params, *args):
        return -np.sum([np.log(-der_pfunc(x, args[0], args[1], params[0], params[1])) for x in args[2]])

    initial_guess = np.array([1, abs(threshold)])
    result = optimize.minimize(cost_func,
                               x0=initial_guess,
                               args=(p_val_gt_threshold, threshold, values[values > threshold]),
                               method='Nelder-Mead',
                               options={'maxiter': 1e4})
    if result.success:
        alpha, sigma = result.x[0], result.x[1]
    else:
        raise ValueError(result.message)

    new_values = pfunc_inverse(Tfreqs, p_val_gt_threshold, threshold, sigma, alpha)
    if inverse:
        new_values = -new_values
    return pd.DataFrame(data={'{}_Tfreq'.format(col): Tfreqs,
                              col: new_values}).sort_values(by='{}_Tfreq'.format(col), ascending=False)


def filter_with_threshold(df_raw: pd.DataFrame,
                          df_filtered: pd.DataFrame,
                          threshold: float,
                          inverse: bool = False) -> pd.DataFrame:
    if inverse:
        return pd.concat([df_raw[df_raw['values'] >= threshold],
                          df_filtered[df_filtered['values'] < threshold]], axis=0).sort_index()
    else:
        return pd.concat([df_raw[df_raw['values'] <= threshold],
                          df_filtered[df_filtered['values'] > threshold]], axis=0).sort_index()


def detect_peaks(df: pd.DataFrame,   prominence: int = 10, inverse: bool = False):
    df = df.copy()
    if inverse:
        df['values'] = -1*df['values']
    peak_indices = signal.find_peaks(df['values'].values, prominence=prominence)[0]
    df_peaks = pd.DataFrame(data={'values': df['values'].iloc[peak_indices]},
                            index=df.iloc[peak_indices].index.values)
    threshold = determine_threshold(values=df['values'].values, peak_indices=peak_indices)
    return df_peaks, threshold, peak_indices


def determine_threshold(values: np.ndarray, peak_indices: np.ndarray) -> float:
    w = signal.peak_widths(values, peak_indices)[0]
    for threshold in reversed(range(int(np.floor(values.min())),
                                    int(np.ceil(values.max())))):
        _t = w[values[peak_indices] > threshold]
        if len(_t[_t <= 3]) > (0.1*len(_t)):  # min of 3 tidal periods and at least more than 10%
            break
    return threshold


def get_total_years(df: pd.DataFrame) -> float:
    return (df.index[-1] - df.index[0]).total_seconds() / (3600 * 24 * 365)


def apply_trendanalysis(df: pd.DataFrame, rule_type: str, rule_value: Union[float, dt.datetime]):
    # There are 2 rule types:  - break -> Values before break are removed
    #                          - linear -> Values are increased/lowered based on value in value/year. It is assumes
    #                                      that there is no linear trend at the latest time (so it works its way back
    #                                      in the past). rule_value should be entered as going forward in time
    if rule_type == 'break':
        if not isinstance(rule_value,dt.datetime):
            raise Exception('rule_value should be of instance dt.datetime')
        #return df[dt.datetime.strptime(rule_value, '%d-%M-%Y'):].copy()
        return df[rule_value:].copy()
    elif rule_type == 'linear':
        df, rule_value = df.copy(), float(rule_value)
        dx = np.array([rule_value*x.total_seconds()/(365*24*3600) for x in (df.index[-1] - df.index)])
        df['values'] = df['values'] + dx
        return df
    elif (rule_type == '') or (rule_type is None):
        return df.copy()
    elif isinstance(rule_type, np.float):
        if np.isnan(rule_type):
            return df.copy()
        else:
            raise ValueError('Incorrect rule_type passed to function. Only break or linear are supported')
    else:
        raise ValueError('Incorrect rule_type passed to function. Only break or linear are supported')


def blend_distributions(df_trend: pd.DataFrame, df_weibull: pd.DataFrame, df_hydra: pd.DataFrame = None) -> pd.DataFrame:
    df_trend = df_trend.sort_values(by='values_Tfreq', ascending=False)
    df_weibull = df_weibull.sort_values(by='values_Tfreq', ascending=False)

    # Trend to weibull
    df_blended1 = df_trend.iloc[:-100].copy()
    df_weibull = df_weibull.loc[df_weibull['values_Tfreq'] < df_blended1['values_Tfreq'].iloc[-1]].copy()

    # Weibull to Hydra
    if df_hydra is not None:
        df_hydra = df_hydra.sort_values(by='values_Tfreq', ascending=False)

        Tfreqs_combined = np.unique(np.concatenate((df_weibull['values_Tfreq'].values, df_hydra['values_Tfreq'].values)))
        vals_weibull = np.interp(Tfreqs_combined,
                                 np.flip(df_weibull['values_Tfreq'].values),
                                 np.flip(df_weibull['values'].values))
        vals_hydra = np.interp(Tfreqs_combined,
                               np.flip(df_hydra['values_Tfreq'].values),
                               np.flip(df_hydra['values'].values))

        Tfreq0, TfreqN = df_hydra['values_Tfreq'].values[0], 1/50
        Tfreqs = np.logspace(np.log10(TfreqN), np.log10(Tfreq0), int(1e5))
        vals_weibull = np.interp(np.log10(Tfreqs),
                                 np.log10(np.flip(df_weibull['values_Tfreq'].values)),
                                 np.flip(df_weibull['values'].values))
        vals_hydra = np.interp(np.log10(Tfreqs),
                               np.log10(np.flip(df_hydra['values_Tfreq'].values)),
                               np.flip(df_hydra['values'].values))
        indices = np.arange(len(Tfreqs))
        grads = np.flip(np.arange(len(indices))) / len(indices) * np.pi

        vals_blend = 0.5*(np.cos(grads)+1)*vals_weibull[indices] + (1-0.5*(np.cos(grads)+1))*vals_hydra[indices]

        df_blended2 = pd.DataFrame(data={'values': vals_blend,
                                         'values_Tfreq': Tfreqs}).sort_values(by='values_Tfreq', ascending=False)

        df_blended = pd.concat([df_blended1,
                                df_weibull.loc[(df_weibull['values_Tfreq'] > df_blended2['values_Tfreq'].iloc[0]) &
                                               (df_weibull['values_Tfreq'] < df_blended1['values_Tfreq'].iloc[-1])],
                                df_blended2,
                                df_hydra.loc[df_hydra['values_Tfreq'] < df_blended2['values_Tfreq'].iloc[-1]]], axis=0)
        df_blended = df_blended.drop_duplicates(subset='values_Tfreq').sort_values(by='values_Tfreq', ascending=False)
    else:
        df_blended = pd.concat([df_blended1,
                                df_weibull.loc[(df_weibull['values_Tfreq'] < df_blended1['values_Tfreq'].iloc[-1])]],
                               axis=0).drop_duplicates(subset='values_Tfreq').sort_values(by='values_Tfreq',
                                                                                         ascending=False)

    return df_blended


def plot_distributions(dist: dict, name: str,
                       keys: List[str] = None,
                       color_map: dict = None,
                       xlabel: str = 'Exceedance frequency [1/yrs]',
                       ylabel: str = 'Waterlevel [m]',
                       legend_loc: str = 'lower right'):
    fig, ax = plt.subplots(figsize=(8, 6))
    if keys is None:
        keys = list(dist.keys())
    for k in keys:
        c = color_map[k] if (color_map is not None) and (k in color_map.keys()) else None
        ax.plot(dist[k]['values_Tfreq'], dist[k]['values'], label=k, c=c)
    ax.set_title(name)
    ax.set_xlabel(xlabel), ax.set_xscale('log'), ax.set_xlim([1e-5, 1e3]), ax.invert_xaxis()
    ax.set_ylabel(ylabel)
    ax.legend(fontsize='medium', loc=legend_loc)
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=tuple(i / 10 for i in range(1, 10)), numticks=12))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter()),
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1)) #this was 10, but now meters instead of cm
    ax.yaxis.set_minor_formatter(ticker.NullFormatter()),
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) #to force 2 decimal places
    ax.grid(b=True, which='major'), ax.grid(b=True, which='minor', ls=':')
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig,ax


def interpolate_interested_Tfreqs_to_csv(df: pd.DataFrame, Tfreqs: List[float],
                                         id: str, csv_dir: os.PathLike, prefix: str,) -> pd.DataFrame:
    df_interp = pd.DataFrame(data={'values': np.interp(Tfreqs,
                                                      np.flip(df['values_Tfreq'].values),
                                                      np.flip(df['values'].values)),
                                   'values_Tfreq': Tfreqs}).sort_values(by='values_Tfreq', ascending=False)
    #prefix = os.path.basename(csv_dir)
    df_interp.to_csv(os.path.join(csv_dir, '{}_{}.csv'.format(prefix, id)), index=False, sep=';')
    return df_interp
