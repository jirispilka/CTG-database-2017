"""
    Utilities for data load and selection

    B(E)3M33UI - Support script for the first semester task

    Jiri Spilka, 2019
"""
import numpy as np
import pandas as pd
import pathlib

from sklearn.metrics import confusion_matrix

FEAT_DIR = pathlib.Path(__file__).parent / '..' / 'features'

CSV_STAT = 'Features_CTU_stat_spectral_figo_mfhurst_20190329.csv'


PH_THR = 7.05


def load_data_binary():

    df = pd.read_csv(str(FEAT_DIR / CSV_STAT))
    df['y'] = (df.pH <= 7.05).astype(int).ravel()

    return df


def load_data_stage_last_k_segments(select_stage=0, nr_seg=1):
    """
    Load k last segments from data

    :param select_stage: 0 - all, 1 - first stage, 2 - second stage
    :param nr_seg: number of last segments to load
    :return:
    """

    df = load_data_binary()

    if select_stage == 0:
        return df.loc[df.segIndex <= nr_seg, :]

    elif select_stage == 1:
        ind = np.logical_and(df.segStageI_index > 0, df.segStageI_index <= nr_seg)
        return df.loc[ind, :]

    elif select_stage == 2:
        ind = np.logical_and(df.segStageII_index > 0, df.segStageII_index <= nr_seg)
        return df.loc[ind, :]

    else:
        raise Exception(f'Unknown value select_stage={select_stage}')


def get_X_y_from_dataframe(df):
    """
    Get feature matrix and labels
    :return:
    """

    y = (df.pH <= 7.05).astype(int).ravel()

    df = df.drop(columns=['name', 'pH', 'year', 'segStart_samp', 'segEnd_samp', 'segIndex', 'segStageI_index',
                          'segStageII_index', 'y'])

    # the stage information might be useful
    df = df.drop(columns=['segStage'])

    # other features that are probably not very useful (correlated to the other ones or irrelevant)
    df = df.drop(columns=['bslnMean', 'bslnSD', 'decDeltaMedian', 'decDeltaMad', 'decDtrdPlus', 'decDtrdMinus',
                          'decDtrdMedian', 'bslnAllBeta0', 'bslnAllBeta1',
                          'MF_hmin_noint', 'H310', 'MF_c1', 'MF_c2', 'MF_c3', 'MF_c4'])

    X = df.get_values()
    feat_names = list(df)

    return X, y, feat_names


def split_train_test_based_on_year(df):

    df_train = df.loc[df.year < 2017, :]
    df_test = df.loc[df.year >= 2017, :]

    return df_train, df_test


def g_mean(estimator, X, y):
    y_pred = estimator.predict(X)
    return g_mean_score(y, y_pred)


def g_mean_score(y, y_pred):
    """Return a modified accuracy score with larger weight of false positives."""

    cm = confusion_matrix(y, y_pred)
    if cm.shape != (2, 2):
        raise ValueError('The ground truth values and the predictions may contain at most 2 values (classes).')

    tn = cm[0, 0]
    fn = cm[1, 0]
    tp = cm[1, 1]
    fp = cm[0, 1]

    se = tp/float(tp + fn)
    sp = tn/float(tn + fp)
    # metrics.
    g = np.sqrt(se * sp)

    return g
