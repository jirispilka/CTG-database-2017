# -*- coding: utf-8 -*-
"""
    Dummy template for classification
    B(E)3M33UI - Support script for the first semester task

    Jiri Spilka, 2019
"""

from sklearn.dummy import DummyClassifier
from sklearn.metrics import make_scorer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import utils


our_scorer = make_scorer(utils.g_mean_score, greater_is_better=True)


def train_model(X, y):
    """
    Return a trained model.

    Please keep the same arguments: X, y (to be able to import this function for evaluation)
    """
    assert 'X' in locals().keys()
    assert 'y' in locals().keys()
    assert len(locals().keys()) == 2

    sc = StandardScaler()
    clf = DummyClassifier(strategy='stratified')
    pipe = Pipeline(steps=[
        ('sc', sc),
        ('clf', clf)])
    pipe.fit(X, y)
    return pipe


def predict(model1, X):
    """
    Produce predictions for X using given filter.
    Please keep the same arguments: X, y (to be able to import this function for evaluation)
    """
    assert len(locals().keys()) == 2

    return model1.predict(X)


if __name__ == '__main__':

    SELECT_STAGE = 0
    NR_SEG = 1

    df = utils.load_data_stage_last_k_segments(select_stage=SELECT_STAGE, nr_seg=NR_SEG)

    # Demonstration how the model will be used
    df_train, df_test = utils.split_train_test_based_on_year(df)

    X_train, y_train, _ = utils.get_X_y_from_dataframe(df_train)

    print("All data")
    print(f"y == 0: {sum(df.y == 0)}")
    print(f"y == 1: {sum(df.y == 1)}")

    print("\nTraining data")
    print(f"y == 0: {sum(df_train.y == 0)}")
    print(f"y == 1: {sum(df_train.y == 1)}")

    print("\nTest data")
    print(f"y == 0: {sum(df_test.y == 0)}")
    print(f"y == 1: {sum(df_test.y == 1)}")

    # Get test data
    X_test, y_test, _ = utils.get_X_y_from_dataframe(df_test)

    # or you can make a custom train/test split (or CV)
    # X = X_train.copy()
    # X.extend(X_test)
    # y = np.hstack((y_train, y_test))
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    # Train the model
    filter1 = train_model(X_train, y_train)

    # Compute predictions for training data and report our accuracy
    y_tr_pred = predict(filter1, X_train)
    print('\ng-mean on training data: ', utils.g_mean_score(y_train, y_tr_pred))

    # Compute predictions for testing data and report our accuracy
    y_tst_pred = predict(filter1, X_test)
    print('g-mean on test data: ', utils.g_mean_score(y_test, y_tst_pred))
