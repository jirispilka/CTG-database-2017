# -*- coding: utf-8 -*-
"""
    Univariate feature analysis

    Jiri Spilka, 2019
"""

import numpy as np
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import roc_curve, roc_auc_score
import seaborn as sns
from matplotlib import pyplot as plt

import utils

SELECT_STAGE = 0
NR_SEG = 3

print('STAGE: ', SELECT_STAGE)
print('NR_WINS: ', NR_SEG)

df = utils.load_data_stage_last_k_segments(select_stage=SELECT_STAGE, nr_seg=NR_SEG)

print(f"y == 0: {sum(df.y == 0)}")
print(f"y == 1: {sum(df.y == 1)}")
print(df.info())


df = df.fillna(df.median())

X, y, selected_features = utils.get_X_y_from_dataframe(df)

print('\nList column names: ')
print('\n'.join(list(df)))

print('\nX and y: ')
print(X.shape)
print(y.shape)

print('# positive windows', sum(y == +1))
print('# negative windows', sum(y == 0))

lauc = list()

lmi = mutual_info_classif(X, y)
for i in range(X.shape[1]):
    a = roc_auc_score(y, X[:, i])
    if a < 0.5:
        a = 1 - a

    lauc.append(a)

print(lmi)
print(lauc)

fig, _ = plt.subplots()
ax = plt.subplot(221)
plt.stem(lmi)
ax.set_xticks(range(0, len(selected_features)))
ax.set_xticklabels(selected_features, rotation=90)
fig.tight_layout()
plt.title('Mutual information')

ax = plt.subplot(223)
plt.stem(lauc)
ax.set_xticks(range(0, len(selected_features)))
ax.set_xticklabels(selected_features, rotation=90)
plt.ylim([0.5, 1])
fig.tight_layout()
plt.title('AUC')

plt.subplot(122)
plt.plot(lmi, lauc, 'or')
fig.tight_layout()
plt.title('Mutual information vs AUC')
plt.grid(True)

ind = np.argsort(lauc)
ind = ind[::-1]

print('\n\n ===== Sorted by importance ====')
for i in range(0, len(ind)):
    print('{0:20}: {1:2.2f}'.format(selected_features[ind[i]], lauc[ind[i]]))

names = [s for i, s in enumerate(selected_features) if i in ind[0:10]]

corr = df[names].corr(method='spearman')
f, ax = plt.subplots(figsize=(11, 9))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.heatmap(corr, cmap=cmap, vmin=-1, vmax=1,
            square=True, xticklabels=True, yticklabels=True,
            linewidths=.5, cbar_kws={"shrink": .5}, ax=ax, annot=True)

plt.figure()
for name in names:
    a = roc_auc_score(df['y'], df[name])
    if a > 0.5:
        fpr, tpr, thrs = roc_curve(df['y'], df[name], pos_label=1)
        roc_auc = a
    else:
        fpr, tpr, thrs = roc_curve(df['y'], df[name], pos_label=0)
        roc_auc = 1 - a

    plt.plot(fpr, tpr, lw=2, label='%s = %0.2f' % (name, roc_auc))

plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")

plt.show()

