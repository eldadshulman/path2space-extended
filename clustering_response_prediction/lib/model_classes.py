"""Model classes used by the response-prediction notebooks.

These are the exact estimators saved into the committed `.joblib` model files,
lifted from `scr/new_tras_clusters/auc_final_chemo_v3.py` and
`scr/new_tras_clusters/auc_final_trastuzumab_new.py`. The classes must live in
an importable module (not a notebook) so that `joblib.load(...)` can reconstruct
the pickled pipelines.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin


class EnsembleModel(BaseEstimator, ClassifierMixin):
    """Average of `predict_proba(X)[:, 1]` across the 5 fold-trained estimators."""

    def __init__(self, models):
        self.models = models

    def predict_proba(self, X):
        return np.vstack(
            [m.predict_proba(X)[:, 1] for m in self.models]
        ).mean(axis=0)

    def predict(self, X):
        return (self.predict_proba(X) > 0.5).astype(int)


class FeatureAligner(BaseEstimator, TransformerMixin):
    """Make sure the query matrix has the columns the model was trained on.

    Missing columns are filled with `fill_values[col]` if provided, else 0.
    """

    def __init__(self, expected_features, fill_values=None):
        self.expected_features = expected_features
        self.fill_values = fill_values

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        out = pd.DataFrame(index=X.index)
        for f in self.expected_features:
            if f in X.columns:
                out[f] = X[f]
            elif self.fill_values is not None and f in self.fill_values:
                out[f] = self.fill_values[f]
            else:
                out[f] = 0
        return out


class CombinedModel(BaseEstimator, ClassifierMixin):
    """Mean of the cluster-based ensemble proba and the MinMax-scaled HER2 SPAND.

    Used by the trastuzumab notebook to combine spatial-cluster signal with
    HER2 heterogeneity. The HER2 SPAND input is scaled once per fold during
    training (the list of fitted MinMaxScalers is passed in) and the test-time
    scaled values are averaged across folds.
    """

    def __init__(self, ensemble_pipeline, scalers_het, spand_feature="her2_spand"):
        self.ensemble_pipeline = ensemble_pipeline
        self.scalers_het = scalers_het
        self.spand_feature = spand_feature

    def predict_proba(self, X, return_spand_only=False):
        if self.spand_feature not in X.columns:
            raise ValueError(f"Missing required feature: {self.spand_feature}")
        X_spand = X[[self.spand_feature]]
        scaled_all = np.array([s.transform(X_spand) for s in self.scalers_het])
        scaled_mean = scaled_all.mean(axis=0).squeeze()
        scaled_mean = np.clip(scaled_mean, 0, 1)

        if return_spand_only:
            return scaled_mean

        # Cluster ensemble prediction
        X_clust = self.ensemble_pipeline.named_steps["aligner"].transform(X)
        proba_clust = self.ensemble_pipeline.named_steps["ensemble"].predict_proba(X_clust)
        return (proba_clust + scaled_mean) / 2

    def predict(self, X):
        return (self.predict_proba(X) > 0.5).astype(int)
