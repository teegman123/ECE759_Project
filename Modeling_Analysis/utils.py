import os
import json
import numpy as np
import pandas as pd
from biomart import BiomartServer
import re
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, KFold, cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB, MultinomialNB, BernoulliNB
from sklearn.ensemble import RandomForestClassifier


root = os.getcwd() # gives path to git clone
data_dir = os.path.join(root, "Data")
fig_dir = os.path.join(root, "Figures")


def load_json(filepath):
    """
    Loads JSON data from a file.

    Args:
        filepath (str): The path to the JSON file.

    Returns:
        dict or list: The JSON data as a Python dictionary or list, or None if an error occurs.
    """
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in: {filepath}")
        return None


def motif_data_to_json(input_dir, output_dir):
    with open(input_dir, "r") as file:
        result = file.readlines()

    motifs = pd.Series(result[0][:-1].split(',')[1:])

    data = {}
    for i in range(1, len(result)):
        gene = result[i][:-1].split(',')[0]
        related_motifs = motifs[np.array(result[i][:-1].split(',')[1:], dtype=int) == 1]
        data[gene] = list(related_motifs)

    # Write the dictionary to a JSON file
    with open(output_dir, "w") as file:
        json.dump(data, file, indent=4) # indent for better readability

    print(f"Dictionary written to {output_dir}")


def balance_response_vector(vec_, abs_boundary=0.5, prop=1.0, seed=None):
    """
    Randomly balances a vector by trimming the majority class, with classes
    defined as elements with absolute value above or below a boundary.

    Parameters:
        vec (list or np.ndarray or pd.Series): Input binary vector (0s and 1s).
        seed (int, optional): Seed for reproducibility.

    Returns:
        np.ndarray: Balanced binary vector.
    """
    if seed is not None:
        np.random.seed(seed)

    vec = np.array(vec_)
    # if not np.all(np.isin(vec, [0, 1])):
    #     raise ValueError("Input vector must contain only 0s and 1s.")

    indices_0 = np.where(np.abs(vec) < abs_boundary)[0]
    indices_1 = np.where(np.abs(vec) >= abs_boundary)[0]

    min_len = min(len(indices_0), len(indices_1))
    sampled_0 = np.random.choice(
        indices_0, min(len(indices_0), int(min_len/prop)), replace=False)
    sampled_1 = np.random.choice(
        indices_1, min(len(indices_1), int(min_len/prop)), replace=False)

    balanced_indices = np.concatenate([sampled_0, sampled_1])
    balanced_indices = np.sort(balanced_indices)

    return vec_.iloc[balanced_indices]


def overwrite_excel_sheet(file_path, sheet_name, data):
    """
    Overwrites an existing sheet in an Excel file with new data.

    Args:
        file_path (str): The path to the Excel file.
        sheet_name (str): The name of the sheet to overwrite.
        data (pd.DataFrame): The DataFrame containing the data to write.
    """
    if os.path.exists(file_path):
      try:
        with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            data.to_excel(writer, sheet_name=sheet_name, index=False)
      except Exception as e:
        print(f"An error occurred: {e}")
    else:
        with pd.ExcelWriter(file_path, engine='openpyxl', mode='w') as writer:
            data.to_excel(writer, sheet_name=sheet_name, index=False)

def download_motif_info(motifs_list):
    values = []

    for motif in motifs_list:
        # Look for Arabidopsis gene ID (e.g., AT1G12630)
        gene_match = re.search(r'(AT\dG\d+)', motif.upper())
        if gene_match:
            values.append(np.array([gene_match.group(1), None]))
        else:
            # If no AGI, take the part after the last dot before _col
            symbol_match = re.search(r'\.([A-Za-z0-9]+)_col', motif)
            if symbol_match:
                values.append(np.array([None, symbol_match.group(1).upper()]))

    values = np.array(values)
    motifs = pd.DataFrame(index=range(len(motifs_list)))
    motifs["Motif Description"] = motifs_list
    motifs["AGI"] = values[:, 0]
    motifs["GeneSymbol"] = values[:, 1]
    motifs["Description"] = None

    # Connect to Ensembl Plants BioMart
    server = BiomartServer("http://plants.ensembl.org/biomart/martservice")
    dataset = server.datasets['athaliana_eg_gene']

    # Query BioMart for gene name to AGI mapping
    symbols = motifs["GeneSymbol"].tolist()
    agis = motifs["AGI"].tolist()

    response = dataset.search({
        'filters': {
            'external_gene_name': [*symbols, *agis]
        },
        'attributes': [
            'external_gene_name',  # Input symbol
            'ensembl_gene_id',     # AGI code
            'description'          # Gene description
        ]
    })

    # Parse the results
    mapped = [line.split("\t") for line in response.text.strip().split("\n")]
    mapped = pd.DataFrame(mapped, columns=["Input Symbol", "AGI Code", "Description"])

    mapped.set_index("Input Symbol", inplace=True, drop=True)
    for i in motifs.index:
        symbol = motifs.at[i, "GeneSymbol"]
        agi = motifs.at[i, "AGI"]
        if symbol in mapped.index:
            motifs.at[i, "Description"] = mapped.loc[symbol, "Description"]
            motifs.at[i, "AGI"] = mapped.loc[symbol, "AGI Code"]
        elif agi in mapped.index:
            motifs.at[i, "Description"] = mapped.loc[agi, "Description"]

    return motifs


def plot_hyperparameter_selection(
        values,
        scores,
        accepting_criterion,
        selected_ind,
        x_label,
        caption=None,
        headline=None,
        save=True
):
    max_score = np.max(scores)
    max_score_ind = np.argmax(scores)
    best_value = values[selected_ind]

    # Create subplots (1 row, 2 columns)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    if accepting_criterion.endswith("jump"):
        threshold = float(accepting_criterion.split(' ')[1])
    elif accepting_criterion == "max":
        ...
    elif accepting_criterion.endswith("margin"):
        threshold = float(accepting_criterion.split(' ')[2])
        axes[0].axhspan(
            max_score - threshold, max_score,
            facecolor='green', alpha=0.2, label="acceptable region")
        axes[0].scatter(
            [values[max_score_ind]], [max_score], marker='x',
            color='green', label=f"max accuracy = %{100 * max_score:0.1f}")
        axes[1].axhspan(
            max_score - threshold, max_score,
            facecolor='green', alpha=0.2, label="acceptable region")
        axes[1].scatter(
            [np.log10(values[max_score_ind])], [max_score], marker='x',
            color='green', label=f"max accuracy = %{100 * max_score:0.1f}")
    else:
        ...

    n = len(values)
    rng = np.array(sorted(list({*range(0, n, int(np.sqrt(n))), n - 1, selected_ind})))

    # Add the main title
    if headline is not None:
        fig.suptitle(headline, fontsize=16)
        fig.text(0.5, 0.92, f'selection criterion: {accepting_criterion}',
                 ha="center", fontsize=8, color="gray")
    if caption is not None:
        fig.text(0.5, -0.1, caption, ha="center", fontsize=12, color="gray")

    # Plotting with Seaborn in subplots
    sns.lineplot(x=np.array(values), y=scores, ax=axes[0])
    axes[0].vlines(
        x=best_value, ymin=0.5, ymax=1.0,
        colors='r', linestyles='dashed', label=f"{x_label} = {best_value:0.2f}")
    sns.scatterplot(
        x=np.array(values)[rng],
        y=np.array(scores)[rng],
        ax=axes[0])
    axes[0].set_title(f'linear x scale')
    axes[0].set_xlabel(f'{x_label}')
    axes[0].set_ylabel('accuracy')
    axes[0].legend()

    sns.lineplot(x=np.log10(values), y=scores, ax=axes[1])
    axes[1].vlines(
        x=np.log10(best_value), ymin=0.5, ymax=1.0,
        colors='r', linestyles='dashed', label=f"{x_label} = {best_value:0.2f}")
    sns.scatterplot(
        x=np.log10(values)[rng],
        y=np.array(scores)[rng],
        ax=axes[1])
    axes[1].set_title('logarithmic x scale')
    axes[1].set_xlabel(f'log({x_label})')
    axes[1].set_ylabel('accuracy')
    axes[1].legend()

    plt.tight_layout()

    if save:
        output_dir = os.path.join(fig_dir, headline + ".png")
        plt.savefig(output_dir, bbox_inches='tight', dpi=300)

    plt.show()


def tune_logistic_regression(
        X_, y_, c_values,
        k, accepting_criterion,
        penalty, l1_ratio=None, random_state=None):
    """
    Tunes the hyperparameter C for L1-penalized logistic regression using K-fold CV.

    Parameters:
        X_ (pd.DataFrame): Feature matrix.
        y_ (pd.Series or np.array): Binary labels.
        c_values (list or np.array): List of C values to try.
        k (int): Number of folds for cross-validation.
        accepting_criterion (str): Since the accuracy increases with larger C values, take the best C that is accepted based on some accepting_acc method.
            max
            first in 0.01 margin
            last 0.002 jump
        penalty (str): penalty
        l1_ratio (float): l1_ratio
        random_state (int, optional): Random seed for reproducibility.

    Returns:
        dict: Contains best C, final test accuracy, and the trained model.
    """

    kf = KFold(n_splits=k, shuffle=True, random_state=random_state)

    scores = []
    for c in c_values:
        val_scores = []

        for train_idx, val_idx in kf.split(X_):
            X_fold_train, X_fold_val = X_.iloc[train_idx], X_.iloc[val_idx]
            y_fold_train, y_fold_val = y_.iloc[train_idx], y_.iloc[val_idx]

            model = LogisticRegression(
                penalty=penalty,
                C=c,
                solver='saga' if penalty == 'elasticnet' else 'liblinear',
                max_iter=10000,
                l1_ratio=l1_ratio)
            model.fit(X_fold_train, y_fold_train)
            y_pred = model.predict(X_fold_val)
            val_scores.append(accuracy_score(y_fold_val, y_pred))

        avg_score = np.mean(val_scores)
        scores.append(avg_score)

    if accepting_criterion.endswith("jump"):
        threshold = float(accepting_criterion.split(' ')[1])
        ind = np.max(np.argwhere((np.diff(scores))/np.max(scores) > threshold))
    elif accepting_criterion == "max":
        ind = np.argmax(scores)
    elif accepting_criterion.endswith("margin"):
        max_score = np.max(scores)
        threshold = float(accepting_criterion.split(' ')[2])
        ind = np.min(np.argwhere(scores - max_score >= -threshold))
    else:
        ind = 0

    return ind, c_values, np.array(scores)


def tune_random_forest(
        X_, y_, k, accepting_criterion,
        n_estimators=None, max_depths=None, random_state=None):
    """
    Step 1: Find the smallest subset of features that gives near-optimal accuracy using Random Forest.
    Step 2: Tune n_estimators and max_depth on that subset.

    Parameters:
        X (pd.DataFrame): Sparse feature matrix.
        y (pd.Series): Binary target.
        k: k used for cross validation.
        accepting_criterion (str): ...
        n_estimators: list of n_estimator values for Step 2.
        max_depths: list of max_depth values for Step 2.
        random_state (int): Random seed.

    Returns:
        best_model (RandomForestClassifier): Fine-tuned classifier.
        best_score (float): Best accuracy score.
        selected_features (list): Names of selected features.
    """

    if max_depths is None:
        max_depths = [None, 10, 20, 30]
    if n_estimators is None:
        n_estimators = [50, 100, 200, 350, 500]

    # Step 1: Fit baseline Random Forest
    baseline_rf = RandomForestClassifier(
        n_estimators=100, random_state=random_state, n_jobs=-1)
    baseline_cv = cross_val_score(baseline_rf, X_, y_, cv=k, scoring='accuracy')
    baseline_score = baseline_cv.mean()
    print(f"Baseline accuracy (all {X_.shape[1]} features): {baseline_score:.4f}")

    baseline_rf.fit(X_, y_)
    features_importance = baseline_rf.feature_importances_
    feature_ranks = np.argsort(features_importance)[::-1]

    step = int(np.sqrt(X_.shape[1]))
    features_scores, features_n = [], []

    for top_n in sorted(list({int(_) for _ in np.logspace(0,np.log10(X_.shape[1]),step)})):
        top_indices = feature_ranks[:top_n]
        X_subset = X_.iloc[:, top_indices]
        scores = cross_val_score(
            RandomForestClassifier(
                n_estimators=100, random_state=random_state, n_jobs=-1),
            X_subset, y_, cv=k, scoring='accuracy')
        score = scores.mean()
        features_scores.append(score)
        features_n.append(top_n)

    max_score = np.max(features_scores)
    if accepting_criterion.endswith("jump"):
        threshold = float(accepting_criterion.split(' ')[1])
        ind = np.max(np.argwhere((np.diff(features_scores))/max_score > threshold))
    elif accepting_criterion == "max":
        ind = np.argmax(features_scores)
    elif accepting_criterion.endswith("margin"):
        threshold = float(accepting_criterion.split(' ')[2])
        ind = np.min(np.argwhere(features_scores - max_score >= -threshold))
    else:
        ind = 0

    best_top_n = features_n[ind]
    sufficient_features = X_.columns[feature_ranks[:best_top_n]]

    # Step 2: Hyperparameter tuning on the selected subset
    X_selected = X_[sufficient_features]
    param_grid = {
        'n_estimators': n_estimators,
        'max_depth': max_depths
    }
    rf = RandomForestClassifier(random_state=random_state, n_jobs=-1)
    grid = GridSearchCV(rf, param_grid, cv=k, scoring='accuracy', n_jobs=-1)
    grid.fit(X_selected, y_)

    best_model = grid.best_estimator_
    best_score = grid.best_score_

    print(f"Best tuned RF accuracy: {best_score:.4f}")
    print(f"Best params: {grid.best_params_}")

    return ind, features_n, features_scores, best_model, list(sufficient_features)
