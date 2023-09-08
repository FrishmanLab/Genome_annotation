import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import time
import numpy as np
from sklearn import datasets, cluster
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.ensemble import AdaBoostClassifier, GradientBoostingClassifier, RandomForestClassifier, StackingClassifier, VotingClassifier

def load_files(first_file, second_file):
    # Load the true protein sequences
    sequences1 = []
    with open(first_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            sequences1.append(line.strip())

    # Load the decoy protein sequences
    sequences2 = []
    with open(second_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            sequences2.append(line.strip())



    # remove any sequences with not usual amino acids
    allowed_letters = set('acgtn*')


    filtered_sequences1 = []
    for seq in sequences1:
        if set(seq.lower()) <= allowed_letters:
            filtered_sequences1.append(seq)

    filtered_sequences2 = []
    for seq in sequences2:
        if set(seq.lower()) <= allowed_letters:
            filtered_sequences2.append(seq)

    # Combine the true and decoy sequences and create labels
    sequences = filtered_sequences1 + filtered_sequences2
    labels = np.concatenate([np.ones(len(filtered_sequences1)), np.zeros(len(filtered_sequences2))])
    return sequences, labels

def train_model(first_file, second_file):

    sequences, labels = load_files(first_file, second_file)
    k = 3  # k-mer length
    # Divide each element in the list into kmers of length k and join into a string with spaces
    seqs = [' '.join([sequences[i][j:j+k] for j in range(0, len(sequences[i]), k)]) for i in range(len(sequences))]

    # split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(seqs, labels, test_size=0.2, random_state=42)

    # convert text data to feature vectors
    vectorizer = CountVectorizer()
    X_train = vectorizer.fit_transform(X_train)
    X_test = vectorizer.transform(X_test)

    print("X_train shape:", X_train.shape) 
    print("X_test shape:", X_test.shape) 

    # train the model
    SGD_svm = SGDClassifier(loss='hinge', penalty='l2', alpha=0.0001, max_iter=1000, random_state=42)
    SGD_start_time = time.time()
    SGD_svm_model = SGD_svm.fit(X_train, y_train)
    SGD_end_time = time.time()

    # Prediction
    start_time_pred = time.time()
    SGD_svm_pred = SGD_svm.predict(X_test)
    end_time_pred = time.time()

    # Evaluate the classifier on test setacc = accuracy_score(y_test, SGD_svm_pred)
    acc = accuracy_score(y_test, SGD_svm_pred)
    prec = precision_score(y_test, SGD_svm_pred)
    rec = recall_score(y_test, SGD_svm_pred)
    f1 = f1_score(y_test, SGD_svm_pred)
    print("SGD model performance:")
    print("Accuracy:", acc)
    print("Precision:", prec)
    print("Recall:", rec)
    print("F1 score:", f1)


def main(argv):
    first_file = argv[0]
    second_file = argv[1]
    train_model(first_file, second_file)

if __name__ == "__main__":
    main(sys.argv[1:])
