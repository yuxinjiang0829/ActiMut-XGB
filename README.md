# ActiMut-XGB
This repository presents an XGBoost-based model for predicting protein thermostability using sequence-based features. Integrating ESM2 embeddings, physicochemical, evolutionary, and positional features, the model is enhanced by transfer learning and validated on CALB mutants, offering insights for protein engineering.

# File Organization
./data: Contains data used in ActiMut-XGB.
./data/seq: Contains sequence data for wild types and mutants used by fireprot,calb,myoglobin, including sequences of different lengths
./features: Contains features used by fireprot,calb,myoglobin, including Integrating ESM2 embeddings, physicochemical, evolutionary, and positional features
./weight: Contains saved model weights.
./model: Contains the code needed to train the model
