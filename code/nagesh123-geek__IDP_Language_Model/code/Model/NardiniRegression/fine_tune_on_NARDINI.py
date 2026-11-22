import numpy as np
import pandas as pd
import json

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

# ---------------------------------------------------
# Load embeddings from JSON
# ---------------------------------------------------

#X = np.load("embeddings.npy")
with open("../../../seq_500/local_embeddings.json", "r") as f:
    X = np.array(json.load(f), dtype=np.float32)

print("Embeddings shape:", X.shape,flush=True) # shape: (N, 512)

# ---------------------------------------------------
# Load regression targets
# ---------------------------------------------------
y_df = pd.read_csv(
    "../../features_len_le_500.txt",
    sep="\t"
)

feature_names = y_df.columns.tolist()

y = y_df.values.astype(np.float32)

print("Embeddings shape:", X.shape,flush=True)
print("Targets shape:", y.shape,flush=True)

# ---------------------------------------------------
# Train / Test split
# ---------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42
)

# ---------------------------------------------------
# Scale embeddings
# ---------------------------------------------------
x_scaler = StandardScaler()

X_train = x_scaler.fit_transform(X_train)
X_test = x_scaler.transform(X_test)

# ---------------------------------------------------
# Scale targets (important for MLP regression)
# ---------------------------------------------------
y_scaler = StandardScaler()

y_train = y_scaler.fit_transform(y_train)
y_test_scaled = y_scaler.transform(y_test)

# ---------------------------------------------------
# Convert to tensors
# ---------------------------------------------------
X_train_t = torch.tensor(X_train, dtype=torch.float32)
X_test_t = torch.tensor(X_test, dtype=torch.float32)

y_train_t = torch.tensor(y_train, dtype=torch.float32)
y_test_t = torch.tensor(y_test_scaled, dtype=torch.float32)

# ---------------------------------------------------
# DataLoader
# ---------------------------------------------------
batch_size = 1024

train_loader = DataLoader(
    TensorDataset(X_train_t, y_train_t),
    batch_size=batch_size,
    shuffle=True
)

# ---------------------------------------------------
# MLP Model
# ---------------------------------------------------
class MLPRegressor(nn.Module):

    def __init__(self, input_dim, output_dim):
        super().__init__()

        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, 1024),
            nn.ReLU(),
            nn.Dropout(0.2),

            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Dropout(0.2),

            nn.Linear(512, 256),
            nn.ReLU()
        )

        self.regressor = nn.Linear(256, output_dim)

    def forward(self, x, return_embedding=False):

        emb = self.feature_extractor(x)

        if return_embedding:
            return emb

        out = self.regressor(emb)

        return out
# ---------------------------------------------------
# Device
# ---------------------------------------------------
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

model = MLPRegressor(
    input_dim=X.shape[1],
    output_dim=y.shape[1]
).to(device)

# ---------------------------------------------------
# Loss + Optimizer
# ---------------------------------------------------
criterion = nn.MSELoss()

optimizer = torch.optim.Adam(
    model.parameters(),
    lr=1e-4
)

# ---------------------------------------------------
# Training
# ---------------------------------------------------
epochs = 100000

for epoch in range(epochs):

    model.train()

    total_loss = 0

    for xb, yb in train_loader:

        xb = xb.to(device)
        yb = yb.to(device)

        optimizer.zero_grad()

        pred = model(xb)

        loss = criterion(pred, yb)

        loss.backward()

        optimizer.step()

        total_loss += loss.item()

    avg_loss = total_loss / len(train_loader)

    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1:3d} | Loss: {avg_loss:.6f}",flush=True)

# ---------------------------------------------------
# Evaluation
# ---------------------------------------------------
model.eval()

with torch.no_grad():

    pred_scaled = model(X_test_t.to(device)).cpu().numpy()

# inverse transform
y_pred = y_scaler.inverse_transform(pred_scaled)

# ---------------------------------------------------
# Metrics
# ---------------------------------------------------
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print("\nOverall Results",flush=True)
print("MSE :", mse,flush=True)
print("R2  :", r2,flush=True)

# ---------------------------------------------------
# Per-feature R2
# ---------------------------------------------------
print("\nPer-feature R2:\n",flush=True)

for i, name in enumerate(feature_names):

    r2_i = r2_score(y_test[:, i], y_pred[:, i])

    print(f"{name:25s} : {r2_i:.4f}",flush=True)

# ---------------------------------------------------
# Save predictions
# ---------------------------------------------------
pred_df = pd.DataFrame(
    y_pred,
    columns=feature_names
)

pred_df.to_csv(
    "mlp_predicted_features.txt",
    sep="\t",
    index=False
)

print("\nSaved predictions:",flush=True)
print("mlp_predicted_features.txt",flush=True)


# ---------------------------------------------------
# Save trained model
# ---------------------------------------------------
torch.save(model.state_dict(), "mlp_regressor_50k.pth")

# ---------------------------------------------------
# Save scalers
# ---------------------------------------------------
import joblib

joblib.dump(x_scaler, "x_scaler.pkl")
joblib.dump(y_scaler, "y_scaler.pkl")

# ---------------------------------------------------
# Save feature names
# ---------------------------------------------------
with open("feature_names.txt", "w") as f:
    for name in feature_names:
        f.write(name + "\n")

print("\nSaved:")
print(" - mlp_regressor.pth")
print(" - x_scaler.pkl")
print(" - y_scaler.pkl")
print(" - feature_names.txt")
