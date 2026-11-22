import numpy as np
import torch
import torch.nn as nn
import joblib
import json

# ---------------------------------------------------
# Load embeddings
# ---------------------------------------------------
#X = np.load("embeddings.npy")

#X = np.load("embeddings.npy")
with open("../../../seq_500/local_embeddings.json", "r") as f:
    X = np.array(json.load(f), dtype=np.float32)

print("Embeddings shape:", X.shape,flush=True) # shape: (N, 512)


# ----------------------------
# Load embeddings
# ----------------------------
#with open("embeddings.json", "r") as f:
#    data = json.load(f)

#X = np.array(data)
#print("Shape:", X.shape)



# ---------------------------------------------------
# Load scaler
# ---------------------------------------------------
x_scaler = joblib.load("x_scaler.pkl")

X_scaled = x_scaler.transform(X)

# ---------------------------------------------------
# Convert to tensor
# ---------------------------------------------------
X_t = torch.tensor(X_scaled, dtype=torch.float32)
#X_t = torch.tensor(X, dtype=torch.float32)

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
# Create model
# ---------------------------------------------------
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# output_dim must match training
output_dim =  y_dim = 90

model = MLPRegressor(
    input_dim=X.shape[1],
    output_dim=output_dim
).to(device)

# ---------------------------------------------------
# Load weights
# ---------------------------------------------------
model.load_state_dict(
    torch.load("mlp_regressor_50k.pth", map_location=device)
)

model.eval()

# ---------------------------------------------------
# Extract embeddings
# ---------------------------------------------------
with torch.no_grad():

    extracted_embeddings = model(
        X_t.to(device),
        return_embedding=True
    ).cpu().numpy()

print(extracted_embeddings.shape)

# ---------------------------------------------------
# Save embeddings
# ---------------------------------------------------
np.save(
    "local_nardini_FT.npy",
    extracted_embeddings
)

print("Saved embeddings:")
print("mlp_extracted_embeddings.npy")
