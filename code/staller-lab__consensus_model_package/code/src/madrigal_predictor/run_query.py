import torch
from parrot.brnn_architecture import BRNN_MtM
from importlib.resources import files, as_file

# Select your hardware device (GPU if available, otherwise CPU)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = BRNN_MtM(20, hidden_size=10, num_layers=1, num_classes=1, device=device)

weights = files("madrigal_predictor") / "consensus_model_all_six.pt"
with as_file(weights) as model_path:
    state = torch.load(model_path, map_location=device, weights_only=True) # Best practice: use weights_only=True for security when unpickling
model.load_state_dict(state)
model = model.to(device)

# Set the model to evaluation mode
model.eval()

# 5. Function That Takes Input and Returns Output
def run(query_tensor):
    with torch.no_grad():
        predictions = model(query_tensor.to(device))
    return predictions.cpu()


