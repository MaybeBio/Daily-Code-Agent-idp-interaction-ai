import torch
from tqdm import tqdm
from torch.nn.functional import sigmoid


# Test Loops
def predict(model, dataloader, criterion, device, verbose=True):
    model.eval()
    output = []
    running_loss = 0.0
    with torch.no_grad():
        for batch in tqdm(dataloader):
            x = batch[:-1]
            labels = batch[-1]
            for i in range(len(x)):
                if not isinstance(x[i][0], str):
                    x[i] = x[i].to(device).float()
            labels = labels.to(device).float()
            logits = model(*x)
            prob = sigmoid(logits)
            output.append(prob.detach())
            loss = criterion(logits, labels)
            running_loss += loss.item()
    output = torch.cat(output, dim=0)
    output = output.detach().cpu().numpy()
    loss_test = running_loss / len(dataloader)
    if verbose:
        print(f"Loss: {loss_test:.2f}")
    return output
