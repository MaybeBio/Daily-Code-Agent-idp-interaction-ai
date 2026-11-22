import torch
from tqdm import tqdm


SQN_EMBED_MODEL_PARTS = ('embed_model.model', 'embed_model1.model', 'embed_model2.model')


def get_partial_state_dict(model):
    return {
        k: v
        for k, v in model.state_dict().items()
        if not k.startswith(SQN_EMBED_MODEL_PARTS)
    }


def save_partial_model(model, checkpoint_path):
    torch.save({'model': get_partial_state_dict(model)}, checkpoint_path)


# Training Loops
def train(model, dataloader, optimizer, criterion, epochs, checkpoint_path, device, verbose=True):
    losses_train = []
    for epoch in range(epochs):
        if verbose:
            print('-'*10)
            print(f'Epoch {epoch + 1}/{epochs}:')
        model.train()
        running_loss = 0.0
        for batch in tqdm(dataloader):
            x = batch[:-1]
            labels = batch[-1]
            for x0 in x:
                if not isinstance(x0[0], str):
                    x0 = x0.to(device).float()
            labels = labels.to(device).float()
            model.zero_grad()
            logits = model(*x)
            loss = criterion(logits, labels)
            loss.backward()
            running_loss += loss.item()
            optimizer.step()
        
        epoch_loss = running_loss / len(dataloader)
        losses_train.append(epoch_loss)
        if verbose:
            print(f"Loss: {epoch_loss:.2f}")
    
    torch.save(
        {'epoch': epoch,
         'model': get_partial_state_dict(model),
         'optimizer': optimizer.state_dict(),
         'loss': epoch_loss},
        checkpoint_path
    )
    return losses_train
