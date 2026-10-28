import torch
import torch.nn as nn

class UniDiffCriterions(nn.Module):
    def __init__(self, args, reduction="none") -> None:
        super().__init__()
        self.sequence_loss = nn.CrossEntropyLoss(reduction=reduction)

        self.args = args
        self.seq_loss_coeff = self.args.seq_loss_coeff
        self.img_diff_loss_coeff = self.args.img_diff_loss_coeff
        self.img_recon_loss_coeff = self.args.img_recon_loss_coeff

    def forward(self, batch_data, model_output):
        with torch.no_grad():
            protein_seq_mask = batch_data["protein_seq_mask"].bool()
            protein_seq = batch_data["protein_seq"][protein_seq_mask]
            zm_label = batch_data['zm_label'].squeeze(-1).bool()

            mask_num = torch.sum(protein_seq_mask.float(), dim=1, keepdim=True)
            seq_loss_weights = torch.ones_like(protein_seq_mask, dtype=torch.float32)
            seq_loss_weights = seq_loss_weights / mask_num
            seq_loss_weights = seq_loss_weights * protein_seq_mask.float()
            seq_loss_weights = seq_loss_weights[protein_seq_mask.bool()]

        protein_seq_output, img_diff_loss, img_recon_loss = model_output

        # Protein sequence loss
        protein_seq_output = protein_seq_output[protein_seq_mask]
        sequence_loss = (
            self.sequence_loss(
                protein_seq_output.to(torch.float32).view(-1, protein_seq_output.size(-1)), 
                protein_seq.view(-1), 
            )
            * self.seq_loss_coeff
        )
        sequence_loss = sequence_loss * seq_loss_weights

        if zm_label.all():
            sequence_loss = torch.zeros(1, device=sequence_loss.device).sum()
        else:
            sequence_loss = sequence_loss.sum() / (~zm_label).float().sum()
        
        # Protein image diffusion loss
        img_diff_loss = img_diff_loss.mean() * self.img_diff_loss_coeff

        # Protein image reconstruction loss
        img_recon_loss = img_recon_loss.mean() * self.img_recon_loss_coeff

        loss = sequence_loss + img_diff_loss + img_recon_loss

        loss_items = {
            'sequence_loss': (sequence_loss, self.seq_loss_coeff),
            'image_diff_loss': (img_diff_loss, self.img_diff_loss_coeff),
            'image_recon_loss': (img_recon_loss, self.img_recon_loss_coeff),
        }

        log_loss = {'total_loss': loss.item()}
        for name, (val, coeff) in loss_items.items():
            if coeff > 0:
                log_loss[name] = val.item() / coeff

        return loss, log_loss