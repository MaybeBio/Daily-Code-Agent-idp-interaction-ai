# torch
import torch
from torch.nn.functional import sigmoid


def get_llps_score(model, *args):
    with torch.no_grad():
        args = ([a] if isinstance(a, str) else a.unsqueeze(0) for a in args)
        logits = model(*args)
        score = sigmoid(logits)
    return score.item()


def get_batch_llps_score(model, *args):
    with torch.no_grad():
        logits = model(*args)
        score = sigmoid(logits)
    return score.detach().cpu().numpy()


if __name__ == '__main__':
    # LLPS models
    from src.lib.llps_model import get_llps_model
    # input
    from src.utils.force_field import force_field_transform_factory
    from src.utils.saprotseq import get_saprotseq

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    force_field_transform = force_field_transform_factory(img_size=320)
    # model
    model = get_llps_model('ensemble', local_files_only=True)
    model = model.to(device)
    # input
    # sequence: str
    # saprotseq: str
    # force_field: Tensor, shape (1, D, D)
    sequence = 'MEELSADEIRRRRLARLAGGQTSQPTTPLTSPQRENPPGPPIAASAPGPSQSLGLNVHNMTPATSPIGASGVAHRSQSSEGVSSLSSSPSNSLETQSQSLSRSQSMDIDGVSCEKSMSQVDVDSGIENMEVDENDRREKRSLSDKEPSSGPEVSEEQALQLVCKIFRVSWKDRDRDVIFLSSLSAQFKQNPKEVFSDFKDLIGQILMEVLMMSTQTRDENPFASLTATSQPIAAAARSPDRNLLLNTGSNPGTSPMFCSVASFGASSLSSLYESSPAPTPSFWSSVPVMGPSLASPSRAASQLAVPSTPLSPHSAASGTAAGSQPSSPRYRPYTVTHPWASSGVSILSSSPSPPALASSPQAVPASSSRQRPSSTGPPLPPASPSATSRRPSSLRISPSLGASGGASNWDSYSDHFTIETCKETDMLNYLIECFDRVGIEEKKAPKMCSQPAVSQLLSNIRSQCISHTALVLQGSLTQPRSLQQPSFLVPYMLCRNLPYGFIQELVRTTHQDEEVFKQIFIPILQGLALAAKECSLDSDYFKYPLMALGELCETKFGKTHPVCNLVASLRLWLPKSLSPGCGRELQRLSYLGAFFSFSVFAEDDVKVVEKYFSGPAITLENTRVVSQSLQHYLELGRQELFKILHSILLNGETREAALSYMAAVVNANMKKAQMQTDDRLVSTDGFMLNFLWVLQQLSTKIKLETVDPTYIFHPRCRITLPNDETRVNATMEDVNDWLTELYGDQPPFSEPKFPTECFFLTLHAHHLSILPSCRRYIRRLRAIRELNRTVEDLKNNESQWKDSPLATRHREMLKRCKTQLKKLVRCKACADAGLLDESFLRRCLNFYGLLIQLLLRILDPAYPDITLPLNSDVPKVFAALPEFYVEDVAEFLFFIVQYSPQALYEPCTQDIVMFLVVMLCNQNYIRNPYLVAKLVEVMFMTNPAVQPRTQKFFEMIENHPLSTKLLVPSLMKFYTDVEHTGATSEFYDKFTIRYHISTIFKSLWQNIAHHGTFMEEFNSGKQFVRYINMLINDTTFLLDESLESLKRIHEVQEEMKNKEQWDQLPRDQQQARQSQLAQDERVSRSYLALATETVDMFHILTKQVQKPFLRPELGPRLAAMLNFNLQQLCGPKCRDLKVENPEKYGFEPKKLLDQLTDIYLQLDCARFAKAIADDQRSYSKELFEEVISKMRKAGIKSTIAIEKFKLLAEKVEEIVAKNARAEIDYSDAPDEFRDPLMDTLMTDPVRLPSGTIMDRSIILRHLLNSPTDPFNRQTLTESMLEPVPELKEQIQAWMREKQNSDH'
    saprotseq = get_saprotseq(sequence)
    force_field = force_field_transform(sequence).to(device)
    # inference
    score = get_llps_score(model, sequence, saprotseq)
    print("Predicted score:", score)
    # batch inference, provide input based on the model type
    # sequence: list of str
    # saprotseq: list of str
    # force_field: Tensor, e.g., torch.stack([force_field, force_field])
    scores = get_batch_llps_score(
        model, [sequence, sequence], [saprotseq, saprotseq]
    )
    print("Predicted scores:", scores)
