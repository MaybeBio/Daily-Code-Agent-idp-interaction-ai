import numpy as np
import sys, os
import re
import torch
from transformers import T5Tokenizer, T5EncoderModel


dir_name = os.path.dirname(__file__)
sys.path.append(dir_name)


def read_fasta(file_path):
    with open(file_path, 'r') as f:
        data = f.read()
    
    dic = {}
    for block in data.strip(">").split("\n>"):
        ac, seq = block.strip().split("\n", 1)
        ac = ac.strip()
        seq = seq.replace("\n", "").replace(" ", "")
        seq = re.sub(r"[UBJOZ]", "X", seq)
        dic[ac] = seq
    return dic


class ProtT5Embedder:
    def __init__(
        self,
        model_name="Rostlab/prot_t5_xl_half_uniref50-enc",
        device="cuda",
        chunk_size=2048,
        overlap=512,
        return_cpu=False,
    ):
        assert overlap < chunk_size
        assert overlap % 2 == 0, "overlap must be even"

        self.device = torch.device(device)
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.trim = overlap // 2
        self.return_cpu = return_cpu

        self.tokenizer = T5Tokenizer.from_pretrained(
            model_name, 
            do_lower_case=False,
        )

        self.model = (
            T5EncoderModel.from_pretrained(
                model_name,
            )
            .half()
            .to(self.device)
            .eval()
        )

        if self.device.type == "cpu":
            self.model = self.model.float()

    def embed(self, seq: str) -> torch.Tensor:
        seq = re.sub(r"[UZOB]", "X", seq)
        aa = list(seq)
        L = len(aa)

        step = self.chunk_size - self.overlap
        out_chunks = []

        with torch.inference_mode():
            start = 0
            idx = 0
            while start < L:
                end = min(start + self.chunk_size, L)
                chunk = " ".join(aa[start:end])

                tok = self.tokenizer(
                    chunk,
                    add_special_tokens=True,
                    return_tensors="pt",
                    padding=False,
                )
                input_ids = tok["input_ids"].to(self.device, non_blocking=True)
                attention_mask = tok["attention_mask"].to(self.device, non_blocking=True)

                hidden = self.model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                ).last_hidden_state

                residue = hidden[0, : (end - start), :]

                if idx == 0:
                    take_s = 0
                    take_e = (end - start) - self.trim if end < L else (end - start)
                else:
                    take_s = self.trim
                    take_e = (end - start) - self.trim if end < L else (end - start)

                residue = residue[take_s:take_e]

                if self.return_cpu:
                    residue = residue.detach().cpu()
                else:
                    residue = residue.detach()

                out_chunks.append(residue)

                del hidden, input_ids, attention_mask, tok

                if end == L:
                    break

                start += step
                idx += 1

        emb = torch.cat(out_chunks, dim=0)
        return emb
    
    def __call__(self, seq: str) -> torch.Tensor:
        return self.embed(seq)


class SequenceFeaturizer:
    def __init__(self, feature_type, device):
        self.feature_type = feature_type
        self.device = device
        self.embedding_method = ProtT5Embedder(device=self.device, chunk_size=2048, overlap=512)

    def featurize(self, seq):
        if self.feature_type == "ProtT5":
            return self.embedding_method(seq)
        else:
            raise ValueError(f"Unknown feature_type: {self.feature_type}")
 

class ResultWriter:
    def __init__(self,output_path,threthold=0.26):
        self.path = output_path
        with open(self.path,"w") as f:
            f.write(f"# Line1. accession\n")
            f.write(f"# Line2. sequence\n")
            f.write(f"# Line3. binding probability (float)\n")
            f.write(f"# Line4. binding prediction (1:binding, 0:non-binding), threshold={threthold}\n")
            f.write(f"# Line5. pLDDT (float)\n")
            f.write(f"# Line6. RSA (float)\n")
            f.write(f"# Line7. SS3 (H:Helix, S:Strand, C:Coil)\n")
            f.write(f"# Line8. SS8 (H:Alpha helix, G:3-10 helix, I:pi helix, E:Extended strand, B:Isolated beta bridge, T:Turn, S:Bend, -:Other (coil))\n")

    def append_write(self,ac,result,threthold=0.26):
        with open(self.path,"a") as f:
            f.write(f">{ac}\n")
            f.write(f"{result['seq']}\n")
            f.write(f"{','.join([f'{x:.4f}' for x in result['binding']])}\n")
            f.write(f"{','.join([str(x) for x in result['binary_binding']])}\n")
            f.write(f"{','.join([f'{x:.4f}' for x in result['plddt']])}\n")
            f.write(f"{','.join([f'{x:.4f}' for x in result['rsa']])}\n")
            f.write(f"{','.join(result['ss3_str'])}\n")
            f.write(f"{','.join(result['ss8_str'])}\n")

