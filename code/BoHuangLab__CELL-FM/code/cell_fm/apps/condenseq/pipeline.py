"""CELL-FM CondenSeq: sequence -> condensate titration curve.

The demo runs the same three stages as the offline analysis, in memory:

  1. generate   CELL-FM CS samples one image per protein-intensity level, over a
                log-spaced concentration ladder (evaluate_seq2img_dict.py)
  2. classify   a ViT labels each image condensed / diffuse (evaluate_single_img.py)
  3. summarise  a centred moving average over the binary calls gives the condensate
                probability curve (ma_plot.py), and AUC / AAC integrate it
                (analysis/ana_all_mutation_log_scale.py)

Weights come from a HF model repo by default; set CELLFM_LOCAL_WEIGHTS to a
directory holding the three .bin files to run against a local copy instead.
"""

import os
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd
import torch

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.models.cell_fm.cell_fm_cs_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_cs_model import CELLFMCSModel
from cell_fm.models.vit_cls_condenseq_img.config import ViTConfig
from cell_fm.models.vit_cls_condenseq_img.model import ViTModel

from . import metrics

MODEL_REPO = os.environ.get("CELLFM_MODEL_REPO", "BoHuangLab/CELL-FM")
LOCAL_WEIGHTS = os.environ.get("CELLFM_LOCAL_WEIGHTS", "")

# The weight repo holds one subfolder per model family, so paths are prefixed.
# CELLFM_LOCAL_WEIGHTS points at a flat directory, hence the basename() below.
GENERATOR_WEIGHTS = "condenseq/cellfm_seq2img.bin"
VAE_WEIGHTS = "condenseq/vae.bin"
CLASSIFIER_WEIGHTS = "condenseq/vit_cls.bin"
REFERENCE_NUCLEUS = "condenseq/reference_nucleus.npy"

# The concentration ladder. The offline runs scanned 4096 levels; the demo scans
# fewer for latency, and the moving-average window scales with it (see WINDOW_DIVISOR).
MIN_INTENSITY = 20.0
MAX_INTENSITY = 2048.0
WINDOW_DIVISOR = 8  # window = num_images // 8, matching ana_all_mutation_log_scale.py

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# CondenSeq is a fixed-length library: all 14,578 training sequences are exactly
# 66 residues. Anything else is off-distribution, so the demo requires this length
# rather than silently padding or cropping.
SEQUENCE_LENGTH = 66

# Generation starts from a random latent, so the curve moves a little between runs.
# The offline scripts fixed this at 6; keeping the same default makes the demo's
# numbers repeatable and comparable to them.
DEFAULT_SEED = 6

# Fixed sampling settings. The concentration ladder length and ODE step count match
# the offline runs, so demo numbers stay comparable to the published ones.
NUM_IMAGES = 512
NUM_STEPS = 100
BATCH_SIZE = 64


@lru_cache(maxsize=1)
def count_download() -> None:
    """Fetch the model repo's config.json, once per process.

    The Hub counts this repo's downloads by requests for that file, not for the .bin
    files, so this is what puts a Space start or a notebook run on the model page. A
    mirror set through CELLFM_MODEL_REPO may not carry the file; that is not an error.
    """
    from huggingface_hub import hf_hub_download
    from huggingface_hub.utils import EntryNotFoundError

    try:
        hf_hub_download(repo_id=MODEL_REPO, filename="config.json")
    except EntryNotFoundError:
        pass


def resolve_weights(filename: str) -> str:
    """Local directory if CELLFM_LOCAL_WEIGHTS is set, else pull from the model repo."""
    if LOCAL_WEIGHTS:
        path = os.path.join(LOCAL_WEIGHTS, os.path.basename(filename))
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found (CELLFM_LOCAL_WEIGHTS={LOCAL_WEIGHTS})")
        return path

    from huggingface_hub import hf_hub_download

    count_download()
    return hf_hub_download(repo_id=MODEL_REPO, filename=filename)


def generator_config() -> CELLFMConfig:
    """Hyperparameters from scripts/cell_fm_cs/evaluate_seq2img.sh."""
    return CELLFMConfig(
        img_resize=160,
        img_crop_size=160,
        img_type="GFP",
        cell_image="nucl",
        seq_zero_mask_ratio=0.0,
        path_type="Linear",
        prediction="velocity",
        # VAE
        num_down_blocks=3,
        latent_channels=4,
        vae_block_out_channels="128,256,512",
        # CELL-FM
        img_mask_ratio=0,
        cond_out_channels="32,64",
        sample_size=40,
        esm_embedding="esmc_600m",
        encoder_hidden_size=1152,
        max_protein_sequence_len=2048,
        encoder_num_hidden_layers=4,
        num_heads=8,
        dim_head=64,
        dropout=0,
        final_dropout=0,
        encoder_patch_size=4,
        # image generator
        img_generator_num_layers=4,
        img_generator_patch_size=2,
        attention_head_dim=64,
        num_attention_heads=8,
        # image decoder
        img_decoder_num_hidden_layers=2,
        img_decoder_hidden_size=512,
        img_decoder_num_heads=8,
        img_decoder_dim_head=64,
        cell_image_ratio=0.5,
        # checkpoints
        vae_loadcheck_path=resolve_weights(VAE_WEIGHTS),
        loadcheck_path=resolve_weights(GENERATOR_WEIGHTS),
        infer=True,
    )


def classifier_config() -> ViTConfig:
    """Hyperparameters from scripts/vit_cls_condenseq_img/pretrain.sh."""
    return ViTConfig(
        image_size=160,
        in_channels=2,
        patch_size=4,
        num_layers=8,
        num_heads=8,
        hidden_dim=512,
        mlp_dim=2048,
        num_classes=2,
        loadcheck_path=resolve_weights(CLASSIFIER_WEIGHTS),
        infer=True,
    )


def get_device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


@lru_cache(maxsize=1)
def load_models():
    """Download the weights and build both models, on CPU.

    Kept off the GPU deliberately: under ZeroGPU a CUDA device only exists inside a
    @spaces.GPU call, and the 2.1 GB download plus model construction would eat that
    budget. Call this once at startup, then models_on(device) inside the GPU call.
    """
    generator = CELLFMCSModel(config=generator_config(), loss_fn=UniDiffCriterions)
    generator.eval()

    classifier = ViTModel(config=classifier_config())
    classifier.eval()

    from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

    return generator, classifier, EsmSequenceTokenizer()


def models_on(device):
    """Move the cached models to device. A no-op once they are already there."""
    generator, classifier, vocab = load_models()
    return generator.to(device), classifier.to(device), vocab


@lru_cache(maxsize=1)
def reference_nucleus() -> torch.Tensor:
    """The DAPI channel every generated image is conditioned on, as (1, 1, 160, 160) in [-1, 1].

    Taken from CondenSeq protein index 12626, image 0 -- the same one the offline
    seq2img runs used, so demo curves are comparable to the published ones.
    """
    arr = np.load(resolve_weights(REFERENCE_NUCLEUS))
    return torch.from_numpy(arr).float().unsqueeze(0)


def intensity_ladder(num_images: int) -> torch.Tensor:
    """Log-spaced protein intensity levels, the x axis of the titration curve."""
    return torch.logspace(
        start=float(np.log10(MIN_INTENSITY)),
        end=float(np.log10(MAX_INTENSITY)),
        steps=num_images,
    )


def clean_sequence(sequence: str) -> str:
    """Strip FASTA headers, whitespace and case; reject anything not a standard residue."""
    lines = [ln for ln in str(sequence).splitlines() if not ln.startswith(">")]
    seq = "".join(lines).strip().upper()
    seq = "".join(seq.split())

    if not seq:
        raise ValueError("Enter a protein sequence.")
    bad = sorted(set(seq) - VALID_AA)
    if bad:
        raise ValueError(
            f"Sequence contains non-standard residues: {', '.join(bad)}. "
            "Use the 20 standard amino acids."
        )
    if len(seq) != SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence is {len(seq)} residues; CELL-FM CondenSeq takes exactly "
            f"{SEQUENCE_LENGTH}. The whole training library is {SEQUENCE_LENGTH}-mers, "
            "so other lengths are outside what the model has seen."
        )
    return seq


@dataclass
class Result:
    sequence: str
    intensities: np.ndarray          # (N,) protein intensity levels
    predictions: np.ndarray          # (N,) 0 = diffuse, 1 = condensed
    curve: np.ndarray                # (N,) moving-averaged condensate probability
    window: int
    auc: float                       # area under the curve
    aac: float                       # area above the curve (reentrant dissolution)
    no_reentrant_auc: float
    c_sat: float
    images: np.ndarray = field(repr=False, default=None)   # (N, 160, 160) every generated image
    image_intensities: np.ndarray = field(repr=False, default=None)


@torch.no_grad()
def generate_images(sequence, num_images, num_steps, batch_size, seed=DEFAULT_SEED, progress=None):
    """Stage 1: one generated image per intensity level, returned as uint16 (N, 2, 160, 160).

    uint16 is not incidental -- the offline pipeline wrote a TIFF between generation
    and classification, so the classifier saw quantised pixels. Keeping the same
    round-trip here keeps demo numbers comparable to the published ones.
    """
    device = get_device()
    generator, _, vocab = models_on(device)

    from esm.utils import encoding

    if seed is not None:
        torch.manual_seed(int(seed))
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(int(seed))

    tokens = encoding.tokenize_sequence(sequence, vocab, True).unsqueeze(0).to(device)
    nucleus = reference_nucleus().to(device)
    levels = intensity_ladder(num_images).to(device)

    out = []
    for start in range(0, num_images, batch_size):
        end = min(start + batch_size, num_images)
        n = end - start

        if progress is not None:
            progress(start / num_images, desc=f"Generating images {start}-{end} of {num_images}")

        sample = generator.sequence_to_image(
            tokens.repeat(n, 1),
            nucleus.repeat(n, 1, 1, 1),
            levels[start:end].unsqueeze(1),
            num_steps=num_steps,
        )

        # channel 0 = nucleus (conditioning), channel 1 = generated protein
        combined = torch.cat([nucleus.repeat(n, 1, 1, 1), sample], dim=1)
        combined = ((combined + 1) / 2.0).clamp(0.0, 1.0)
        out.append((combined * 65535).cpu().numpy().astype(np.uint16))

    return np.concatenate(out, axis=0)


@torch.no_grad()
def classify_images(images, batch_size, progress=None):
    """Stage 2: condensed (1) vs diffuse (0) for each generated image."""
    device = get_device()
    _, classifier, _ = models_on(device)

    preds = []
    for start in range(0, images.shape[0], batch_size):
        end = min(start + batch_size, images.shape[0])

        if progress is not None:
            progress(start / images.shape[0], desc=f"Classifying {start}-{end} of {images.shape[0]}")

        batch = torch.from_numpy(images[start:end].astype(np.float32)).to(device)
        batch = batch / 65535.0 * 2.0 - 1.0        # back to [-1, 1]
        logits = classifier.predict(batch)
        preds.append(torch.argmax(logits, dim=1).cpu().numpy())

    return np.concatenate(preds).astype(int)


def summarise(intensities, predictions, num_images):
    """Stage 3: moving average, then AUC / AAC / c_sat over the smoothed curve."""
    df = pd.DataFrame(
        {"protein_intensity_level": intensities, "predicted_class": predictions}
    )
    window = max(1, num_images // WINDOW_DIVISOR)
    kw = dict(window_size=window, use_log_scale=True)

    return df, window, dict(
        curve=metrics.moving_average(predictions, window),
        auc=metrics.compute_area_under_curve(df, **kw),
        aac=metrics.compute_area_above_curve(df, **kw),
        no_reentrant_auc=metrics.compute_area_if_no_reentrant(df, **kw),
        c_sat=metrics.compute_concentration(df, window_size=window),
    )


def sample_panel(images, intensities, k=6):
    """Pick k images spread evenly along the ladder, for the strip shown under the curve."""
    idx = np.linspace(0, len(intensities) - 1, k).astype(int)
    return images[idx, 1].astype(np.float32) / 65535.0, intensities[idx]


def run(sequence, num_images=NUM_IMAGES, num_steps=NUM_STEPS, batch_size=BATCH_SIZE,
        seed=DEFAULT_SEED, progress=None) -> Result:
    """The whole pipeline: sequence in, titration curve and its areas out."""
    sequence = clean_sequence(sequence)

    images = generate_images(sequence, num_images, num_steps, batch_size, seed, progress)
    predictions = classify_images(images, batch_size, progress)

    intensities = intensity_ladder(num_images).numpy()
    _, window, stats = summarise(intensities, predictions, num_images)
    # every image, in ladder order, so the gallery shows all of them
    panel, panel_x = sample_panel(images, intensities, k=num_images)

    return Result(
        sequence=sequence,
        intensities=intensities,
        predictions=predictions,
        window=window,
        images=panel,
        image_intensities=panel_x,
        **stats,
    )
