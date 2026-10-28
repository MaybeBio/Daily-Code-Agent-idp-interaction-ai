from pathlib import Path
import argparse
import subprocess
import yaml
from idr_diverge.utils.helpers import resolve_config_paths



#Command-line interface for generating embeddings pipeline.
#Handles argument parsing and loading config file, and makes calls to embedding scripts in idr_diverge/embed
#This script was developed with assistance from ChatGPT (OpenAI) and has been reviewed and modified by Denise Le.



def parse_args():
    # First parser only reads --config
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument(
        "--config",
        type=Path,
        help="Path to YAML config file."
    )

    config_args, remaining_argv = config_parser.parse_known_args()

    # Load config values if provided
    config = {}
    if config_args.config:
        with open(config_args.config) as f:
            config = yaml.safe_load(f) or {}
        
        config = resolve_config_paths(
            config,
            config_args.config,
            path_keys={
                "esm_script",
                "idr_lm_script",
                "bert_config",
                "pretrain_config",
                "model_file",
                "data_dir",
                "result_dir",
            },
        )
        
        #print(config)

    parser = argparse.ArgumentParser(
        description="Run embedding scripts for protein sequence files.",
        parents=[config_parser]
    )
    
    parser.add_argument(
        "--embed-type",
        choices=["esm", "IDR_LM", "IDR_LM_random"],
        default=config.get("embed_type"),
        help="Embedding type to generate.",
    )
    
    parser.add_argument(
        "--esm-script",
        default=config.get("esm_script"),
        help="Script to get esm embeddings.",
    )
    
    parser.add_argument(
        "--idr-lm-script",
        default=config.get("idr_lm_script"),
        help="Script to get esm embeddings.",
    )

    parser.add_argument(
        "--model-file",
        type=Path,
        default=config.get("model_file"),
        help="Path to the trained model checkpoint. Required for IDR_LM.",
    )
    
    parser.add_argument(
        "--pretrain-config",
        type=Path,
        default=config.get("pretrain_config"),
        help="Path to the pretraining configuration file. Required for IDR_LM.",
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        default=config.get("data_dir"),
        help="Directory containing input FASTA files or path to a single FASTA file.",
    )

    parser.add_argument(
        "--result-dir",
        type=Path,
        default=config.get("result_dir"),
        help="Directory where embedding outputs will be written.",
    )

    args = parser.parse_args(remaining_argv)

    # Validation
    if args.embed_type is None:
        parser.error("--embed-type must be provided either in the config or on the command line.")

    if args.embed_type == "IDR_LM" and args.model_file is None:
        parser.error("--model-file is required when --embed-type is IDR_LM.")

    return args


args = parse_args()
ESM_SCRIPT = args.esm_script
IDR_LM_SCRIPT = args.idr_lm_script
PRETRAIN_CONFIG = args.pretrain_config


def get_input_files(data_dir: Path) -> list[Path]:
    """
    Return sequence files from the input directory.
    """
    if not data_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {data_dir}")
    if data_dir.is_file():
        if not data_dir.suffix.lower() in {".fa", ".fasta", ".faa", ".txt"}:
            raise ValueError(f"Input path must be a fasta file or directory with fasta files")
        else:
            files = [data_dir]
            return files

    files = sorted(
        f for f in data_dir.iterdir()
        if f.is_file() and f.suffix.lower() in {".fa", ".fasta", ".faa", ".txt"}
    )

    if not files:
        raise FileNotFoundError(f"No FASTA-like files found in: {data_dir}")

    return files


def ensure_output_dir(result_dir: Path) -> None:
    """
    Create the output directory if it does not already exist.
    """
    result_dir.mkdir(parents=True, exist_ok=True)


def extract_model_tag(model_file: Path) -> str:
    """
    Create a short tag from the model filename for output naming.
    Example:
        checkpoint_epoch_10.pt -> pt
        model_lowlr_epoch10 -> epoch10
    """
    stem = model_file.stem
    parts = stem.split("_")
    return parts[-1] if parts else stem


def run_command(cmd: list[str]) -> None:
    """
    Run a subprocess command and print debugging output if it fails.
    """
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("Command failed:")
        print(" ".join(cmd))
        if result.stdout:
            print("\nSTDOUT:")
            print(result.stdout)
        if result.stderr:
            print("\nSTDERR:")
            print(result.stderr)
        raise subprocess.CalledProcessError(
            returncode=result.returncode,
            cmd=cmd,
            output=result.stdout,
            stderr=result.stderr,
        )

    if result.stdout:
        print(result.stdout.strip())


def run_esm_embeddings(input_files: list[Path], result_dir: Path) -> None:
    """
    Run the ESM embedding script on each input file.
    """
    for input_file in input_files:
        output_prefix = result_dir / f"{input_file.stem}_esm"
        print(f"Processing {input_file.name} -> {output_prefix}")

        cmd = [
            "python",
            str(ESM_SCRIPT),
            str(input_file),
            str(output_prefix),
        ]
        run_command(cmd)


def run_idr_lm_embeddings(
    input_files: list[Path],
    result_dir: Path,
    model_file: Path,
    random_init: bool = False,
) -> None:
    """
    Run the IDR_LM embedding script on each input file.

    If random_init=True, call the script in its random/uninitialized mode.
    """
    if not random_init and model_file is None:
        raise ValueError(
            "--model-file is required when --embed-type is 'IDR_LM'."
        )

    if model_file is not None and not model_file.exists():
        raise FileNotFoundError(f"Model file does not exist: {model_file}")

    model_tag = extract_model_tag(model_file) if model_file else "random"

    for input_file in input_files:
        output_name = f"{input_file.stem}_{model_tag}_IDRLM"
        
        if random_init:
            output_name = f"randominit_{input_file.stem}_{model_tag}_IDRLM"
            
            cmd = ["python",str(IDR_LM_SCRIPT),
            "--seq_input_path", str(input_file),
            "--out-dir", str(result_dir),
            "--pretrain-config", str(PRETRAIN_CONFIG),
            "--output-name", output_name]
         

        else:
            cmd = ["python",str(IDR_LM_SCRIPT),
                    "--seq_input_path",str(input_file),
                    "--model-path", str(model_file),
                    "--pretrain-config", str(PRETRAIN_CONFIG),
                    "--out-dir", str(result_dir),
                    "--output-name", output_name]
        output_prefix = result_dir / output_name
        print(f"Processing {input_file.name} -> {output_prefix}")
        print(cmd)
        run_command(cmd)


def main():
    args = parse_args()
    
    ensure_output_dir(args.result_dir)
    input_files = get_input_files(args.data_dir)

    if args.embed_type == "esm":
        run_esm_embeddings(input_files, args.result_dir)

    elif args.embed_type == "IDR_LM":
        run_idr_lm_embeddings(
            input_files=input_files,
            result_dir=args.result_dir,
            model_file=args.model_file,
            random_init=False,
        )

    elif args.embed_type == "IDR_LM_random":
        run_idr_lm_embeddings(
            input_files=input_files,
            result_dir=args.result_dir,
            model_file=None,
            random_init=True,
        )


if __name__ == "__main__":
    main()