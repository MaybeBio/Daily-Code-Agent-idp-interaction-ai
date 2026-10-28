import wandb
import os
from pathlib import Path
import numpy as np
import time
from datasets import Dataset, DatasetDict
from typing import Dict, List
import argparse

from idr_diverge.utils.helpers import resolve_config_paths

from transformers import (
    BertConfig,
    BertForMaskedLM,
    BertTokenizer,
    DataCollatorForLanguageModeling,
    Trainer,
    TrainingArguments,
)

from idr_diverge.IDR_LM_train.config import load_config

# Portions of this code are adapted from Hugging Face Transformers documentation and example scripts for language model pretraining. 
# Structural modifications were made with assistance from ChatGPT (OpenAI) and further refined by Denise Le.


#DEFAULT_CONFIG_PATH = './pretrain_args.yaml'

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run script using a YAML config file."
    )

    parser.add_argument(
        "--config",
        type=str,
        required=True,
        #default=DEFAULT_CONFIG_PATH,
        help="Path to YAML config file for pretraining arguments."
    )

    args = parser.parse_args()

    # Convert to absolute path and check it exists
    args.config = Path(args.config).resolve()

    if not args.config.exists():
        parser.error(f"Config file does not exist: {args.config}")

    return args

CONFIG_PATH = parse_args().config

args = load_config(CONFIG_PATH)

args = resolve_config_paths(config=args, config_path=CONFIG_PATH,
                                  path_keys={'train_data_path', 'model_config_path', 'resume_checkpoint_path', 'final_model_dir'})


# ---------------------------
# Helpers
# ---------------------------

def read_fasta(file_path) -> List[Dict]:
    sequences = []
    sequence_id = None
    sequence_data = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id is not None:
                    sequences.append({"id": sequence_id, "sequences": ''.join(sequence_data)})
                sequence_id = line[1:]  # Remove the '>' character
                sequence_data = []
            else:
                sequence_data.append(line)

        # Add the last sequence
        if sequence_id is not None:
            sequences.append({"id": sequence_id, "sequences": ''.join(sequence_data)})

    return sequences


def combine_sequences(sequences, tokenizer, max_length=512):
    """
    Combines multiple sequences into a list of tokenized sequences with a specified maximum length.
    
    Args:
        sequences (list of str): List of input sequences to be tokenized and combined.
        tokenizer (PreTrainedTokenizer): The tokenizer to use for encoding the sequences.
        max_length (int): The maximum length of the combined token sequences.
    
    Returns:
        list of list of int: List of tokenized sequences, each within the specified maximum length.
    """
    combined_sequences = []  # List to store combined tokenized sequences
    combined_tokens = [tokenizer.cls_token_id]  # Start with the CLS token
    combined_attention = []
    current_length = 1  # Account for the CLS token

    for seq in sequences:
        seq = ' '.join(seq)
        tokens = tokenizer.encode(seq, add_special_tokens=False, truncation=True, max_length=max_length-2)  # Tokenize without special tokens
        tokens.append(tokenizer.sep_token_id)  # Append SEP token

        if current_length + len(tokens) > max_length:
            # If adding the current tokens exceeds max_length, pad and save the current sequence
            remaining_space = max_length - current_length
            combined_tokens.extend([tokenizer.pad_token_id] * remaining_space)  # Pad the remaining space
            #combined_sequences.append(combined_tokens)
            
            attention_mask = [1] * current_length + [0] * remaining_space
            combined_attention.append(attention_mask)
            combined_sequences.append(combined_tokens)

            # Reset for the next sequence
            combined_tokens = [tokenizer.cls_token_id] + tokens
            current_length = len(combined_tokens)
        else:
            combined_tokens.extend(tokens)
            current_length += len(tokens)

    # Append the last combined sequence
    remaining_space = max_length - current_length
    combined_tokens.extend([tokenizer.pad_token_id] * remaining_space)  # Pad the remaining space

    attention_mask = [1] * current_length + [0] * remaining_space

    combined_sequences.append(combined_tokens)
    combined_attention.append(attention_mask)
    combined_dict = {
            'input_ids': combined_sequences,
            'attention_mask': combined_attention
            
    }
    return combined_dict

def tokenize_dataset(dataset:DatasetDict,tokenizer:BertTokenizer,max_length=512) -> DatasetDict:
    #Appends sequences together to fill up the context length
    #Assumes test/train split
    
    train_input = [' '.join(seq) for seq in dataset['train']['sequences']]
    test_input = [' '.join(seq) for seq in dataset['test']['sequences']]
    
    train_tokens  = combine_sequences(train_input, tokenizer, max_length)
    test_tokens  = combine_sequences(test_input, tokenizer, max_length)
    
    train_tokenized_dataset = Dataset.from_dict(train_tokens)
    test_tokenized_dataset = Dataset.from_dict(test_tokens)
    
    tokenized_dataset = DatasetDict({
    'train': train_tokenized_dataset,
    'test': test_tokenized_dataset
    })
    
    return tokenized_dataset

def load_fasta_as_dataset(file_path):
    sequences = read_fasta(file_path)
    return Dataset.from_list(sequences)

def initialize_wandb():
    os.environ["WANDB_LOG_MODEL"] = "checkpoint"
    wandb.init(project=args.wandb_project, name=args.wandb_run_name)

def load_and_tokenize_dataset(data_file: Path, tokenizer: BertTokenizer) -> DatasetDict:
    dataset = load_fasta_as_dataset(str(data_file))
    dataset = dataset.shuffle(seed=args.seed)#.select(range(10)) #TODO comment this out later
    dataset = dataset.train_test_split(test_size=args.test_size)

    print("Tokenizing dataset...") 
    tokenized_dataset = tokenize_dataset(dataset, tokenizer, max_length=args.max_seq_length)
    
    return tokenized_dataset

def load_model(config_file: Path, checkpoint_path: Path | None = None) -> BertForMaskedLM:
    config = BertConfig.from_json_file(str(config_file))

    if checkpoint_path is not None and checkpoint_path.exists():
        print(f"Loading model from checkpoint: {checkpoint_path}")
        model = BertForMaskedLM.from_pretrained(str(checkpoint_path), config=config)
    else:
        print("Initializing model from config") # Initialize the model from a config without pretrained weights
        model = BertForMaskedLM(config=config)

    return model

def build_training_args(logging_steps: int) -> TrainingArguments:
    return TrainingArguments(
        output_dir=str(args.checkpoint_output_dir),
        overwrite_output_dir=True,
        do_train=True,
        evaluation_strategy="epoch",
        save_strategy="epoch",
        num_train_epochs=args.num_train_epochs,
        learning_rate=args.learning_rate,
        no_cuda=args.no_cuda,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        save_total_limit=10,
        logging_dir= args.checkpoint_output_dir / "logs",
        logging_steps=logging_steps,
        warmup_ratio=args.warmup_proportion,
        
        load_best_model_at_end=True, 
        report_to="wandb"
    )


def validate_tokenized_dataset(tokenized_dataset):
    assert len(tokenized_dataset["train"][0]["input_ids"]) == args.max_seq_length
    assert all(
        len(input_ids) <= args.max_seq_length for input_ids in tokenized_dataset["train"]["input_ids"]
    ), "Input sequence exceeds max length"


def main():
    start = time.time()
    print(f"Start time: {start}")

    initialize_wandb()
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert")
    tokenized_dataset = load_and_tokenize_dataset(args.train_data_path , tokenizer)

    train_dataset = tokenized_dataset["train"]
    eval_dataset = tokenized_dataset["test"]

    validate_tokenized_dataset(tokenized_dataset)

    print(f"Dataset prep time: {time.time() - start:.2f} seconds")

    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer,
        mlm=True,
        mlm_probability=args.masked_lm_prob,
    )

    model = load_model(args.model_config_path, args.resume_checkpoint_path)

    logging_steps = max(1, len(train_dataset) // args.batch_size)
    training_args = build_training_args(logging_steps=logging_steps)

    trainer = Trainer(
        model=model,
        args=training_args,
        data_collator=data_collator,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset
    )

    trainer.train()

    eval_results = trainer.evaluate()
    print("Evaluation results:")
    print(eval_results)

    trainer.save_model(str(args.final_model_dir))

    end = time.time()
    print(f"Total runtime: {end - start:.2f} seconds")


if __name__ == "__main__":
    main()
