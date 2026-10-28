import os
import torch
import numpy as np
import pandas as pd

from torch.utils.data import DataLoader, SequentialSampler, TensorDataset
from transformers import BertConfig, BertModel, BertTokenizer, BertForMaskedLM

from idr_diverge.IDR_LM_train.config import load_config
from idr_diverge.utils.helpers import resolve_config_paths
import argparse
from pathlib import Path


#Core logic adapted from 
#Pang & Liu (2023) "IDP-LM: Prediction of protein intrinsic disorder and disorder functions based on language models"
#https://github.com/YihePang/IDP-LM
#Modified by Denise Le



os.environ['CUDA_VISIBLE_DEVICES'] = '0'

#MODEL_CONFIG_PATH = 
#OUT_DIR = '/home/moseslab/denise/embeddings/IDR_LM'

#DEFAULT_CONFIG = 'Paper/configs/idrlm_pretrain.yaml'
#args = load_config(DEFAULT_CONFIG)

def load_file_2_data(file_path):
    """Load FASTA file and return list of [id, sequence] pairs."""
    with open(file_path, 'r') as loadfile:
        load_f = [line.strip('\n') for line in loadfile]
    
    load_data = []
    for i in range(len(load_f)):
        if i % 2 == 0:
            load_data.append(load_f[i:i+2])
    
    return load_data


def convert_sequences_to_inputs(sequences,tokenizer,vocab_list, max_seq_length):
    """
    Convert protein sequences to model input format with tokenization and padding.
    """
    each_max = max_seq_length - 2
	
    new_data = []
    for i, sequence in enumerate(sequences):
        one_seq  = [r for r in sequence]
        s = 0
        
        for j in range(int(-(-len(one_seq)//each_max))):  #向上取整
            if s + each_max >= len(one_seq):
                end = len(one_seq) - s
                new_data.append(one_seq[s:s+end])
            elif s + each_max < len(one_seq):
                new_data.append(one_seq[s:s+each_max])
            s = s + each_max

    new_new_data = []
    for one_data in new_data:
        new_new_data.append(["[CLS]"] + one_data + ["[SEP]"])

    inputs = []
    for one_data in new_new_data:
        input_ids = tokenizer.convert_tokens_to_ids(one_data)

		# zero_padding
        input_array = np.zeros(max_seq_length, dtype=int)
        input_array[:len(input_ids)] = input_ids

		# mask
        mask_array = np.zeros(max_seq_length, dtype=int)
        mask_array[:len(input_ids)] = 1


        one_input = []
        one_input.append(input_array)
        one_input.append(mask_array)
        inputs.append(one_input)
	
	# print("covert finished, got totally:",len(inputs))

    return inputs

def resolve_path(config_file, rel_path):
    base = Path(config_file).parent  # directory where config lives
    return (base / rel_path).resolve()


def get_encoding_from_model(sequences, inputs, model_path, config_path,args,pretrain_config_path):
    """
    Generates embeddings from a pre-trained BERT model.

    Args:
        sequences (list): List of input sequences for which embeddings are to be generated.
        inputs (list): List of input tensors containing token IDs and attention masks for each sequence.
        model_path (str): Path to the pre-trained BERT model checkpoint.

    Returns:
        list: A list of embeddings for each input sequence, after removing padding and restoring original sequence lengths.
    """
    # Device configuration (use GPU if available)
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    
    if not(model_path):
        #Randomly initialized model
        #print(config_path)
        #print(Path.cwd())
        config_path = resolve_path(pretrain_config_path,config_path)
        #config_path = resolve_config_paths(args.model_config_path,config_path,path_keys={"bert_config"})
        print(config_path)
        config = BertConfig.from_json_file(config_path)
        model = BertForMaskedLM(config=config)
    
    else:
        # Load the pre-trained BERT model
        model = BertModel.from_pretrained(model_path)
    
    model.to(device)
    model.eval()

    model.config.output_hidden_states = True #MODIFIED

    # Convert inputs to tensors
    all_input_ids = torch.tensor([f[0] for f in inputs], dtype=torch.long)
    all_input_mask = torch.tensor([f[1] for f in inputs], dtype=torch.long)
    
    # Create DataLoader
    test_data = TensorDataset(all_input_ids, all_input_mask)
    test_sampler = SequentialSampler(test_data)
    test_dataloader = DataLoader(test_data, sampler=test_sampler, batch_size=args.batch_size)

    get_encodings = []

    # Encode sequences using BERT
    for input_ids, input_mask in test_dataloader:
        input_ids, input_mask = input_ids.to(device), input_mask.to(device)
        with torch.no_grad():
            #**********MODIFIED*************
            outputs = model(input_ids=input_ids, attention_mask=input_mask)
            last_layer_embeddings = outputs.hidden_states[-1]
            get_encodings.append(last_layer_embeddings.cpu().numpy())

            #encoded_outputs, pooled_output = model(input_ids=input_ids, attention_mask=input_mask)
            #get_encodings.append(encoded_outputs.cpu().numpy())
    # Method 2: Average the hidden states across all tokens to get a protein-level embedding
        
    # Combine batch encodings
    all_batch_encodings = np.concatenate(get_encodings, axis=0)

    org_lens_encoding = []
    each_max = args.max_seq_length - 2
    begin_idx = 0

    # Restore sequences from slices
    for seq in sequences:
        slice_num = int(np.ceil(len(seq) / each_max))
        if slice_num > 1:
            seq_encodings = all_batch_encodings[begin_idx:begin_idx + slice_num].reshape(-1, all_batch_encodings.shape[-1])
        else:
            seq_encodings = all_batch_encodings[begin_idx].reshape(-1, all_batch_encodings.shape[-1])
        begin_idx += slice_num
        org_lens_encoding.append(seq_encodings)

    final_encodings = []

    # Remove padding and restore original sequence lengths
    for i, seq in enumerate(sequences):
        slice_num = int(np.ceil(len(seq) / each_max))
        org_len = len(seq) + 2 * slice_num
        final_encodings.append(org_lens_encoding[i][:org_len])

    cleaned_encodings = []

    for i, seq in enumerate(sequences):
        slice_num = int(np.ceil(len(seq) / each_max))
        if slice_num <= 1:
            cleaned_encoding_1 = final_encodings[i][1:-1]
        else:
            del_idx = [j * args.max_seq_length for j in range(slice_num)] + \
                      [((j + 1) * args.max_seq_length) - 1 for j in range(slice_num - 1)] + \
                      [len(seq) + 2 * slice_num - 1]
            cleaned_encoding_1 = np.delete(final_encodings[i], del_idx, axis=0)
        
        mean_embedding = np.mean(cleaned_encoding_1, axis=0)
        cleaned_encodings.append(mean_embedding)

    return cleaned_encodings

    

def split_into_batches(lst, batch_size):
    for i in range(0, len(lst), batch_size):
        yield lst[i:i + batch_size]


def main(data, encoding_path, model_path,args,out_name='encodings', pretrain_config_path=None):
    """
    Processes dataset and generates embeddings, saving them to files.

    Args:
        data (list): List of tuples where each tuple contains a sequence name and the corresponding sequence data.
        encoding_path (str): Directory where the generated embeddings will be saved.
        model_path (str): Path to the pre-trained BERT model checkpoint.

    Returns:
        None
    """
    # Create output directory if it doesn't exist
    os.makedirs(encoding_path, exist_ok=True)
    batch_size =100

    csv_file_path = os.path.join(encoding_path, out_name+".csv")
    
    # Filter sequences that have not been encoded
    #for seq_name, seq_data in data:
    #    seq_name = seq_name.replace('>', '')
    #    if not os.path.exists(os.path.join(encoding_path, f"{seq_name}.npy")):
    #        inputs_sequences.append(seq_data)
    #        final_data.append((seq_name, seq_data))
    
    #if inputs_sequences:
    # Load vocabulary and tokenizer
    
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert")
    #vocab_list = [line.strip("\n") for line in open(args.vocab_file, 'r')]
    vocab_list = list(tokenizer.vocab.keys())

    #Split data into batches
    data_batches = list(split_into_batches(data, batch_size))

    for i, data_batch in enumerate(data_batches):
        print(str(i),'/',str(len(data_batches)))
        #print(data_batch)
        
        inputs_sequences = []
        final_data = []
        
        for seq_name, seq_data in data_batch:
            seq_name = seq_name.replace('>', '')
            print(seq_name)
            inputs_sequences.append(seq_data)
            final_data.append((seq_name, seq_data))
            
        # Convert sequences to model inputs
        inputs = convert_sequences_to_inputs(inputs_sequences, tokenizer, vocab_list, args.max_seq_length)

        # Get model encodings
        encodings = get_encoding_from_model(inputs_sequences, inputs, model_path, config_path = args.model_config_path, args=args, pretrain_config_path=pretrain_config_path)

            # Save encodings to files
            #for seq_name, encoding in zip([fd[0] for fd in final_data], encodings):
            #    np.save(os.path.join(encoding_path, f"{seq_name}.npy"), encoding)
            
            # Prepare data for CSV
        csv_data = []
        for seq_name, encoding in zip([fd[0] for fd in final_data], encodings):
            flattened_encoding = encoding.flatten()  # Flatten the encoding array
            csv_data.append([seq_name] + flattened_encoding.tolist())
     
        # Convert to DataFrame
        df = pd.DataFrame(csv_data)

        # Save to CSV
        if os.path.exists(csv_file_path):
            df.to_csv(csv_file_path, mode='a',index=False, header=False)
        else:
            df.to_csv(csv_file_path, mode='w',index=False, header=False)

    print(f"Encodings saved to {csv_file_path}")

        

def get_embedding_IDP(data_file, file_path, model_path,output_name,args, pretrain_config_path):
    """
    Loads data and generates embeddings for IDP sequences using a pre-trained BERT model.

    Args:
        data_file (str): Path to the data file containing the sequences.
        file_path (str): Directory where the generated embeddings will be saved.
        model_path (str): Path to the pre-trained BERT model checkpoint.

    Returns:
        None
    """
    # Load data and process sequences
    test_data = load_file_2_data(data_file)
    print("IDP-BERT processing sequences:", len(test_data))
    main(test_data, file_path, model_path,args,output_name, pretrain_config_path)
    print("Done")

DEFAULT_CONFIG = 'Paper/configs/idrlm_pretrain.yaml' #TODO

# Argument parser structure adapted with assistance from ChatGPT (OpenAI) and was reviewed and modified by the author.
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate IDR_LM embeddings for protein sequences in a FASTA file."
    )

    parser.add_argument(
        "--seq_input_path",
        type=Path,
        help="Path to the input FASTA file containing protein sequences."
    )

    parser.add_argument(
        "--model-path",
        type=Path,
        default=None,
        help=(
            "Optional path to a trained IDR_LM checkpoint. "
            "If omitted, a randomly initialized model will be used."
        )
    )
    
    parser.add_argument(
        "--pretrain-config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=(
            "Config file for pretrained model settings (Default should be /configs/idrlm_pretrain.yaml). "
        )
    )

    parser.add_argument(
        "--output-name",
        type=str,
        default=None,
        help=(
            "Base name for the output file. "
            "Defaults to 'randominit_<input_filename>' if not provided."
        )
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Directory where embedding files will be written."
    )

    return parser.parse_args()


if __name__ == "__main__":
    args_parse = parse_args()
    
    
    args_config = load_config(args_parse.pretrain_config)

    seq_input_path = args_parse.seq_input_path
    model_path = args_parse.model_path

    if args_parse.output_name is not None:
        output_name = args_parse.output_name
    else:
        if model_path is None:
            output_name = f"randominit_{seq_input_path.stem}"
        else:
            output_name = seq_input_path.stem

    #output_path = args_parse.out_dir / output_name

    get_embedding_IDP(
        data_file=str(seq_input_path),
        file_path=str(args_parse.out_dir),
        model_path=str(model_path) if model_path is not None else None,
        output_name=output_name,
        pretrain_config_path = str(args_parse.pretrain_config),
        args=args_config
    )