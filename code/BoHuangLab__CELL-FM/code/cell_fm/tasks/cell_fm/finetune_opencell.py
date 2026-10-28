# -*- coding: utf-8 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.opencell_data.dataset import OpenCellDataset
from cell_fm.models.cell_fm.cell_fm_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_model import CELLFMModel
from cell_fm.utils.cli_utils import cli
from transformers import TrainingArguments
from cell_fm.pipeline.cd_accelerator.accelerator import CDTrainer
from cell_fm.logging.loggers import CELLFMLoggingCallback
from cell_fm.logging import logger


@cli(CELLFMConfig)
def main(args) -> None:
    trainset = OpenCellDataset(args, split_key=args.split_key)
    model = CELLFMModel(config=CELLFMConfig(**vars(args)), loss_fn=UniDiffCriterions)

    logger.info(args)

    training_args = TrainingArguments(
        output_dir=args.output_dir, 
        learning_rate=args.learning_rate, 
        weight_decay=args.weight_decay, 
        gradient_accumulation_steps=args.gradient_accumulation_steps, 
        per_device_train_batch_size=args.per_device_train_batch_size, 
        per_device_eval_batch_size=args.per_device_eval_batch_size, 
        num_train_epochs=args.num_train_epochs, 
        fp16=args.fp16, 
        bf16=args.bf16, 
        logging_dir=args.logging_dir, 
        logging_steps=args.logging_steps, 
        max_steps=args.max_steps, 
        warmup_steps=args.warmup_steps, 
        save_steps=args.save_steps, 
        seed=args.seed, 
        dataloader_num_workers=args.dataloader_num_workers, 
        report_to='wandb' if args.wandb else 'none', 
        disable_tqdm=True, 
        remove_unused_columns=False, 
        overwrite_output_dir=True, 
        log_level='debug', 
        include_inputs_for_metrics=False, 
        save_safetensors=False, 
    )

    trainer = CDTrainer(
        model=model, 
        args=training_args, 
        train_dataset=trainset, 
        data_collator=trainset.collate, 
        callbacks=[CELLFMLoggingCallback()], 
    )
    trainer.train(resume_from_checkpoint=args.ifresume)

if __name__ == "__main__":
    main()

