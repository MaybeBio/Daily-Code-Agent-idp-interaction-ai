# -*- coding: utf-8 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from cell_fm.data.opencell_3d_crop_data.dataset import OpenCell3DCropImageOnlyDataset
from cell_fm.models.vit_cls_3d.config import ViT3DConfig
from cell_fm.models.vit_cls_3d.model import ViT3DModel
from cell_fm.utils.cli_utils import cli
from transformers import TrainingArguments, Trainer
from cell_fm.logging.loggers import VanillaLoggingCallback
from cell_fm.logging import logger


@cli(ViT3DConfig)
def main(args) -> None:
    trainset = OpenCell3DCropImageOnlyDataset(args, split_key=args.split_key)

    config = ViT3DConfig(**vars(args))
    config.num_classes = len(trainset.gene_names)
    model = ViT3DModel(config=config)

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

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=trainset,
        data_collator=trainset.collate,
        callbacks=[VanillaLoggingCallback()],
    )
    trainer.train(resume_from_checkpoint=args.ifresume)


if __name__ == "__main__":
    main()
