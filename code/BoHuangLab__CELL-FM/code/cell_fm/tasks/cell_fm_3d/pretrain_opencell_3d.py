# -*- coding: utf-8 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from cell_fm.data.opencell_3d_crop_data.dataset import OpenCell3DCropDataset
from cell_fm.models.cell_fm_3d.cell_fm_3d_config import CELLFM3DConfig
from cell_fm.models.cell_fm_3d.cell_fm_3d_model import CELLFM3DModel
from cell_fm.utils.cli_utils import cli
from transformers import TrainingArguments
from cell_fm.pipeline.cd_accelerator.accelerator import CDTrainer
from cell_fm.pipeline.ema import EMACallback, add_ema_callback
from cell_fm.logging.loggers import CELLFMLoggingCallback
from cell_fm.logging import logger


@cli(CELLFM3DConfig)
def main(args) -> None:
    trainset = OpenCell3DCropDataset(args, split_key=args.split_key)
    model = CELLFM3DModel(config=CELLFM3DConfig(**vars(args)))

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

    if args.use_ema:
        add_ema_callback(
            trainer,
            EMACallback(
                decay=args.ema_decay,
                use_ema_warmup=args.ema_warmup,
                inv_gamma=args.ema_inv_gamma,
                power=args.ema_power,
                update_after_step=args.ema_update_after_step,
            ),
        )

    trainer.train(resume_from_checkpoint=args.ifresume)


if __name__ == "__main__":
    main()
