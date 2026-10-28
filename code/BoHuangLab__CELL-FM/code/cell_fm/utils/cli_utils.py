# -*- coding: utf-8 -*-
import inspect
import os
from argparse import ArgumentParser
from functools import wraps

from cell_fm.logging import logger
from cell_fm.pipeline.accelerator.trainer import seed_everything
from cell_fm.utils import arg_utils, dist_utils, env_init

import wandb  # isort:skip


def cli(*cfg_classes_and_funcs):
    def decorator(main):
        @wraps(main)
        def wrapper():
            parser = ArgumentParser()
            parser.add_argument(
                "--wandb",
                action="store_true",
                default=False,
                help="Enable Weights & Biases logging. Disabled by default.",
            )
            cfg_classes = []
            cfg_funcs = []
            for cfg in cfg_classes_and_funcs:
                if inspect.isclass(cfg):
                    cfg_classes.append(cfg)
                else:
                    cfg_funcs.append(cfg)
            for cfg_func in cfg_funcs:
                parser = cfg_func(parser)
            parser = arg_utils.add_dataclass_to_parser(cfg_classes, parser)
            args = parser.parse_args()

            logger.info(args)

            seed_everything(args.seed)

            env_init.set_env(args)

            if dist_utils.is_master_node():
                if args.wandb:
                    wandb_api_key = os.getenv("WANDB_API_KEY")
                    if wandb_api_key:
                        wandb.login(key=wandb_api_key)

                    wandb_project = os.getenv("WANDB_PROJECT")
                    wandb_run_name = os.getenv("WANDB_RUN_NAME")
                    wandb_team = os.getenv("WANDB_TEAM")
                    wandb_group = os.getenv("WANDB_GROUP")

                    args.wandb_team = getattr(args, "wandb_team", wandb_team)
                    args.wandb_group = getattr(args, "wandb_group", wandb_group)
                    args.wandb_project = getattr(args, "wandb_project", wandb_project)

                    wandb.init(
                        project=wandb_project,
                        group=wandb_group,
                        name=wandb_run_name,
                        config=args,
                    )
                else:
                    logger.info("Wandb disabled; pass --wandb to enable it")

            logger.success(
                "====================================Start!===================================="
            )
            try:
                main(args)
            except Exception as e:
                logger.exception(e)
                logger.error(
                    "====================================Fail!===================================="
                )
                exit()

            logger.success(
                "====================================Done!===================================="
            )

        return wrapper

    return decorator
