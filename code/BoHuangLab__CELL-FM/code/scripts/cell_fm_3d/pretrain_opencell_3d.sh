ulimit -c unlimited
[ -z "${n_gpu}" ] && n_gpu=$(nvidia-smi -L | wc -l)

# Wandb
export WANDB_RUN_NAME=PT_OC3D_CELLFM3D_R19
export WANDB_PROJECT=CELL-FM-3D
[ -z "${output_dir}" ] && output_dir=pretrain_opencell_3d/$WANDB_RUN_NAME

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/opencell'
[ -z "${split_key}" ] && split_key='all'
[ -z "${phase}" ] && phase='train'
[ -z "${seq_zero_mask_ratio}" ] && seq_zero_mask_ratio=0.0

# Transport parameters
[ -z "${path_type}" ] && path_type=Linear
[ -z "${prediction}" ] && prediction=velocity
[ -z "${timestep_sampler}" ] && timestep_sampler=logit_normal
[ -z "${logit_mean}" ] && logit_mean=0.0
[ -z "${logit_std}" ] && logit_std=1.0

# VAE (3D, pretrained)
[ -z "${num_down_blocks}" ] && num_down_blocks=2
[ -z "${latent_channels}" ] && latent_channels=4
[ -z "${vae_block_out_channels}" ] && vae_block_out_channels='64,128'
[ -z "${input_spatial_size}" ] && input_spatial_size='48,192,192'
[ -z "${norm_num_groups}" ] && norm_num_groups=32
[ -z "${layers_per_block}" ] && layers_per_block=2

# SD3 patch at the UNet's lowest resolution (12x48x48 -> 3x24x24 = 1728 tokens)
[ -z "${patch_d}" ] && patch_d=4
[ -z "${patch_size}" ] && patch_size=2

# Sequence / ESM
[ -z "${esm_embedding}" ] && esm_embedding='esmc_600m'
[ -z "${encoder_hidden_size}" ] && encoder_hidden_size=1152
[ -z "${max_protein_sequence_len}" ] && max_protein_sequence_len=2048
[ -z "${seq_length_control}" ] && seq_length_control=crop

# Cell image conditioning
[ -z "${cell_image}" ] && cell_image='nucl'
[ -z "${cond_out_channels}" ] && cond_out_channels='32'

# UNet around SD3: width per level (latent 24x96x96 -> 12x48x48), ResBlocks per level
[ -z "${unet_block_out_channels}" ] && unet_block_out_channels='64,128'
[ -z "${unet_layers_per_block}" ] && unet_layers_per_block=2

# EMA
[ -z "${ema_decay}" ] && ema_decay=0.9999
[ -z "${ema_inv_gamma}" ] && ema_inv_gamma=1.0
[ -z "${ema_power}" ] && ema_power=0.6666666666666666
[ -z "${ema_update_after_step}" ] && ema_update_after_step=0

# Image generator (SD3-3D)
[ -z "${img_generator_num_layers}" ] && img_generator_num_layers=24
[ -z "${attention_head_dim}" ] && attention_head_dim=64
[ -z "${num_attention_heads}" ] && num_attention_heads=18

# Checkpoint paths
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_opencell_3d/PT_VAE3D_OC_192_KL1e-4/checkpoint-50000/pytorch_model.bin
[ -z "${loadcheck_path}" ] && loadcheck_path=pretrain_opencell_3d/PT_OC3D_CELLFM3D_R18/checkpoint-60000/ema_pytorch_model.bin

# Training
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=16
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=1

[ -z "${num_train_epochs}" ] && num_train_epochs=5000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=100000
[ -z "${save_steps}" ] && save_steps=10000

[ -z "${MASTER_PORT}" ] && MASTER_PORT=12325
[ -z "${MASTER_ADDR}" ] && MASTER_ADDR=127.0.0.1

DISTRIBUTED_ARGS="--nproc_per_node $n_gpu \
                  --master_port $MASTER_PORT \
                  --master_addr $MASTER_ADDR"

python -m torch.distributed.run $DISTRIBUTED_ARGS cell_fm/tasks/cell_fm_3d/pretrain_opencell_3d.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --phase $phase \
            --data_aug \
            --seq_zero_mask_ratio $seq_zero_mask_ratio \
            --path_type $path_type \
            --prediction $prediction \
            --timestep_sampler $timestep_sampler \
            --logit_mean $logit_mean \
            --logit_std $logit_std \
            --num_down_blocks $num_down_blocks \
            --latent_channels $latent_channels \
            --vae_block_out_channels $vae_block_out_channels \
            --input_spatial_size $input_spatial_size \
            --norm_num_groups $norm_num_groups \
            --layers_per_block $layers_per_block \
            --patch_d $patch_d \
            --patch_size $patch_size \
            --cell_image $cell_image \
            --cond_out_channels $cond_out_channels \
            --unet_block_out_channels $unet_block_out_channels \
            --unet_layers_per_block $unet_layers_per_block \
            --use_ema \
            --ema_decay $ema_decay \
            --ema_inv_gamma $ema_inv_gamma \
            --ema_power $ema_power \
            --ema_update_after_step $ema_update_after_step \
            --esm_embedding $esm_embedding \
            --encoder_hidden_size $encoder_hidden_size \
            --max_protein_sequence_len $max_protein_sequence_len \
            --seq_length_control $seq_length_control \
            --img_generator_num_layers $img_generator_num_layers \
            --attention_head_dim $attention_head_dim \
            --num_attention_heads $num_attention_heads \
            --vae_loadcheck_path $vae_loadcheck_path \
            --loadcheck_path $loadcheck_path \
            --learning_rate $learning_rate \
            --weight_decay $weight_decay \
            --gradient_accumulation_steps $gradient_accumulation_steps \
            --per_device_train_batch_size $per_device_train_batch_size \
            --per_device_eval_batch_size $per_device_eval_batch_size \
            --num_train_epochs $num_train_epochs \
            --logging_dir $logging_dir \
            --logging_steps $logging_steps \
            --warmup_steps $warmup_steps \
            --max_steps $max_steps \
            --save_steps $save_steps \
            --seed 42 \
            --wandb \
            --ft \
            --bf16 \

            # --ft \
            # --ifresume \
            # --bf16 \
