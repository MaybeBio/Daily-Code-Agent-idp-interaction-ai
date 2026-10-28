ulimit -c unlimited
[ -z "${n_gpu}" ] && n_gpu=$(nvidia-smi -L | wc -l)

# Wandb
export WANDB_RUN_NAME=PT_VAE3D_OC_192_nd3_lc16_KL1e-5_R1
export WANDB_PROJECT=CELL-FM-3D
[ -z "${output_dir}" ] && output_dir=pretrain_opencell_3d/$WANDB_RUN_NAME

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/opencell'
[ -z "${split_key}" ] && split_key='all'

# Model
[ -z "${num_down_blocks}" ] && num_down_blocks=3
[ -z "${latent_channels}" ] && latent_channels=16
[ -z "${vae_block_out_channels}" ] && vae_block_out_channels='64,128,256'
[ -z "${input_spatial_size}" ] && input_spatial_size='48,192,192'
[ -z "${norm_num_groups}" ] && norm_num_groups=32
[ -z "${layers_per_block}" ] && layers_per_block=2

# Loss
[ -z "${recon_loss_coeff}" ] && recon_loss_coeff=1.0
[ -z "${kl_loss_coeff}" ] && kl_loss_coeff=1e-5

# Training
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_opencell_3d/PT_VAE3D_OC_192_nd3_lc16_KL1e-4_R1/checkpoint-100000/pytorch_model.bin
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=4
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=2

[ -z "${num_train_epochs}" ] && num_train_epochs=5000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=100000
[ -z "${save_steps}" ] && save_steps=10000

[ -z "${MASTER_PORT}" ] && MASTER_PORT=21829
[ -z "${MASTER_ADDR}" ] && MASTER_ADDR=127.0.0.1

DISTRIBUTED_ARGS="--nproc_per_node $n_gpu \
                  --master_port $MASTER_PORT \
                  --master_addr $MASTER_ADDR"

python -m torch.distributed.run $DISTRIBUTED_ARGS cell_fm/tasks/vae_3d/pretrain_vae_3d.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --data_aug \
            --num_down_blocks $num_down_blocks \
            --latent_channels $latent_channels \
            --vae_block_out_channels $vae_block_out_channels \
            --input_spatial_size $input_spatial_size \
            --norm_num_groups $norm_num_groups \
            --layers_per_block $layers_per_block \
            --recon_loss_coeff $recon_loss_coeff \
            --kl_loss_coeff $kl_loss_coeff \
            --vae_loadcheck_path $vae_loadcheck_path \
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
            --seed 6 \
            --wandb \
            --ft \

            # --ft \
            # --bf16 \
            # --ifresume \
