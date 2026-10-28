ulimit -c unlimited
[ -z "${n_gpu}" ] && n_gpu=$(nvidia-smi -L | wc -l)

# Wandb
export WANDB_RUN_NAME=PT_VAE_HPA
export WANDB_PROJECT=CELL-FM
[ -z "${output_dir}" ] && output_dir=pretrain_hpa/$WANDB_RUN_NAME

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/human_protein_atlas'
[ -z "${split_key}" ] && split_key=cellfm_train

[ -z "${img_crop_size}" ] && img_crop_size=1024
[ -z "${img_resize}" ] && img_resize=256

# Model
[ -z "${num_down_blocks}" ] && num_down_blocks=3
[ -z "${latent_channels}" ] && latent_channels=4
[ -z "${vae_block_out_channels}" ] && vae_block_out_channels='128,256,512'

# Loss
[ -z "${recon_loss_coeff}" ] && recon_loss_coeff=1.0
[ -z "${kl_loss_coeff}" ] && kl_loss_coeff=1e-4

# Training
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=.
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=64
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=128

[ -z "${num_train_epochs}" ] && num_train_epochs=5000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=50000
[ -z "${save_steps}" ] && save_steps=1000

[ -z "${MASTER_PORT}" ] && MASTER_PORT=21828
[ -z "${MASTER_ADDR}" ] && MASTER_ADDR=127.0.0.1

DISTRIBUTED_ARGS="--nproc_per_node $n_gpu \
                  --master_port $MASTER_PORT \
                  --master_addr $MASTER_ADDR"

python -m torch.distributed.run $DISTRIBUTED_ARGS cell_fm/tasks/vae/pretrain_vae.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --data_aug \
            --img_crop_size $img_crop_size \
            --img_resize $img_resize \
            --num_down_blocks $num_down_blocks \
            --latent_channels $latent_channels \
            --vae_block_out_channels $vae_block_out_channels \
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
            --seed 666666 \
            --ft \

            # --fp16 \
            # --bf16 \
            # --ifresume \