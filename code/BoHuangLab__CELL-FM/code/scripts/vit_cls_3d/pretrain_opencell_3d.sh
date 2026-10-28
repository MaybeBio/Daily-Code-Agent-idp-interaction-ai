ulimit -c unlimited
[ -z "${n_gpu}" ] && n_gpu=$(nvidia-smi -L | wc -l)

# Wandb
export WANDB_RUN_NAME=PT_OC3D_ViT_cls
export WANDB_PROJECT=CELL-FM
[ -z "${output_dir}" ] && output_dir=pretrain_opencell_3d/$WANDB_RUN_NAME

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/opencell'
[ -z "${split_key}" ] && split_key='all'

# Input
[ -z "${input_spatial_size}" ] && input_spatial_size='48,192,192'
[ -z "${in_channels}" ] && in_channels=2

# Model
[ -z "${patch_d}" ] && patch_d=8
[ -z "${patch_size}" ] && patch_size=8
[ -z "${num_layers}" ] && num_layers=12
[ -z "${num_heads}" ] && num_heads=8
[ -z "${hidden_dim}" ] && hidden_dim=512
[ -z "${mlp_dim}" ] && mlp_dim=2048

# Training
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=64
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=4

[ -z "${num_train_epochs}" ] && num_train_epochs=500000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=100000
[ -z "${save_steps}" ] && save_steps=10000

[ -z "${MASTER_PORT}" ] && MASTER_PORT=12324
[ -z "${MASTER_ADDR}" ] && MASTER_ADDR=127.0.0.1

DISTRIBUTED_ARGS="--nproc_per_node $n_gpu \
                  --master_port $MASTER_PORT \
                  --master_addr $MASTER_ADDR"

python -m torch.distributed.run $DISTRIBUTED_ARGS cell_fm/tasks/vit_cls_3d/pretrain_opencell_3d.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --data_aug \
            --input_spatial_size $input_spatial_size \
            --in_channels $in_channels \
            --patch_d $patch_d \
            --patch_size $patch_size \
            --num_layers $num_layers \
            --num_heads $num_heads \
            --hidden_dim $hidden_dim \
            --mlp_dim $mlp_dim \
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

            # --loadcheck_path $loadcheck_path \
            # --ft \
            # --ifresume \
            # --bf16 \
