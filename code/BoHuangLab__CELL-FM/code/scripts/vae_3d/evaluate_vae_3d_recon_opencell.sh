ulimit -c unlimited

[ -z "${output_dir}" ] && output_dir='./output/opencell_3d/PT_VAE3D_OC_192_KL1e-4'

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/opencell'
[ -z "${split_key}" ] && split_key='all'

# Model
[ -z "${num_down_blocks}" ] && num_down_blocks=2
[ -z "${latent_channels}" ] && latent_channels=4
[ -z "${vae_block_out_channels}" ] && vae_block_out_channels='64,128'
[ -z "${input_spatial_size}" ] && input_spatial_size='48,192,192'
[ -z "${norm_num_groups}" ] && norm_num_groups=32
[ -z "${layers_per_block}" ] && layers_per_block=2

# Loss
[ -z "${recon_loss_coeff}" ] && recon_loss_coeff=1.0
[ -z "${kl_loss_coeff}" ] && kl_loss_coeff=1e-4

# Training (required by VAE3DConfig / cli but not used during inference)
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_opencell_3d/PT_VAE3D_OC_192_KL1e-4/checkpoint-50000/pytorch_model.bin
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=2
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=2

[ -z "${num_train_epochs}" ] && num_train_epochs=5000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=100000
[ -z "${save_steps}" ] && save_steps=10000

python cell_fm/tasks/vae_3d/evaluate_vae_3d_recon_opencell.py \
    --output_dir $output_dir \
    --data_path $data_path \
    --split_key $split_key \
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
    --infer \

    # --bf16 \
