ulimit -c unlimited

# Output
[ -z "${output_dir}" ] && output_dir=/hpc/reference/opencell/opencell/3d_crop_virtual_staining

# Dataset
[ -z "${data_path}" ] && data_path='/hpc/reference/opencell/opencell'
[ -z "${split_key}" ] && split_key='all'
[ -z "${phase}" ] && phase='test'

# Transport parameters
[ -z "${path_type}" ] && path_type=Linear
[ -z "${prediction}" ] && prediction=velocity

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

# Cell image conditioning
[ -z "${cell_image}" ] && cell_image='nucl'
[ -z "${cond_out_channels}" ] && cond_out_channels='32'

# UNet around SD3: width per level (latent 24x96x96 -> 12x48x48), ResBlocks per level.
# Must match the run that wrote the checkpoint.
[ -z "${unet_block_out_channels}" ] && unet_block_out_channels='64,128'
[ -z "${unet_layers_per_block}" ] && unet_layers_per_block=2

# Image generator (SD3-3D)
[ -z "${img_generator_num_layers}" ] && img_generator_num_layers=24
[ -z "${attention_head_dim}" ] && attention_head_dim=64
[ -z "${num_attention_heads}" ] && num_attention_heads=18

# Checkpoint paths
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_opencell_3d/PT_VAE3D_OC_192_KL1e-4/checkpoint-50000/pytorch_model.bin
# [ -z "${loadcheck_path}" ] && loadcheck_path=pretrain_opencell_3d/PT_OC3D_CELLFM3D_R14/checkpoint-100000/pytorch_model.bin
[ -z "${loadcheck_path}" ] && loadcheck_path=pretrain_opencell_3d/PT_OC3D_CELLFM3D_R19/checkpoint-50000/ema_pytorch_model.bin

# Evaluation
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=1
[ -z "${num_steps}" ] && num_steps=100


python cell_fm/tasks/cell_fm_3d/virtual_staining_opencell_3d.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --phase $phase \
            --path_type $path_type \
            --prediction $prediction \
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
            --esm_embedding $esm_embedding \
            --encoder_hidden_size $encoder_hidden_size \
            --max_protein_sequence_len $max_protein_sequence_len \
            --img_generator_num_layers $img_generator_num_layers \
            --attention_head_dim $attention_head_dim \
            --num_attention_heads $num_attention_heads \
            --vae_loadcheck_path $vae_loadcheck_path \
            --loadcheck_path $loadcheck_path \
            --per_device_eval_batch_size $per_device_eval_batch_size \
            --seed 6 \
            --num_steps $num_steps \
            --infer \
