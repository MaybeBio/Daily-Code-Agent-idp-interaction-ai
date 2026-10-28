ulimit -c unlimited
[ -z "${output_dir}" ] && output_dir=output/condenseq/seq2img/

# Dataset
[ -z "${data_path}" ] && data_path=/hpc/reference/opencell/condenseq
[ -z "${split_key}" ] && split_key=cellfm_test
[ -z "${phase}" ] && phase='test'  # 'train', 'test'

[ -z "${img_resize}" ] && img_resize=160
[ -z "${img_crop_size}" ] && img_crop_size=160
[ -z "${seq_zero_mask_ratio}" ] && seq_zero_mask_ratio=0.0
[ -z "${img_type}" ] && img_type='GFP'
[ -z "${cell_image}" ] && cell_image='nucl'

# Transport parameters
[ -z "${path_type}" ] && path_type=Linear
[ -z "${prediction}" ] && prediction=velocity

# Model
## VAE
[ -z "${num_down_blocks}" ] && num_down_blocks=3
[ -z "${latent_channels}" ] && latent_channels=4
[ -z "${vae_block_out_channels}" ] && vae_block_out_channels='128,256,512'

## CELL-FM
[ -z "${img_mask_ratio}" ] && img_mask_ratio=0
[ -z "${cond_out_channels}" ] && cond_out_channels='32,64'
[ -z "${sample_size}" ] && sample_size=40
[ -z "${esm_embedding}" ] && esm_embedding='esmc_600m'
[ -z "${encoder_hidden_size}" ] && encoder_hidden_size=1152
[ -z "${max_protein_sequence_len}" ] && max_protein_sequence_len=2048
[ -z "${encoder_num_hidden_layers}" ] && encoder_num_hidden_layers=4
[ -z "${num_heads}" ] && num_heads=8
[ -z "${dim_head}" ] && dim_head=64
[ -z "${dropout}" ] && dropout=0
[ -z "${final_dropout}" ] && final_dropout=0
[ -z "${encoder_patch_size}" ] && encoder_patch_size=4

### Image generator
[ -z "${img_generator_num_layers}" ] && img_generator_num_layers=4
[ -z "${img_generator_patch_size}" ] && img_generator_patch_size=2
[ -z "${attention_head_dim}" ] && attention_head_dim=64
[ -z "${num_attention_heads}" ] && num_attention_heads=8

### Image decoder
[ -z "${img_decoder_num_hidden_layers}" ] && img_decoder_num_hidden_layers=2
[ -z "${img_decoder_hidden_size}" ] && img_decoder_hidden_size=512
[ -z "${img_decoder_dim_head}" ] && img_decoder_dim_head=64
[ -z "${img_decoder_num_heads}" ] && img_decoder_num_heads=8

[ -z "${cell_image_ratio}" ] && cell_image_ratio=0.5

# Training
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_condenseq/vae/checkpoint-50000/pytorch_model.bin
[ -z "${loadcheck_path}" ] && loadcheck_path=pretrain_condenseq/cellfm_seq2img/checkpoint-50000/pytorch_model.bin

# Evaluation
[ -z "${num_steps}" ] && num_steps=100


python cell_fm/tasks/cell_fm_cs/generate_img.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --phase $phase \
            --img_resize $img_resize \
            --img_crop_size $img_crop_size \
            --img_type $img_type \
            --seq_zero_mask_ratio $seq_zero_mask_ratio \
            --cell_image $cell_image \
            --path_type $path_type \
            --prediction $prediction \
            --num_down_blocks $num_down_blocks \
            --latent_channels $latent_channels \
            --vae_block_out_channels $vae_block_out_channels \
            --img_mask_ratio $img_mask_ratio \
            --cond_out_channels $cond_out_channels \
            --sample_size $sample_size \
            --esm_embedding $esm_embedding \
            --encoder_hidden_size $encoder_hidden_size \
            --max_protein_sequence_len $max_protein_sequence_len \
            --encoder_num_hidden_layers $encoder_num_hidden_layers \
            --num_heads $num_heads \
            --dim_head $dim_head \
            --dropout $dropout \
            --final_dropout $final_dropout \
            --encoder_patch_size $encoder_patch_size \
            --img_generator_num_layers $img_generator_num_layers \
            --img_generator_patch_size $img_generator_patch_size \
            --attention_head_dim $attention_head_dim \
            --num_attention_heads $num_attention_heads \
            --img_decoder_num_hidden_layers $img_decoder_num_hidden_layers \
            --img_decoder_hidden_size $img_decoder_hidden_size \
            --img_decoder_num_heads $img_decoder_num_heads \
            --img_decoder_dim_head $img_decoder_dim_head \
            --cell_image_ratio $cell_image_ratio \
            --vae_loadcheck_path $vae_loadcheck_path \
            --loadcheck_path $loadcheck_path \
            --num_steps $num_steps \
            --seed 6 \
            --infer \

            # --ft \
            # --ifresume \
            # --bf16 \
            # --fp16 \