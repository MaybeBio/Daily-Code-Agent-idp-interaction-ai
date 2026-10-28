ulimit -c unlimited
[ -z "${n_gpu}" ] && n_gpu=$(nvidia-smi -L | wc -l)

# Wandb
export WANDB_RUN_NAME=PT_CondenSeq_CELLFM_S1
export WANDB_PROJECT=CELL-FM
[ -z "${output_dir}" ] && output_dir=pretrain_condenseq/$WANDB_RUN_NAME

# Dataset
[ -z "${data_path}" ] && data_path=/hpc/reference/opencell/condenseq
[ -z "${split_key}" ] && split_key=cellfm_train
[ -z "${phase}" ] && phase='train'  # 'train', 'test'

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
[ -z "${img_mask_ratio}" ] && img_mask_ratio=0.5
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

# Loss
[ -z "${seq_loss_coeff}" ] && seq_loss_coeff=1.0
[ -z "${img_diff_loss_coeff}" ] && img_diff_loss_coeff=1.0
[ -z "${img_recon_loss_coeff}" ] && img_recon_loss_coeff=1.0

# Training
[ -z "${vae_loadcheck_path}" ] && vae_loadcheck_path=pretrain_condenseq/vae/checkpoint-50000/pytorch_model.bin
[ -z "${loadcheck_path}" ] && loadcheck_path=.
[ -z "${learning_rate}" ] && learning_rate=3e-4
[ -z "${weight_decay}" ] && weight_decay=0.0
[ -z "${gradient_accumulation_steps}" ] && gradient_accumulation_steps=1
[ -z "${per_device_train_batch_size}" ] && per_device_train_batch_size=512
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=128

[ -z "${num_train_epochs}" ] && num_train_epochs=5000
[ -z "${logging_dir}" ] && logging_dir=$output_dir
[ -z "${logging_steps}" ] && logging_steps=100
[ -z "${warmup_steps}" ] && warmup_steps=1000
[ -z "${max_steps}" ] && max_steps=50000
[ -z "${save_steps}" ] && save_steps=10000

[ -z "${dataloader_num_workers}" ] && dataloader_num_workers=32

[ -z "${MASTER_PORT}" ] && MASTER_PORT=29182
[ -z "${MASTER_ADDR}" ] && MASTER_ADDR=127.0.0.1

DISTRIBUTED_ARGS="--nproc_per_node $n_gpu \
                  --master_port $MASTER_PORT \
                  --master_addr $MASTER_ADDR"


python -m torch.distributed.run $DISTRIBUTED_ARGS cell_fm/tasks/cell_fm_cs/pretrain_condenseq.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --phase $phase \
            --data_aug \
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
            --seq_loss_coeff $seq_loss_coeff \
            --img_diff_loss_coeff $img_diff_loss_coeff \
            --img_recon_loss_coeff $img_recon_loss_coeff \
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
            --dataloader_num_workers $dataloader_num_workers \
            --seed 666666 \

            # --ft \
            # --ifresume \
            # --bf16 \
            # --fp16 \