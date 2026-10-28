ulimit -c unlimited

[ -z "${output_dir}" ] && output_dir=./output/opencell_3d/embedding/PT_OC3D_ViT_cls_R2_90k

# Dataset
[ -z "${data_path}" ] && data_path=/hpc/reference/opencell/opencell
[ -z "${split_key}" ] && split_key=all

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

# Inference
[ -z "${loadcheck_path}" ] && loadcheck_path=pretrain_opencell_3d/PT_OC3D_ViT_cls_R2/checkpoint-90000/pytorch_model.bin
[ -z "${per_device_eval_batch_size}" ] && per_device_eval_batch_size=8

python cell_fm/tasks/vit_cls_3d/embed_opencell_3d.py \
            --output_dir $output_dir \
            --data_path $data_path \
            --split_key $split_key \
            --input_spatial_size $input_spatial_size \
            --in_channels $in_channels \
            --patch_d $patch_d \
            --patch_size $patch_size \
            --num_layers $num_layers \
            --num_heads $num_heads \
            --hidden_dim $hidden_dim \
            --mlp_dim $mlp_dim \
            --loadcheck_path $loadcheck_path \
            --per_device_eval_batch_size $per_device_eval_batch_size \
            --seed 6 \
            --infer \
