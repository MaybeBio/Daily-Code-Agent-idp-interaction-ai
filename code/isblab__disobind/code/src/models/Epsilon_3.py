###########################################################
############Attention coupled with Convolutions############
###########################################################
import numpy as np

import torch
from torch import nn

from src.models.get_activation import get_activation
from src.models.get_layers import ( create_projection_layers, 
                                    create_upsampling_layers, 
                                    create_downsampling_layers )


class Epsilon_3( nn.Module ):
    def __init__( self, emb_size, projection_layer, output_dim, 
                    activation1, activation2, input_layer,
                    num_samples, num_hid_layers, bias,
                    max_len, dropouts, norm,
                    temperature, output_layer, objective, device ):
        super( Epsilon_3, self ).__init__()
        self.emb_size = emb_size
        self.projection_dim = projection_layer[0]
        self.input_layer_params = input_layer
        self.num_samples = num_samples
        self.num_upsample_layers = num_hid_layers[0]
        self.num_downsample_layers = num_hid_layers[1]
        self.num_blocks = num_hid_layers[2]
        self.scale_factor = num_hid_layers[3]
        self.hidden_block = num_hid_layers[4:]
        self.bias = bias
        self.norm = norm
        self.max_len = max_len
        self.dropout1 = nn.Dropout( p = dropouts[0] )
        self.dropout2 = nn.Dropout( p = dropouts[1] )
        self.us_dropout_prob = dropouts[2]
        self.ds_dropout_prob = dropouts[3]
        self.mc_dropout = nn.Dropout( p = dropouts[4] )
        self.activ1_name, self.activ1_param = activation1
        self.activation1 = get_activation( self.activ1_name, device, self.activ1_param )
        self.activ2_name, self.activ2_param = activation2
        self.activation2 = get_activation( self.activ2_name, device )
        self.apply_activation2 = activation2[1]
        self.output_dim = output_dim
        self.out_layer_name = output_layer
        self.device = device
        self.objective = objective

        # No. of features after concatenation.
        self.state0 = 2*self.projection_dim
        # No. of features after downsampling.
        self.state1 = self.state0//self.scale_factor**( self.num_downsample_layers )
        # # No. of features after upsampling.
        self.state2 = self.state1*self.scale_factor**( self.num_upsample_layers )

        if self.hidden_block[0] == "residual":
            if self.norm[1] == "LN":
                self.norm_layer = nn.LayerNorm( self.state2 )
            if self.norm[1] == "IN":
                self.norm_layer = nn.InstanceNorm1d( self.state2 )
            if self.norm[1] == "BN":
                self.norm_layer = nn.BatchNorm1d( self.state2 )

        self.projection_layer1 = create_projection_layers( proj_type = projection_layer[1],
                                                            emb_size = self.emb_size, 
                                                            projection_dim = self.projection_dim, 
                                                            bias = projection_layer[2], 
                                                            multiplier = projection_layer[3],
                                                            activation = self.activation1, 
                                                            device = self.device )
        if projection_layer[4] == "separate":
            self.projection_layer2 = create_projection_layers( proj_type = projection_layer[1],
                                                                emb_size = self.emb_size, 
                                                                projection_dim = self.projection_dim, 
                                                                bias = projection_layer[2], 
                                                                multiplier = projection_layer[3],
                                                                activation = self.activation1, 
                                                                device = self.device )
        else:
            self.projection_layer2 = None


        if self.input_layer_params[2] == "lin":
            if "bin" in self.objective[0]:
                in_feat = self.max_len//self.objective[1]
            else:
                in_feat = self.max_len
            self.interface = nn.Sequential( 
                                            nn.Linear( in_features = in_feat, out_features = 1, bias = True )
                                            )

        self.downsampling_layers =  create_downsampling_layers( num_downsample_layers = self.num_downsample_layers, 
                                                                state = self.state0, 
                                                                scale = self.scale_factor,
                                                                bias = self.bias, 
                                                                activation = [self.activ1_name, self.activation1], 
                                                                ds_dropout_prob = self.ds_dropout_prob, 
                                                                norm = self.norm,
                                                                device = self.device )

        self.upsampling_layers =  create_upsampling_layers( num_upsample_layers = self.num_upsample_layers, 
                                                            state = self.state1, 
                                                            scale = self.scale_factor,
                                                            bias = self.bias, 
                                                            activation = [self.activ1_name, self.activation1], 
                                                            us_dropout_prob = self.us_dropout_prob, 
                                                            norm = self.norm,
                                                            device = self.device )

        # if self.objective[2] == "avg":
        #     self.bin_residues = nn.AvgPool1d( kernel_size = self.objective[1], stride = self.objective[1] )
        # elif self.objective[2] == "max":
        #     self.bin_residues = nn.MaxPool1d( kernel_size = self.objective[1], stride = self.objective[1] )


        if temperature != None:
            print( f"Using temperature scaling with Temperature = {temperature}..." )
            self.temperature = nn.Parameter( torch.tensor( [temperature] ) )
            # self.temperature = torch.tensor( [temperature] ).to( self.device )
            # self.temperature = temperature
        else:
            self.temperature = None        
        
        print( f"Output layer = {self.out_layer_name}" )
        self.contact = nn.Linear( in_features = self.state2, 
                                    out_features = self.output_dim, 
                                    bias = self.bias )



    def forward( self, e1, e2, interaction_mask ):
        # Input to projection block.
        z1, z2 = self.projection_block( e1, e2 )

        # Get interaction tensor from projected embeddings z1, z2.
        interaction_tensor = self.get_interaction_tensor( z1, z2 )

        if "interface" in self.objective[0]:
            o = self.interface_block( interaction_tensor, interaction_mask )
        elif "interaction" in self.objective[0]:
            n, l1, l2, c = interaction_tensor.shape
            o = interaction_tensor.view( n, -1, c )

        if self.hidden_block[0] == "vanilla":
            o = self.vanilla_block( o )
        elif self.hidden_block[0] == "residual":
            o = self.residual_block( o )
        
        o = self.output_block( o )

        return o



    def projection_block( self, e1, e2 ):
        if self.projection_layer2 != None:
            z1 = self.projection_layer1( e1 )
            z2 = self.projection_layer2( e2 )
        
        else:
            z1 = self.projection_layer1( e1 )
            z2 = self.projection_layer1( e2 )

        # Coarse grain after projection.
        if self.objective[3]:
            z1 = torch.permute( z1, ( 0, 2, 1 ) )
            z2 = torch.permute( z2, ( 0, 2, 1 ) )
            
            z1 = self.bin_residues( z1 )
            z2 = self.bin_residues( z2 )
            
            z1 = torch.permute( z1, ( 0, 2, 1 ) )
            z2 = torch.permute( z2, ( 0, 2, 1 ) )

        z1 = self.dropout1( z1 )
        z2 = self.dropout1( z2 )

        return z1, z2


    def avg2d( self, interaction_tensor, interaction_mask ):
        """
        Given input tensor [N,H,W,C], compute avg2d along H and W to obtain
            interface tensors:
            I1 -> [N,H,C] and I2 -> [N,W,C]
        Ignore the padding elements from mean computation.
        """
        eps = 1e-4
        masked_tensor = interaction_tensor*interaction_mask

        total_I1 = torch.sum( interaction_mask, dim = 2 )
        total_I2 = torch.sum( interaction_mask, dim = 1 )

        sum_I1 = torch.sum( masked_tensor, dim = 2 )
        sum_I2 = torch.sum( masked_tensor, dim = 1 )

        I1 = sum_I1/( total_I1 + eps )
        I2 = sum_I2/( total_I2 + eps )

        I = torch.cat( [I1, I2], axis = 1 )
        return I



    def interface_block( self, interaction_tensor, interaction_mask ):
        if self.input_layer_params[2] == "avg1d":
            I = torch.mean( interaction_tensor, axis = 2 )

        elif self.input_layer_params[2] == "avg2d":
            # I1 = torch.mean( interaction_tensor, axis = 2 )
            # I2 = torch.mean( interaction_tensor, axis = 1 )
            # I = torch.cat( [I1, I2], axis = 1 )
            I = self.avg2d( interaction_tensor, interaction_mask )

        elif self.input_layer_params[2] == "lin":
            I1 = torch.permute( interaction_tensor, ( 0, 3, 1, 2 ) )
            I1 = self.interface( I1 ).squeeze( -1 )
            I1 = torch.permute( I1, ( 0, 2, 1 ) )
            
            I2 = torch.permute( interaction_tensor, ( 0, 3, 2, 1 ) )
            I2 = self.interface( I2 ).squeeze( -1 )
            I2 = torch.permute( I2, ( 0, 2, 1 ) )
            
            I = torch.cat( [I1, I2], axis = 1 )
        return I


    def capture_interaction( self, z1, z2, aggregate ):
        """
        Combine the input residue features.
        """
        if aggregate == "add":
            z = z1 + z2

        elif aggregate == "substract":
            z = z1 - z2

        elif aggregate == "mag":
            z = torch.abs( z1 - z2 )

        elif aggregate == "multiply":
            z = z1 * z2

        elif aggregate == "os":
            z = z1.unsqueeze( 2 ) + z2.unsqueeze( 1 )
        
        elif aggregate == "od":
            z = torch.abs( z1.unsqueeze( 2 ) - z2.unsqueeze( 1 ) )

        elif aggregate == "op":
            z = z1.unsqueeze( 2 ) * z2.unsqueeze( 1 )

        elif aggregate == "dot":
            z = torch.sum( z1 * z2, dim = 1, keepdim = True )

        elif aggregate == "cosine":
            from torch.nn.functional import cosine_similarity
            z = cosine_similarity( z1, z2, dim = 1, eps = 1e-6 )

        elif aggregate == "concat":
            z = torch.cat( ( z1, z2 ), dim = 1 )

        return z



    def get_interaction_tensor( self, z1, z2 ):
        # Concatenate along the feature dimension.
        if self.input_layer_params[0] == None:
            z_i1 = z1
            z_i2 = z2
            concat_dim = len( z_i1.shape ) - 1
            interaction_tensor = torch.cat( [z_i1, z_i2], concat_dim )

        else:
            op1, op2 = self.input_layer_params[0].split( "-" )

            if self.input_layer_params[1] == "vanilla":
                # Capture interaction for z1-z2 and concat along feature dim.
                # [N, H, C] [N, W, C] --> [N, H, W, 2C]
                # Results in [N, H, W, 2C] tensor.
                z_i1 = self.capture_interaction( z1, z2, op1 )
                z_i2 = self.capture_interaction( z1, z2, op2 )

            elif self.input_layer_params[1] == "concat":
                # Concatenate along the sequence dim before capturing interaction.
                # [N, H, C] [N, W, C] --> [N, H+W, C]
                # Results in [N, H+W, H+W, 2C] tensor.
                z = torch.cat( [z1, z2], axis = 1 )
                z_i1 = self.capture_interaction( z, z, op1 )
                z_i2 = self.capture_interaction( z, z, op2 )

            concat_dim = len( z_i1.shape ) - 1
            interaction_tensor = torch.cat( [z_i1, z_i2], concat_dim )

        return interaction_tensor



    def vanilla_block( self, o ):
        if self.downsampling_layers != None:
            o = self.downsampling_layers( o )

        if self.upsampling_layers != None:
            o = self.upsampling_layers( o )
        return o


    def residual_block( self, o ):
        for i in range( self.num_blocks ):
            o = self.dropout2( o )
            o_ = o
            o = self.downsampling_layers( o )

            o = self.upsampling_layers( o )
            
            if self.hidden_block[1] == None:
                continue
            elif self.hidden_block[1] == "vanilla":
                o = o + o_
            else:
                if self.norm[1] in ["IN", "BN"]:
                    o = torch.permute( o, ( 0, 2, 1 ) )
                    o_ = torch.permute( o_, ( 0, 2, 1 ) )

                if self.hidden_block[1] == "addactivnorm":
                    o = self.norm_layer( self.activation1( o + o_ ) )
                elif self.hidden_block[1] == "addnorm":
                    o = self.norm_layer( o + o_ )

                if self.norm[1] in ["IN", "BN"]:
                    o = torch.permute( o, ( 0, 2, 1 ) )
                    o_ = torch.permute( o_, ( 0, 2, 1 ) )

        return o


    def output_block( self, o ):
        """
        Output block for the neural network.
            returns logits or the sigmoid output as specified.
        """
        logit = self.contact( o )

        if self.temperature != None:
            logit = logit/self.temperature
        
        if self.apply_activation2:
            o = self.activation2( logit )
            return o.squeeze( -1 )
        else:
            o = logit
            return o.squeeze( -1 )

