from collections import OrderedDict

import torch
from torch import nn
import torch.nn.functional as F


class Permute( nn.Module ):
    def __init__( self, *dims ):
        super( Permute, self ).__init__()
        self.dims = dims

    def forward( self, x ):
        return torch.permute( x, *self.dims )


def create_projection_layers( proj_type, emb_size, projection_dim, multiplier, bias, activation, device ):
    """
    Create projection layer.
        Linear --> ReLU
    """
    if proj_type == "vanilla":
        print( "Using vanilla projection layer..." )
        return nn.Sequential( 
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation
                            ).to( device )
    
    elif proj_type == "ln1":
        print( "Using projection layer with Lnorm over the inputs..." )
        return nn.Sequential( 
                            nn.LayerNorm( emb_size ),
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation
                            ).to( device )

    elif proj_type == "ln2":
        print( "Using projection layer with Lnorm post activation..." )
        return nn.Sequential( 
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation,
                            nn.LayerNorm( projection_dim )
                            ).to( device )

    elif proj_type == "in1":
        print( "Using projection layer with Inorm1d over the inputs..." )
        return nn.Sequential( 
                            Permute( ( 0, 2, 1 ) ),
                            nn.InstanceNorm1d( emb_size ),
                            Permute( ( 0, 2, 1 ) ),
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation
                            ).to( device )

    elif proj_type == "in2":
        print( "Using projection layer with Inorm1d post activation..." )
        return nn.Sequential( 
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation,
                            # nn.InstanceNorm1d( 100 ),
                            Permute( ( 0, 2, 1 ) ),
                            nn.InstanceNorm1d( projection_dim ),
                            Permute( ( 0, 2, 1 ) )
                            ).to( device )

    elif proj_type == "bn1":
        print( "Using projection layer with Bnorm1d over the inputs..." )
        return nn.Sequential( 
                            Permute( ( 0, 2, 1 ) ),
                            nn.BatchNorm1d( emb_size ),
                            Permute( ( 0, 2, 1 ) ),
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation
                            ).to( device )

    elif proj_type == "bn2":
        print( "Using projection layer with Bnorm1d post activation..." )
        return nn.Sequential( 
                            nn.Linear( in_features = emb_size, out_features = projection_dim, bias = bias ),
                            activation,
                            Permute( ( 0, 2, 1 ) ),
                            nn.BatchNorm1d( projection_dim ),
                            Permute( ( 0, 2, 1 ) )
                            ).to( device )


def create_upsampling_layers( num_upsample_layers, state, scale, bias, activation, us_dropout_prob, norm, device ):
    """
    Create hidden layers, each time doubling the no. of features.
        Dropout --> Linear --> ReLU
    """
    layers = []
    us_layers = ""

    activ, activation = activation
    for N in range( num_upsample_layers ):
        in_features = ( state )*scale**N
        out_features = ( state )*scale**( N + 1 )
        
        if N == 0:
            us_layers += f"{in_features}"
            us_layers += f" --> {out_features}"
        else:
            us_layers += f" --> {out_features}"

        if norm[0]:
            # BNorm/INorm require C as the 1st dimension.
            # [N, H, C] --> [N, C, H]
            if norm[1] in ["BN", "IN"]:
                layers.append(
                  ( f"US_Permute_{N}1", Permute( ( 0, 2, 1 ) ) )
                  )
                if norm[1] == "IN":
                    layers.append(
                      ( f"US_INorm1d{N}" ,nn.InstanceNorm1d( in_features ) )
                      )
                elif norm[1] == "BN":
                    layers.append(
                      ( f"US_BNorm1d{N}" ,nn.BatchNorm1d( in_features ) )
                      )
                # [N, C, H] --> [N, H, C]
                layers.append(
                  ( f"US_Permute_{N}2", Permute( ( 0, 2, 1 ) ) )
                  )
            else:
                layers.append(
                  ( f"US_LNorm{N}" ,nn.LayerNorm( in_features ) )
                  )

        else:
            print( f"Dropout in US layer: {us_dropout_prob}" )
            layers.append(
                ( f"US_Dropout_{N}", nn.Dropout( p = us_dropout_prob ) )
                )


        layers.append(
            ( f"US_Linear_{N}", nn.Linear( in_features = in_features, out_features = out_features, bias = bias ) )
            )

        layers.append(
            ( f"US_{activ}_{N}", activation )
            )

    print( f"Upsampling Layers: {us_layers}" )
    if num_upsample_layers == 0:
        return None
    else:
        return nn.Sequential( OrderedDict( layers ) ).to( device )


def create_downsampling_layers( num_downsample_layers, state, scale, bias, activation, ds_dropout_prob, norm, device ):
    """
    Create hidden layers, each time halving the no. of features.
        Linear --> ReLU
    """
    layers = []
    ds_layers = ""

    activ, activation = activation
    for N in range( num_downsample_layers ):
        in_features = ( state )//scale**N
        out_features = ( state )//scale**( N + 1 )

        if N == 0:
            ds_layers += f"{in_features}"
            ds_layers += f" --> {out_features}"
        else:
            ds_layers += f" --> {out_features}"

        if norm[0]:
            # BNorm/INorm require C as the 1st dimension.
            # [N, H, C] --> [N, C, H]
            if norm[1] in ["BN", "IN"]:
                layers.append(
                  ( f"DS_Permute_{N}1", Permute( ( 0, 2, 1 ) ) )
                  )
                if norm[1] == "IN":
                    layers.append(
                      ( f"US_INorm1d{N}" ,nn.InstanceNorm1d( in_features ) )
                      )
                elif norm[1] == "BN":
                    layers.append(
                      ( f"DS_BNorm1d{N}" ,nn.BatchNorm1d( in_features ) )
                      )
                # [N, C, H] --> [N, H, C]
                layers.append(
                  ( f"DS_Permute_{N}2", Permute( ( 0, 2, 1 ) ) )
                  )
            else:
                layers.append(
                	( f"DS_LNorm{N}" ,nn.LayerNorm( in_features ) )
                	)
        else:
            print( f"Dropout in DS layer: {ds_dropout_prob}" )
            layers.append(
                ( f"DS_Dropout_{N}", nn.Dropout( p = ds_dropout_prob ) )
                )

        layers.append(
            ( f"DS_Linear_{N}", nn.Linear( in_features = in_features, out_features = out_features, bias = bias ) )
            )

        layers.append(
            ( f"DS_{activ}_{N}", activation )
            )

    print( f"Downsampling Layers: {ds_layers}" )
    if num_downsample_layers == 0:
        return None
    else:
        return nn.Sequential( OrderedDict( layers ) ).to( device )
