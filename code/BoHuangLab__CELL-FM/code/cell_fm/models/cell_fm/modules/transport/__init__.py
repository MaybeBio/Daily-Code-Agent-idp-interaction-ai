from .transport import Transport, ModelType, WeightType, PathType, TimestepSampler, Sampler

def create_transport(
    path_type='Linear',
    prediction="velocity",
    loss_weight=None,
    train_eps=None,
    sample_eps=None,
    timestep_sampler='uniform',
    logit_mean=0.0,
    logit_std=1.0,
):
    """function for creating Transport object
    **Note**: model prediction defaults to velocity
    Args:
    - path_type: type of path to use; default to linear
    - learn_score: set model prediction to score
    - learn_noise: set model prediction to noise
    - velocity_weighted: weight loss by velocity weight
    - likelihood_weighted: weight loss by likelihood weight
    - train_eps: small epsilon for avoiding instability during training
    - sample_eps: small epsilon for avoiding instability during sampling
    - timestep_sampler: how to draw the training timestep; 'uniform' or 'logit_normal' (SD3)
    - logit_mean: mean of the logit-normal draw; SD3 convention, positive favours the noise end
    - logit_std: std of the logit-normal draw
    """

    if prediction == "noise":
        model_type = ModelType.NOISE
    elif prediction == "score":
        model_type = ModelType.SCORE
    else:
        model_type = ModelType.VELOCITY

    if loss_weight == "velocity":
        loss_type = WeightType.VELOCITY
    elif loss_weight == "likelihood":
        loss_type = WeightType.LIKELIHOOD
    else:
        loss_type = WeightType.NONE

    sampler_choice = {
        "uniform": TimestepSampler.UNIFORM,
        "logit_normal": TimestepSampler.LOGIT_NORMAL,
    }

    if timestep_sampler not in sampler_choice:
        raise ValueError(
            f"Unknown timestep_sampler '{timestep_sampler}'; expected one of {sorted(sampler_choice)}"
        )

    timestep_sampler = sampler_choice[timestep_sampler]

    path_choice = {
        "Linear": PathType.LINEAR,
        "GVP": PathType.GVP,
        "VP": PathType.VP,
    }

    path_type = path_choice[path_type]

    if (path_type in [PathType.VP]):
        train_eps = 1e-5 if train_eps is None else train_eps
        sample_eps = 1e-3 if train_eps is None else sample_eps
    elif (path_type in [PathType.GVP, PathType.LINEAR] and model_type != ModelType.VELOCITY):
        train_eps = 1e-3 if train_eps is None else train_eps
        sample_eps = 1e-3 if train_eps is None else sample_eps
    else: # velocity & [GVP, LINEAR] is stable everywhere
        train_eps = 0
        sample_eps = 0
    
    # create flow state
    state = Transport(
        model_type=model_type,
        path_type=path_type,
        loss_type=loss_type,
        train_eps=train_eps,
        sample_eps=sample_eps,
        timestep_sampler=timestep_sampler,
        logit_mean=logit_mean,
        logit_std=logit_std,
    )

    return state