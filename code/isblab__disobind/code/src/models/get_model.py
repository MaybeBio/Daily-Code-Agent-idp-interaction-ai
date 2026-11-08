from src.models.Epsilon_3 import Epsilon_3


def get_model( config ):
	if config.Model == "Epsilon_3":

		model = Epsilon_3(
			emb_size = config.emb_size,
			projection_layer = config.projection_layer,
			input_layer = config.input_layer,
			output_dim = config.output_dim,
			num_samples = config.num_samples,
			num_hid_layers = config.num_hid_layers,
			bias = config.bias,
			dropouts = config.dropouts,
			norm = config.norm,
			max_len = config.max_seq_len,
			activation1 = config.activation1,
			activation2 = config.activation2,
			temperature = config.temperature,
			output_layer = config.output_layer,
			device = config.device,
			objective = config.objective
		)

		return model
