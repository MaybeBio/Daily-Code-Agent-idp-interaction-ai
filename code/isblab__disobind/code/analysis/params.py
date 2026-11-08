####### Store parameter files anme for each kinda model. #######
######## ------>"May the Force serve u well..." <------#########
################################################################
from typing import Dict

def parameter_files( model_version: int ) -> Dict:
	params = {
		"interaction": get_interaction_params_dict( model_version ),
		"interface": get_inteface_params_dict( model_version )
	}
	return params


def get_interaction_params_dict( model_version: int ) -> Dict:
	if model_version == 19:
		print( "Using v_19 model parameters for interaction prediction..." )
		interaction_param = {
			# v_19 params -------------
			"cg_1": [
				"Epsilon_3_6.2",
				"model_global-[256, 'ln2', True, 1, '']_[0, 3, 0, 2, 'vanilla', '']_['elu', None]_0.0004_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.98__0"
			],
			"cg_5": [
				"Epsilon_3_6.1",
				"model_global-[128, 'ln2', True, 1, '']_[0, 3, 0, 2, 'vanilla', '']_['elu', None]_0.0002_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.97__0"
			],
			"cg_10": [
				"Epsilon_3_6",
				"model_global-[128, 'ln2', True, 1, '']_[0, 3, 0, 2, 'vanilla', '']_['elu', None]_0.0001_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.97__0"
			],
		}
	else:
		raise ValueError( f"Unsupported model version = {model_version}. (Supported -> 19)..." )

	return interaction_param


def get_inteface_params_dict( model_version: int ) -> Dict:
	if model_version == 19:
		print( "Using v_19 model parameters for interface prediction..." )
		interface_param = {
			# v_19 params -------------
			"cg_1": [
				"Epsilon_3_16",
				"model_global-[128, 'ln2', True, 1, '']_[0, 0, 0, 0, 'vanilla', '']_['elu', None]_0.0002_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.98__0"
			],
			"cg_5": [
				"Epsilon_3_16.1",
				"model_global-[128, 'ln2', True, 1, '']_[0, 0, 0, 0, 'vanilla', '']_['elu', None]_0.0002_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.91__0"
			],
			"cg_10": [
				"Epsilon_3_16.2",
				"model_global-[128, 'ln2', True, 1, '']_[0, 0, 0, 0, 'vanilla', '']_['elu', None]_0.0002_0.5_[0.9, 3]_[0.2, 0, 0, 0, 0]_0.05_0.87__0"
			]
		}
	else:
		raise ValueError( f"Unsupported model version = {model_version}. (Supported -> 19)..." )

	return interface_param
