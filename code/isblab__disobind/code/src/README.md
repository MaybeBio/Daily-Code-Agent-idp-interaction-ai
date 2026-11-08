# Model training
For training the model, move to `src` directory.  

Specify the configurations in the `model_versions.py` file and run to create a CONFIG_FILE:  
```
python model_versions.py
```

Next, start model training using:
```
python hparams_search.py -f [CONFIG_FILE] -m manual