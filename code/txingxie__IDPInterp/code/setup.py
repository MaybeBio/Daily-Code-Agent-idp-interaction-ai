from setuptools import setup, find_packages

setup(
    name="idpinterp",
    version="0.1.0",
    description="Interpretability of intrinsically disordered proteins via sparse autoencoders",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "torch>=2.0",
        "transformers>=4.30",
        "hydra-core>=1.3",
        "omegaconf>=2.3",
        "wandb>=0.15",
        "pandas>=2.0",
        "pyarrow>=12.0",
        "numpy>=1.24",
        "requests>=2.28",
        "tqdm>=4.65",
        "scikit-learn>=1.3",
        "matplotlib>=3.7",
    ],
)
