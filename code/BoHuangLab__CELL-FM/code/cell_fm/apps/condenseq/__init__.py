"""The CondenSeq titration app: sequence -> condensate probability curve, AUC and AAC.

These two modules used to live only in the (gitignored) huggingface_space/ working copy,
which meant the demo's logic existed solely inside a deployment artefact and was not on
GitHub at all. The Colab notebooks now fetch model code from the repo, so anything they
import has to be tracked -- hence this package. The Space imports the same modules from
here rather than keeping its own copy.
"""

from . import metrics, pipeline  # noqa: F401
